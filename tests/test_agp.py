import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from hapreturn.agp import compare_agp, read_agp, read_source_map, write_tables
from hapreturn.agp_cli import main, prepare_comparison


def agp_text(spec):
    """Spec preserves object/record order. Components use 1-based inclusive spans."""
    lines = []
    for obj, items in spec.items():
        cursor = 1
        for part, item in enumerate(items, 1):
            if item[0] == "gap":
                length = item[1]
                fields = [obj, cursor, cursor + length - 1, part, "N", length, "scaffold", "yes", "proximity_ligation"]
            else:
                contig, start, end, strand = item
                length = end - start + 1
                fields = [obj, cursor, cursor + length - 1, part, "W", contig, start, end, strand]
            lines.append("\t".join(map(str, fields)))
            cursor += length
    return "\n".join(lines) + "\n"


def c(name, size=100, strand="+"):
    return name, 1, size, strand


class AgpTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def compare(self, first, second):
        (self.root / "before.agp").write_text(agp_text(first), encoding="utf-8")
        (self.root / "after.agp").write_text(agp_text(second), encoding="utf-8")
        return prepare_comparison(self.root / "before.agp", self.root / "after.agp")

    def test_renaming_and_coordinate_shift_are_not_moves(self):
        r = self.compare({"old": [c("A"), c("B"), c("C")]},
                         {"new": [c("A"), c("X", 25), c("B"), c("C")]})
        self.assertEqual(r["summary"]["changed_unique_bp"], 0)
        self.assertEqual(r["summary"]["after_only_bp"], 25)
        self.assertFalse(any(a["flags"] for a in r["atoms"] if a["source_contig"] != "X"))
        self.assertEqual(r["summary"]["junction_added"], 2)
        self.assertEqual(r["summary"]["junction_removed"], 1)

    def test_whole_reverse_preserves_oriented_junctions(self):
        r = self.compare({"old": [c("A"), ("gap", 100), c("B", strand="-")]},
                         {"new": [c("B"), ("gap", 100), c("A", strand="-")]})
        s = r["summary"]
        self.assertEqual(s["whole_scaffold_reverse"], 1)
        self.assertEqual(s["equivalent_scaffold_rename"], 1)
        self.assertEqual(s["changed_unique_bp"], 0)
        self.assertEqual(len(r["junction_changes"]), 0)

    def test_internal_inversion_keeps_block_internal_adjacency(self):
        r = self.compare({"chr": [c("A"), c("B"), c("C"), c("D")]},
                         {"chr": [c("A"), c("C", strand="-"), c("B", strand="-"), c("D")]})
        s = r["summary"]
        self.assertEqual(s["whole_scaffold_reverse"], 0)
        self.assertEqual(s["orientation_changed_bp"], 200)
        self.assertEqual(s["junction_added"], 2)
        self.assertEqual(s["junction_removed"], 2)
        self.assertEqual(s["relative_order_changed_bp"], 200)

    def test_cross_scaffold_transfer_uses_overlap_anchor_not_names(self):
        r = self.compare({"b1": [c("A"), c("B", 10), c("C")], "b2": [c("D"), c("E")]},
                         {"a1": [c("A"), c("C")], "a2": [c("D"), c("B", 10), c("E")]})
        self.assertEqual(r["summary"]["cross_scaffold_bp"], 10)
        self.assertEqual(r["summary"]["junction_removed"], 3)
        self.assertEqual(r["summary"]["junction_added"], 3)
        self.assertEqual(r["summary"]["complex_scaffold_reorganization"], 1)

    def test_component_split_maps_atomic_intervals_exactly_on_minus(self):
        r = self.compare({"b": [c("A", strand="-")]},
                         {"x": [("A", 1, 40, "-")], "y": [("A", 41, 100, "-")]})
        self.assertEqual(r["summary"]["unique_bp"], 100)
        first, second = r["atoms"]
        self.assertEqual((first["before"][0]["start"], first["before"][0]["end"]), (61, 100))
        self.assertEqual((second["before"][0]["start"], second["before"][0]["end"]), (1, 60))
        self.assertEqual(r["summary"]["scaffold_split"], 1)
        self.assertEqual(r["boundary_changes"][0]["boundary_after_bp"], 40)
        self.assertEqual(r["summary"]["junction_removed"], 1)

    def test_component_split_without_break_retains_continuity(self):
        r = self.compare({"b": [c("A")]}, {"a": [("A", 1, 40, "+"), ("A", 41, 100, "+")]})
        self.assertEqual(r["summary"]["changed_unique_bp"], 0)
        self.assertEqual(len(r["junction_changes"]), 0)
        self.assertEqual(len(r["boundary_changes"]), 1)

    def test_duplicates_excluded_and_all_occurrences_retained(self):
        r = self.compare({"b": [c("A")]}, {"x": [c("A")], "y": [c("A")]})
        self.assertEqual(r["summary"]["ambiguous_bp"], 100)
        self.assertEqual(r["summary"]["after_component_bp"], 200)
        self.assertEqual(r["summary"]["source_union_bp"], 100)
        self.assertEqual(len(r["atoms"][0]["after"]), 2)
        self.assertEqual(r["summary"]["scaffold_split"], 0)
        write_tables(r, self.root / "out")
        with (self.root / "out/segment_correspondence.tsv").open(encoding="utf-8") as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(row["after_count"], "2")
        self.assertEqual(row["after_object"], "")
        self.assertEqual(len(json.loads(row["after_occurrences"])), 2)

    def test_unknown_direction_is_not_called_inversion(self):
        r = self.compare({"b": [c("A", strand="?"), c("B")]}, {"a": [c("A", strand="-"), c("B")]})
        self.assertEqual(r["summary"]["orientation_changed_bp"], 0)
        self.assertGreater(r["summary"]["review_items"], 0)
        self.assertEqual(r["atoms"][0]["before"][0]["coordinate_status"], "orientation_assumed_plus")

    def test_gap_length_change_is_separate_from_structural_movement(self):
        r = self.compare({"b": [c("A"), ("gap", 100), c("B")]},
                         {"a": [c("A"), ("gap", 150), c("B")]})
        self.assertEqual(r["summary"]["gap_changed"], 1)
        self.assertEqual(r["summary"]["junction_added"], 0)
        self.assertEqual(r["summary"]["changed_unique_bp"], 0)

    def test_scaffold_order_is_separate_from_internal_order(self):
        r = self.compare({"b1": [c("A")], "b2": [c("B")], "b3": [c("C")]},
                         {"a1": [c("C")], "a2": [c("B")], "a3": [c("A")]})
        self.assertEqual(r["summary"]["scaffold_order_changed"], 3)
        self.assertEqual(r["summary"]["changed_unique_bp"], 0)

    def test_equal_overlap_tie_does_not_choose_arbitrary_anchor(self):
        r = self.compare({"b": [c("A"), c("B")]}, {"x": [c("A")], "y": [c("B")]})
        self.assertFalse(any(x["reciprocal_anchor"] for x in r["scaffold_map"]))
        self.assertEqual(r["summary"]["cross_scaffold_bp"], 0)
        self.assertEqual(r["summary"]["unanchored_scaffold_bp"], 200)

    def test_missing_interval_is_not_silently_a_deletion(self):
        r = self.compare({"b": [c("A"), c("B")]}, {"a": [c("A")]})
        self.assertEqual(r["summary"]["before_only_bp"], 100)
        self.assertIn("before_only", {e["type"] for e in r["events"]})

    def test_invalid_tiling_length_or_orientation_has_line_diagnostic(self):
        for bad in ("c\t2\t10\t1\tW\tA\t1\t9\t+\n", "c\t1\t10\t1\tW\tA\t1\t9\t+\n",
                    "c\t1\t10\t1\tW\tA\t1\t10\tbad\n", "c\t1\t10\t1\tU\t10\tscaffold\tyes\tproximity_ligation\n"):
            path = self.root / "bad.agp"
            path.write_text(bad)
            with self.assertRaisesRegex(ValueError, r"bad.agp:1:"):
                read_agp(path)

    def test_source_map_overrides_name_inference(self):
        self.compare({"b": [c("lcl|h1tg000001l")]}, {"a": [c("lcl|h1tg000001l")]})
        path = self.root / "sources.tsv"
        path.write_text("contig_id\tsource\nlcl|h1tg000001l\tparent_A\n", encoding="utf-8")
        r = prepare_comparison(self.root / "before.agp", self.root / "after.agp", source_map=path)
        self.assertEqual(r["atoms"][0]["source_hap"], "parent_A")
        self.assertEqual(r["atoms"][0]["source_method"], "source_map")
        path.write_text("contig_id\tsource\n", encoding="utf-8")
        r = prepare_comparison(self.root / "before.agp", self.root / "after.agp", source_map=path)
        self.assertEqual(r["atoms"][0]["source_hap"], "unresolved")

    def test_cli_review_exit_and_input_overwrite_guard(self):
        self.compare({"b": [c("A")]}, {"a": [c("A", strand="?")]})
        args = ["--before-agp", str(self.root / "before.agp"), "--after-agp", str(self.root / "after.agp"),
                "--out-dir", str(self.root / "out"), "--no-plots"]
        self.assertEqual(main(args), 2)
        self.assertEqual(main(args + ["--force", "--allow-review"]), 0)
        # An arbitrary AGP input named like an output must never be overwritten.
        input_path = self.root / "curation_events.tsv"
        input_path.write_text(agp_text({"b": [c("A")]}))
        saved = input_path.read_bytes()
        self.assertEqual(main(["--before-agp", str(input_path), "--after-agp", str(self.root / "after.agp"),
                               "--out-dir", str(self.root), "--no-plots", "--force"]), 1)
        self.assertEqual(saved, input_path.read_bytes())

    def test_plot_cli_portable_report_escapes_identifiers(self):
        strange = '</script><img src=x onerror="alert(1)">'
        self.compare({"old": [c(strange)]}, {"new": [c(strange)]})
        self.assertEqual(main(["--before-agp", str(self.root / "before.agp"), "--after-agp", str(self.root / "after.agp"),
                               "--out-dir", str(self.root / "out"), "--allow-review"]), 0)
        content = (self.root / "out/curation_report.html").read_text(encoding="utf-8")
        self.assertNotIn(strange, content)
        self.assertIn("\\u003c/script>", content)
        self.assertTrue((self.root / "out/curation_comparison.pdf").read_bytes().startswith(b"%PDF"))

    def test_main_pipeline_agp_integration_and_refresh(self):
        # The tiny existing demo has two 100-bp chromosomes; retain real PAF QC.
        spec = {"chr1": [c("h1tg000001l")], "chr2": [c("h2tg000001l")]}
        self.compare(spec, spec)
        fasta = self.root / "final.fa"
        fasta.write_text(">chr1\n" + "A" * 100 + "\n>chr2\n" + "A" * 100 + "\n")
        for label, target in (("hap1", "chr1"), ("hap2", "chr2")):
            (self.root / (label + ".paf")).write_text(f"q\t100\t0\t100\t+\t{target}\t100\t0\t100\t100\t100\t60\n")
        command = [sys.executable, "-m", "hapreturn", "--final-fasta", str(fasta),
                   "--hap1-paf", str(self.root / "hap1.paf"), "--hap2-paf", str(self.root / "hap2.paf"),
                   "--ploidy", "2", "--min-aln-len", "1", "--before-agp", str(self.root / "before.agp"),
                   "--after-agp", str(self.root / "after.agp"), "--out-dir", str(self.root / "result")]
        completed = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn('Joint evidence matrix', (self.root / "result/report.html").read_text(encoding="utf-8"))
        saved = json.loads((self.root / "result/curation/agp_comparison.json").read_text(encoding="utf-8"))
        self.assertEqual({x["after_assignment"] for x in saved["atoms"]}, {"hap1", "hap2"})
        refreshed = subprocess.run([sys.executable, "-m", "hapreturn.report", "--result-dir", str(self.root / "result"),
                                    "--out-dir", str(self.root / "refreshed"), "--min-aln-len", "1"], capture_output=True, text=True)
        self.assertEqual(refreshed.returncode, 0, refreshed.stderr)
        self.assertTrue((self.root / "refreshed/curation/curation_report.html").is_file())


if __name__ == "__main__":
    unittest.main()
