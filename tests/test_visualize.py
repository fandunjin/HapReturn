import tempfile
import unittest
import subprocess
import sys
import re
from dataclasses import replace
from pathlib import Path

from hapreturn.visualize import make_report, window_support, prepare_data, alignment_matrix
from hapreturn.evidence import Alignment, read_paf
from hapreturn.comparison import ribbon_vertices, common_scale, stacked_positions


class VisualizationTests(unittest.TestCase):
    def test_cli_retains_report_on_review_exit(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            root = Path(directory)
            (root / "final.fa").write_text(">chr1\nAAAAAAAAAA\n")
            (root / "empty.paf").write_text("")
            (root / "out").mkdir()
            (root / "out" / "overview.pdf").write_bytes(b"legacy generated figure")
            result = subprocess.run([sys.executable, "-m", "hapreturn",
                "--final-fasta", str(root / "final.fa"),
                "--hap1-paf", str(root / "empty.paf"), "--hap2-paf", str(root / "empty.paf"),
                "--force", "--out-dir", str(root / "out")], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertTrue((root / "out" / "report.html").is_file())
            self.assertFalse((root / "out" / "overview.pdf").exists())

    def test_window_union_and_partial_end(self):
        self.assertEqual(window_support(25, [(0, 12), (5, 18), (20, 25)], 10), [10, 8, 5])

    def test_report_exports_and_escaping(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            out = Path(directory)
            name = 'chr<&"'
            records = {name: 25, "mixed": 20, "missing": 10}
            intervals = {"hap1": {name: [(0, 25)], "mixed": [(0, 10)]},
                         "hap2": {"mixed": [(10, 20)]}}
            assigned = {name: "hap1", "mixed": "unresolved", "missing": "unresolved"}
            make_report(out, ["hap1", "hap2"], records, intervals, assigned,
                        dict(records, unselected=5), window_size=10)
            self.assertEqual({p.name for p in out.glob("*.pdf")},
                             {"partition.pdf", "alignment_matrix.pdf", "hap_comparison.pdf"})
            for path in out.glob("*.pdf"):
                self.assertTrue(path.read_bytes().startswith(b"%PDF"))
            self.assertEqual(len(re.findall(rb"/Type /Page\b", (out / "partition.pdf").read_bytes())), 2)
            self.assertFalse(list(out.glob("*.png")) + list(out.glob("*.svg")))
            report = (out / "report.html").read_text(encoding="utf-8")
            self.assertIn('chr&lt;&amp;&quot;', report)
            self.assertIn('addEventListener', report)
            self.assertNotIn('<script src=', report)
            rows = (out / "window_provenance.tsv").read_text(encoding="utf-8")
            self.assertIn('mixed\t0\t10\t10\t0\thap1', rows)
            self.assertIn('mixed\t10\t20\t0\t10\thap2', rows)
            self.assertIn('missing\t0\t10\t0\t0\tunresolved', rows)
            ids = re.findall(r'\bid="([^"]+)"', report)
            self.assertEqual(len(ids), len(set(ids)), "Embedded SVG identifiers must be unique")

    def test_matrix_union_is_target_normalized_and_sources_are_distinct(self):
        records = {"chr": 20}
        intervals = {"hap1": {"chr": [(0, 15)]}, "hap2": {"chr": [(0, 5)]}}
        data = prepare_data(["hap1", "hap2"], records, intervals, {"chr": "hap1"}, records)
        first = Alignment("hap1", "same_id", 20, 0, 10, "+", "chr", 20, 0, 10)
        hits = [first, replace(first, start=5, end=15), replace(first, source="hap2", end=5)]
        matrix = alignment_matrix(data, hits, intervals)
        self.assertEqual(matrix[("hap1", "same_id")]["chr"], 0.75)
        self.assertEqual(matrix[("hap2", "same_id")]["chr"], 0.25)

    def test_reverse_strand_ribbon_endpoints_and_shared_scale(self):
        hit = Alignment("hap1", "q", 100, 10, 30, "+", "chr", 200, 40, 80)
        self.assertEqual(ribbon_vertices(hit, (0, 1), (0, 1)), (0.1, 0.3, 0.2, 0.4))
        self.assertEqual(ribbon_vertices(replace(hit, strand="-"), (0, 1), (0, 1)), (0.1, 0.3, 0.4, 0.2))
        scale = common_scale({"q": 100, "tiny": 1}, {"chr": 200})
        left, right = stacked_positions({"q": 100, "tiny": 1}, scale), stacked_positions({"chr": 200}, scale)
        self.assertAlmostEqual(right["chr"][1] / left["q"][1], 2)
        self.assertAlmostEqual(left["q"][1] / left["tiny"][1], 100)

    def test_paf_rejects_bad_query_coordinates(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            path = Path(directory) / "bad.paf"
            path.write_text("q\t10\t0\t11\t+\tchr\t20\t0\t10\t10\t10\t60\n")
            with self.assertRaisesRegex(ValueError, "query coordinates"):
                list(read_paf(path, {"chr": 20}, 0, 1, "hap1"))

    def test_standalone_chromosome_command(self):
        from hapreturn.comparison import main
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            root = Path(directory)
            (root / "final.fa").write_text(">chr\nAAAAAAAAAA\n")
            (root / "one.paf").write_text("q\t10\t0\t10\t-\tchr\t10\t0\t10\t10\t10\t60\n")
            (root / "two.paf").write_text("")
            (root / "assignment.tsv").write_text("seq_id\tassignment\nchr\thap1\n")
            code = main(["--final-fasta", str(root / "final.fa"), "--hap1-paf", str(root / "one.paf"),
                         "--hap2-paf", str(root / "two.paf"), "--assignment", str(root / "assignment.tsv"),
                         "--min-aln-len", "1", "--out-pdf", str(root / "comparison.pdf")])
            self.assertEqual(code, 0)
            self.assertTrue((root / "comparison.pdf").is_file())


if __name__ == "__main__":
    unittest.main()
