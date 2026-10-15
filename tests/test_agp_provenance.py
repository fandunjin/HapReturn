import csv
import json
import tempfile
import unittest
from pathlib import Path

from hapreturn.agp_cli import prepare_comparison, main
from hapreturn.agp_provenance import attach_provenance, export_partition, assignments_from_paf
from test_agp import agp_text, c


class ProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def result(self, first, second=None):
        for filename,spec in (("before.agp",first),("after.agp",second or first)):
            (self.root/filename).write_text(agp_text(spec),encoding="utf-8")
        return prepare_comparison(self.root/"before.agp",self.root/"after.agp")

    def test_actual_agp_bases_exclude_gaps_and_preserve_mixed_minor_source(self):
        r=self.result({"chr":[c("h1tg000001l",20),("gap",100),c("h2tg000001l",80)]})
        p=attach_provenance(r)
        row=p["chromosomes"][0]
        self.assertEqual((row["hap1_bp"],row["hap2_bp"],row["gap_bp"]),(20,80,100))
        self.assertEqual(row["hap2_fraction"],.8)
        self.assertEqual(row["composition_status"],"mixed")
        self.assertEqual(row["assignment"],"hap2")
        self.assertEqual(len(p["components"]),2)
        self.assertEqual([x["start"] for x in p["components"]],[1,121])

    def test_selected_chromosomes_and_unplaced_do_not_sort_by_length(self):
        r=self.result({"tiny":[c("h1tg000001l",10)],"large":[c("h2tg000001l",200)]})
        p=attach_provenance(r,selected=["tiny"])
        self.assertEqual([x["assignment"] for x in p["chromosomes"]],["hap1","unplaced"])

    def test_original_scaffold_and_source_contig_are_distinct(self):
        r=self.result({"group1":[c("h1tg000001l"),c("h2tg000001l")]},
                      {"new1":[c("h1tg000001l")],"new2":[c("h2tg000001l")]})
        p=attach_provenance(r)
        self.assertEqual({x["original_scaffold"] for x in p["components"]},{"group1"})
        self.assertEqual({x["source_contig"] for x in p["components"]},{"h1tg000001l","h2tg000001l"})

    def test_balanced_assignment_retains_independent_labels_and_flags_changes(self):
        spec={"c1":[c("h1tg000001l",100)],"c2":[c("h1tg000002l",60),c("h2tg000001l",40)]}
        r=self.result(spec)
        p=attach_provenance(r,mode="balanced",groups={"g":["c1","c2"]},ploidy=2)
        row=p["chromosomes"][1]
        self.assertEqual(row["independent_assignment"],"hap1")
        self.assertEqual(row["assignment"],"hap2")
        self.assertIn("constraint_changed_source_majority",row["review_flags"])

    def test_no_equal_copy_constraint_is_invented(self):
        r=self.result({"c1":[c("h1tg000001l")],"c2":[c("h1tg000002l")]})
        p=attach_provenance(r)
        self.assertEqual([x["assignment"] for x in p["chromosomes"]],["hap1","hap1"])
        with self.assertRaisesRegex(ValueError,"explicit homolog groups"):
            attach_provenance(r,mode="balanced",ploidy=2)

    def test_unknown_source_bases_are_not_assigned_by_default(self):
        r=self.result({"c":[c("unknown")]})
        p=attach_provenance(r)
        self.assertEqual(p["chromosomes"][0]["assignment"],"unresolved")
        self.assertEqual(p["review_chromosomes"],1)

    def test_actual_duplicate_placements_count_target_bases_but_flag_qc(self):
        r=self.result({"c":[c("h1tg000001l"),c("h1tg000001l")]})
        p=attach_provenance(r)
        self.assertEqual(p["chromosomes"][0]["hap1_bp"],200)
        self.assertEqual(p["chromosomes"][0]["source_contig_count"],1)
        self.assertIn("duplicate_source_interval",p["chromosomes"][0]["review_flags"])

    def test_agp_and_minimap_assignments_can_differ_without_changing_agp_sources(self):
        r=self.result({"c":[c("h1tg000001l")]})
        mm={"c":{"assignment":"hap2","evidence_assignment":"hap2","hap1_covered_bp":"30","hap2_covered_bp":"100"}}
        p=attach_provenance(r,alignment_rows=mm)
        row=p["method_comparison"][0]
        self.assertEqual(row["comparison_status"],"different_source")
        self.assertEqual(row["agp_assignment"],"hap1")
        self.assertEqual(row["minimap_assignment"],"hap2")
        self.assertEqual(p["components"][0]["source_hap"],"hap1")

    def test_missing_minimap_evidence_is_not_called_agreement(self):
        r=self.result({"c":[c("h1tg000001l")]})
        self.assertEqual(r["agp_provenance"]["method_comparison"][0]["comparison_status"],"not_available")

    def test_minimap_independent_and_group_constrained_labels_both_retained(self):
        r=self.result({"c":[c("h1tg000001l")]})
        p=attach_provenance(r,alignment_rows={"c":{"assignment":"hap1","evidence_assignment":"hap2","assignment_mode":"balanced"}})
        row=p["method_comparison"][0]
        self.assertEqual(row["comparison_status"],"agree")
        self.assertFalse(row["independent_labels_agree"])
        self.assertEqual(row["minimap_assignment_mode"],"balanced")

    def test_paf_union_coverage_is_similarity_not_physical_composition(self):
        r=self.result({"c":[c("h1tg000001l",100)]})
        # Both sources can align across the same target bases, despite pure AGP hap1 provenance.
        for name,length in (("h1",30),("h2",100)):
            (self.root/(name+".paf")).write_text(f"q\t{length}\t0\t{length}\t+\tc\t100\t0\t{length}\t{length}\t{length}\t60\n")
        mm=assignments_from_paf((self.root/"h1.paf",self.root/"h2.paf"),r["after_lengths"],min_aln_len=1)
        self.assertEqual(mm["c"]["hap1_covered_bp"]+mm["c"]["hap2_covered_bp"],130)
        p=attach_provenance(r,alignment_rows=mm)
        self.assertEqual(p["method_comparison"][0]["comparison_status"],"different_source")

    def test_partition_exports_whole_current_sequences_without_renaming(self):
        r=self.result({"c1":[c("h1tg000001l",4)],"c2":[c("h2tg000001l",3)],"other":[c("h1tg000003l",2)]})
        attach_provenance(r,selected=["c1","c2"])
        fasta=self.root/"final.fa"
        fasta.write_text(">c1 metadata\nACGT\n>c2\nTTT\n>other\nAA\n")
        directory=export_partition(r,fasta,self.root)
        self.assertEqual((directory/"hap1.fa").read_text(),">c1\nACGT\n")
        self.assertEqual((directory/"hap2.fa").read_text(),">c2\nTTT\n")
        self.assertEqual((directory/"unplaced.fa").read_text(),">other\nAA\n")

    def test_custom_source_labels_cannot_escape_output_directory(self):
        r=self.result({"c":[c("h1tg000001l")]})
        for labels in (("../evil","hap2"),("hap1","NUL"),("hap1","HAP1")):
            with self.assertRaises(ValueError):
                attach_provenance(r,labels=labels)

    def test_custom_source_palette_has_no_ambiguous_unused_hap_legend(self):
        from hapreturn.agp_report import palette, COLORS
        colors=palette({"agp_provenance":{"labels":["parent_A","parent_B"]},
                        "atoms":[{"source_hap":"parent_A"},{"source_hap":"parent_B"}]})
        self.assertEqual(colors["parent_A"],COLORS["hap1"])
        self.assertEqual(colors["parent_B"],COLORS["hap2"])
        self.assertNotIn("hap1",colors)

    def test_standalone_cli_writes_trace_partition_and_unavailable_comparison(self):
        self.result({"c1":[c("h1tg000001l",4)],"c2":[c("h2tg000001l",3)]})
        (self.root/"final.fa").write_text(">c1\nACGT\n>c2\nTTT\n")
        status=main(["--before-agp",str(self.root/"before.agp"),"--after-agp",str(self.root/"after.agp"),
                     "--chromosome-count","1","--final-fasta",str(self.root/"final.fa"),
                     "--out-dir",str(self.root/"out"),"--no-plots"])
        self.assertEqual(status,0)
        self.assertTrue((self.root/"out/agp_partition/hap1.fa").is_file())
        with (self.root/"out/agp_vs_minimap.tsv").open() as handle:
            self.assertEqual(next(csv.DictReader(handle,delimiter="\t"))["comparison_status"],"not_available")


if __name__=="__main__":
    unittest.main()
