import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from hapreturn.cli import load_fasta, paf_intervals, union_length


SCRIPT = Path(__file__).resolve()


def paf(query, target, target_len, start, end, secondary=False):
    span = end - start
    tag = "\ttp:A:S" if secondary else "\ttp:A:P"
    return f"{query}\t{span}\t0\t{span}\t+\t{target}\t{target_len}\t{start}\t{end}\t{span}\t{span}\t60{tag}\n"


class SplitTests(unittest.TestCase):
    def test_overlapping_intervals_count_once(self):
        self.assertEqual(union_length([(0, 70), (20, 80), (90, 100)]), 90)

    def test_empty_header_is_clear_input_error(self):
        with tempfile.TemporaryDirectory(dir=SCRIPT.parent) as directory:
            path = Path(directory) / "bad.fa"
            path.write_text(">\nAAAA\n")
            with self.assertRaisesRegex(ValueError, "Empty FASTA identifier"):
                load_fasta(path)

    def test_mismatched_paf_target_length_rejected(self):
        with tempfile.TemporaryDirectory(dir=SCRIPT.parent) as directory:
            path = Path(directory) / "bad.paf"
            path.write_text(paf("query", "chr1", 101, 0, 90))
            with self.assertRaisesRegex(ValueError, "coordinates/length disagree"):
                paf_intervals(path, {"chr1": "A" * 100}, 10, 1)

    def test_tetraploid_and_unplaced_with_group_qc(self):
        with tempfile.TemporaryDirectory(dir=SCRIPT.parent) as directory:
            root = Path(directory)
            fasta = root / "final.fa"
            fasta.write_text("".join(f">chr{i}\n" + "A" * 100 + "\n" for i in range(1, 5)) + ">unplaced\nTT\n")
            (root / "chromosomes.txt").write_text("chr1\nchr2\nchr3\nchr4\n")
            (root / "groups.tsv").write_text("group_id\tseq_id\nC1\tchr1\nC1\tchr2\nC1\tchr3\nC1\tchr4\n")
            (root / "h1.paf").write_text(paf("h1a", "chr1", 100, 0, 90) + paf("h1b", "chr3", 100, 0, 90)
                                         + paf("repeat", "chr2", 100, 0, 90, secondary=True))
            (root / "h2.paf").write_text(paf("h2a", "chr2", 100, 0, 90) + paf("h2b", "chr4", 100, 0, 90))
            result = subprocess.run([sys.executable, "-m", "hapreturn", "--final-fasta", str(fasta),
                                     "--hap1-paf", str(root / "h1.paf"), "--hap2-paf", str(root / "h2.paf"),
                                     "--chromosome-list", str(root / "chromosomes.txt"), "--ploidy", "4",
                                     "--homolog-groups", str(root / "groups.tsv"), "--hap1-per-group", "2",
                                     "--min-aln-len", "1", "--out-dir", str(root / "out")], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn(">chr1\n", (root / "out/hap1.fa").read_text())
            self.assertIn(">chr3\n", (root / "out/hap1.fa").read_text())
            self.assertIn(">chr2\n", (root / "out/hap2.fa").read_text())
            self.assertIn(">unplaced\n", (root / "out/unplaced.fa").read_text())
            with (root / "out/group_qc.tsv").open() as handle:
                self.assertEqual(next(csv.DictReader(handle, delimiter="\t"))["status"], "pass")

    def test_ambiguous_fails_closed(self):
        with tempfile.TemporaryDirectory(dir=SCRIPT.parent) as directory:
            root = Path(directory)
            (root / "final.fa").write_text(">chr1\n" + "A" * 100 + "\n")
            for label in ("h1", "h2"):
                (root / f"{label}.paf").write_text(paf(label, "chr1", 100, 0, 80))
            result = subprocess.run([sys.executable, "-m", "hapreturn", "--final-fasta", str(root / "final.fa"),
                                     "--hap1-paf", str(root / "h1.paf"), "--hap2-paf", str(root / "h2.paf"),
                                     "--min-aln-len", "1", "--out-dir", str(root / "out")], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn(">chr1\n", (root / "out/unresolved.fa").read_text())

    def test_triploid_unequal_sources(self):
        with tempfile.TemporaryDirectory(dir=SCRIPT.parent) as directory:
            root = Path(directory)
            (root / "final.fa").write_text("".join(f">c{i}\n" + "A" * 100 + "\n" for i in range(1, 4)))
            (root / "groups.tsv").write_text("group_id\tseq_id\nG1\tc1\nG1\tc2\nG1\tc3\n")
            (root / "h1.paf").write_text(paf("a", "c1", 100, 0, 90))
            (root / "h2.paf").write_text(paf("b", "c2", 100, 0, 90) + paf("c", "c3", 100, 0, 90))
            result = subprocess.run([sys.executable, "-m", "hapreturn", "--final-fasta", str(root / "final.fa"),
                                     "--hap1-paf", str(root / "h1.paf"), "--hap2-paf", str(root / "h2.paf"),
                                     "--homolog-groups", str(root / "groups.tsv"), "--ploidy", "3",
                                     "--hap1-per-group", "1", "--min-aln-len", "1", "--out-dir", str(root / "out")],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((root / "out/hap1.fa").read_text().count(">"), 1)
            self.assertEqual((root / "out/hap2.fa").read_text().count(">"), 2)


if __name__ == "__main__":
    unittest.main()
