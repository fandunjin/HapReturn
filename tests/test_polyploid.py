import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from hapreturn.cli import main


def alignment(query, target):
    return f"{query}\t100\t0\t100\t+\t{target}\t100\t0\t100\t100\t100\t60\ttp:A:P\n"


def fixture(root, ploidy, n_sources, n_groups=3, first_source_count=None):
    labels = [f"copy{i + 1:02d}" for i in range(n_sources)]
    expected = {label: [] for label in labels}
    groups, records = ["group_id\tseq_id\n"], []
    for group in range(1, n_groups + 1):
        for copy in range(ploidy):
            name = f"unit{group}_copy{copy + 1:02d}"
            source = (0 if copy < first_source_count else 1) if first_source_count is not None else copy % n_sources
            expected[labels[source]].append(name)
            groups.append(f"unit{group}\t{name}\n")
            records.append(f">{name} chromosome\n" + "ACGT" * 25 + "\n")
    (root / "final.fa").write_text("".join(records), encoding="utf-8")
    (root / "groups.tsv").write_text("".join(groups), encoding="utf-8")
    manifest = ["source\tpath\tformat\n"]
    for label in labels:
        (root / f"{label}.paf").write_text("".join(alignment(f"{label}_{name}", name) for name in expected[label]), encoding="utf-8")
        manifest.append(f"{label}\t{label}.paf\tpaf\n")
    (root / "sources.tsv").write_text("".join(manifest), encoding="utf-8")
    return expected


def run_fixture(root, ploidy, extra=()):
    return subprocess.run([sys.executable, "-m", "hapreturn", "--assignment-mode", "provenance", "--final-fasta", str(root / "final.fa"),
                           "--sources", str(root / "sources.tsv"), "--homolog-groups", str(root / "groups.tsv"),
                           "--ploidy", str(ploidy), "--min-aln-len", "1", "--no-plots", "--out-dir", str(root / "out"),
                           *extra], text=True, capture_output=True)


class PolyploidTests(unittest.TestCase):
    def test_same_command_splits_2n_without_ploidy_options(self):
        for ploidy in (2, 4, 8, 12, 16):
            with self.subTest(ploidy=ploidy), tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
                root = Path(directory).resolve()
                fixture(root, ploidy, 2)
                result = subprocess.run([sys.executable, "-m", "hapreturn", "--assignment-mode", "provenance", "--final-fasta", str(root / "final.fa"),
                                         "--hap1-paf", str(root / "copy01.paf"), "--hap2-paf", str(root / "copy02.paf"),
                                         "--min-aln-len", "1", "--no-plots", "--out-dir", str(root / "out")], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                recovered = []
                for label in ("hap1", "hap2"):
                    names = [line[1:] for line in (root / f"out/{label}.fa").read_text().splitlines() if line.startswith(">")]
                    self.assertEqual(len(names), ploidy * 3 // 2)
                    recovered.extend(names)
                self.assertEqual(len(set(recovered)), ploidy * 3)

    def test_more_than_two_sources_rejected(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            root = Path(directory).resolve()
            fixture(root, 8, 8)
            result = run_fixture(root, 8)
            self.assertEqual(result.returncode, 1)
            self.assertIn("Exactly two original sources", result.stderr)

    def test_two_sources_at_high_ploidy(self):
        for ploidy in (8, 12, 16):
            with self.subTest(ploidy=ploidy), tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
                root = Path(directory).resolve()
                fixture(root, ploidy, 2)
                result = subprocess.run([sys.executable, "-m", "hapreturn", "--assignment-mode", "provenance", "--final-fasta", str(root / "final.fa"),
                                         "--hap1-paf", str(root / "copy01.paf"), "--hap2-paf", str(root / "copy02.paf"),
                                         "--homolog-groups", str(root / "groups.tsv"), "--ploidy", str(ploidy),
                                         "--source-copies", f"hap1={ploidy // 2}", "--source-copies", f"hap2={ploidy // 2}",
                                         "--min-aln-len", "1", "--no-plots", "--out-dir", str(root / "out")], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                for label in ("hap1", "hap2"):
                    self.assertEqual((root / f"out/{label}.fa").read_text().count(">"), ploidy // 2 * 3)

    def test_octoploid_unequal_source_copy_counts(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            root = Path(directory).resolve()
            fixture(root, 8, 2, first_source_count=3)
            result = run_fixture(root, 8, ("--source-copies", "copy01=3", "--source-copies", "copy02=5"))
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_wrong_source_copy_counts_are_reported(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            root = Path(directory).resolve()
            fixture(root, 8, 2, first_source_count=3)
            result = run_fixture(root, 8, ("--source-copies", "copy01=4", "--source-copies", "copy02=4"))
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertIn("groups_to_review=3", result.stdout)

    def test_two_sources_tie_is_unresolved(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            root = Path(directory).resolve()
            fixture(root, 8, 2)
            with (root / "copy02.paf").open("a") as handle:
                handle.write(alignment("duplicate", "unit1_copy01"))
            result = run_fixture(root, 8)
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertIn(">unit1_copy01\n", (root / "out/unresolved.fa").read_text())

    def test_manifest_rejects_duplicate_labels(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            root = Path(directory).resolve()
            fixture(root, 8, 8)
            with (root / "sources.tsv").open("a") as handle:
                handle.write("copy01\tcopy01.paf\tpaf\n")
            result = run_fixture(root, 8)
            self.assertEqual(result.returncode, 1)
            self.assertIn("duplicate source label", result.stderr)

    def test_named_fasta_sources_invoke_alignment(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            root = Path(directory).resolve()
            fixture(root, 8, 2)

            def mock_alignment(command, stdout, check):
                label = Path(command[-1]).stem
                stdout.write((root / f"{label}.paf").read_text())

            args = ["hapreturn", "--assignment-mode", "provenance", "--final-fasta", str(root / "final.fa"),
                    "--source-fasta", f"sourceOne={root / 'copy01.fa'}",
                    "--source-fasta", f"sourceTwo={root / 'copy02.fa'}",
                    "--min-aln-len", "1", "--no-plots", "--out-dir", str(root / "out")]
            for label in ("copy01", "copy02"):
                (root / f"{label}.fa").write_text(">original\nAAAA\n")
            with patch.object(sys, "argv", args), patch("hapreturn.cli.subprocess.run", side_effect=mock_alignment) as runner:
                self.assertEqual(main(), 0)
            self.assertEqual(runner.call_count, 2)
            self.assertIn("--secondary=no", runner.call_args.args[0])
            self.assertTrue((root / "out/sourceOne.fa").is_file())


if __name__ == "__main__":
    unittest.main()
