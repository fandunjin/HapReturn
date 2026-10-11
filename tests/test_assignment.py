import csv
import itertools
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from hapreturn.assignment import balanced_assignment, consecutive_groups
from hapreturn.windows import expand_unresolved


class BalancedTests(unittest.TestCase):
    def test_exact_objective_matches_exhaustive_search(self):
        rng = random.Random(19)
        for ploidy in (2, 4, 8, 12):
            for _ in range(12):
                members = [f"c{i}" for i in range(ploidy)]
                support = {n: {"hap1": rng.randrange(100), "hap2": rng.randrange(100)} for n in members}
                assigned, audit = balanced_assignment({"G": members}, support, ["hap1", "hap2"], ploidy)
                scores = [sum(support[n]["hap1" if n in selected else "hap2"] for n in members)
                          for selected in map(set, itertools.combinations(members, ploidy // 2))]
                self.assertEqual(audit["G"]["total_support_bp"], max(scores))
                self.assertEqual(list(assigned.values()).count("hap1"), ploidy // 2)

    def test_three_to_one_changes_smallest_support_penalty(self):
        support = {"a": {"hap1": 100, "hap2": 0}, "b": {"hap1": 90, "hap2": 0},
                   "c": {"hap1": 80, "hap2": 70}, "d": {"hap1": 0, "hap2": 100}}
        result, audit = balanced_assignment({"G": list(support)}, support, ["hap1", "hap2"], 4)
        self.assertEqual(result, {"a": "hap1", "b": "hap1", "c": "hap2", "d": "hap2"})
        self.assertEqual(audit["G"]["total_support_bp"], 360)

    def test_no_evidence_tie_is_audited_and_stable(self):
        records = {f"c{i}": 10 for i in range(8)}
        support = {n: {"hap1": 0, "hap2": 0} for n in records}
        result, audit = balanced_assignment(consecutive_groups(records, 4), support, ["hap1", "hap2"], 4)
        self.assertEqual([n for n in result if result[n] == "hap1"], ["c0", "c1", "c4", "c5"])
        self.assertTrue(all(row["tied_optimum"] for row in audit.values()))

    def test_invalid_group_sizes_are_rejected(self):
        for ploidy in (None, 1, 3):
            with self.assertRaisesRegex(ValueError, "even --ploidy"):
                consecutive_groups({"a": 1, "b": 1}, ploidy)
        with self.assertRaisesRegex(ValueError, "divisible"):
            consecutive_groups({"a": 1, "b": 1}, 4)

    def test_gray_context_expands_and_preserves_raw_evidence(self):
        rows = [("c", 0, 10, 4, 4, "unresolved"), ("c", 10, 20, 2, 0, "unresolved"),
                ("c", 20, 30, 0, 10, "hap2")]
        revised, notes, attempts = expand_unresolved(rows, ["hap1", "hap2"], .25, .6)
        self.assertEqual([r[-1] for r in revised], ["hap1", "hap1", "hap2"])
        self.assertEqual(revised[0][3:5], (4, 4))
        self.assertEqual(notes[("c", 0)], ("unresolved", "expanded_window", 0, 20))
        self.assertTrue(all(row[4] <= 20 for row in attempts))

    def test_empty_gray_run_stays_gray_and_reaches_full_span(self):
        rows = [("c", i, i + 10, 0, 0, "unresolved") for i in range(0, 50, 10)]
        revised, notes, attempts = expand_unresolved(rows, ["hap1", "hap2"], .25, .6)
        self.assertTrue(all(r[-1] == "unresolved" for r in revised))
        self.assertEqual(max(r[5] for r in attempts), 50)

    def test_cli_balances_and_preserves_evidence_assignment(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            root = Path(directory)
            (root / "final.fa").write_text("".join(f">c{i}\n" + "A" * 100 + "\n" for i in range(4)))
            def paf(query, target, span):
                return f"{query}\t100\t0\t{span}\t+\t{target}\t100\t0\t{span}\t{span}\t{span}\t60\n"
            (root / "h1.paf").write_text(paf("a", "c0", 100) + paf("b", "c1", 90) + paf("c", "c2", 80))
            (root / "h2.paf").write_text(paf("d", "c2", 70) + paf("e", "c3", 100))
            result = subprocess.run([sys.executable, "-m", "hapreturn", "--final-fasta", str(root / "final.fa"),
                                     "--hap1-paf", str(root / "h1.paf"), "--hap2-paf", str(root / "h2.paf"),
                                     "--ploidy", "4", "--min-aln-len", "1", "--no-plots", "--out-dir", str(root / "out")], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2, result.stderr)
            with (root / "out/assignment.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([r["assignment"] for r in rows], ["hap1", "hap1", "hap2", "hap2"])
            self.assertEqual(rows[2]["evidence_assignment"], "unresolved")
            self.assertEqual((root / "out/hap1.fa").read_text().count(">"), 2)
            self.assertTrue((root / "out/hap1_to_final.paf").exists())


if __name__ == "__main__":
    unittest.main()
