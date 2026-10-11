import tempfile
import unittest
import subprocess
import sys
from pathlib import Path

from hapreturn.visualize import make_report, window_support


class VisualizationTests(unittest.TestCase):
    def test_cli_retains_report_on_review_exit(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            root = Path(directory)
            (root / "final.fa").write_text(">chr1\nAAAAAAAAAA\n")
            (root / "empty.paf").write_text("")
            result = subprocess.run([sys.executable, "-m", "hapreturn",
                "--final-fasta", str(root / "final.fa"),
                "--hap1-paf", str(root / "empty.paf"), "--hap2-paf", str(root / "empty.paf"),
                "--out-dir", str(root / "out")], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertTrue((root / "out" / "report.html").is_file())

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
            for plot in ("coverage", "assignment", "tracks", "summary", "overview"):
                for ext in ("png", "svg", "pdf"):
                    self.assertGreater((out / f"{plot}.{ext}").stat().st_size, 100)
            report = (out / "report.html").read_text(encoding="utf-8")
            self.assertIn('chr&lt;&amp;&quot;', report)
            self.assertIn('addEventListener', report)
            self.assertNotIn('<script src=', report)
            rows = (out / "window_provenance.tsv").read_text(encoding="utf-8")
            self.assertIn('mixed\t0\t10\t10\t0\thap1', rows)
            self.assertIn('mixed\t10\t20\t0\t10\thap2', rows)
            self.assertIn('missing\t0\t10\t0\t0\tunresolved', rows)


if __name__ == "__main__":
    unittest.main()
