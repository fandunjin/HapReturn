import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from hapreturn.hybrid import fuse_evidence

LABELS = ['hap1', 'hap2']


def atom(start, end, source='hap1', status='unique', strand='+'):
    return {'source_hap': source, 'source_contig': 'contig', 'mapping_status': status,
            'before': [{'strand': '+'}], 'after': [{'object': 'chr', 'start': start + 1, 'end': end, 'strand': strand}]}


def fuse(atoms, first=(), second=()):
    return fuse_evidence({'atoms': atoms}, {'hap1': {'chr': first}, 'hap2': {'chr': second}}, ['chr'], LABELS)


class HybridTests(unittest.TestCase):
    def test_direct_short_fragment_survives_no_alignment(self):
        support, audit, _ = fuse([atom(0, 4355)])
        self.assertEqual(support['chr'], {'hap1': 4355, 'hap2': 0})
        self.assertEqual(audit['chr']['agp_direct_bp'], 4355)

    def test_overlapping_alignment_never_doubles_direct_support(self):
        support, audit, _ = fuse([atom(0, 100)], [(0, 100), (20, 80)], [(20, 80)])
        self.assertEqual(sum(support['chr'].values()), 100)
        self.assertEqual(audit['chr']['conflict_bp'], 0)

    def test_unknown_fill_excludes_competition_and_gap(self):
        support, audit, segments = fuse([atom(0, 100, 'unresolved'), atom(200, 300, 'unresolved')],
                                       [(0, 250), (0, 50)], [(50, 300)])
        self.assertEqual(support['chr'], {'hap1': 50, 'hap2': 50})
        self.assertEqual(audit['chr']['sequence_bp'], 200)
        self.assertEqual(audit['chr']['alignment_ambiguous_bp'], 100)
        self.assertEqual(sum(s['length_bp'] for s in segments), 200)

    def test_conflict_retains_construction_and_flags(self):
        support, audit, _ = fuse([atom(0, 100)], (), [(0, 80)])
        self.assertEqual(support['chr']['hap1'], 100)
        self.assertEqual(audit['chr']['conflict_bp'], 80)

    def test_ambiguous_and_unknown_orientation_are_quarantined(self):
        for kwargs in ({'status': 'ambiguous'}, {'strand': '?'}):
            support, audit, _ = fuse([atom(0, 100, **kwargs)], [(0, 100)])
            self.assertEqual(sum(support['chr'].values()), 0)
            self.assertEqual(audit['chr']['unsafe_agp_bp'], 100)

    def test_unmatched_agp_uses_alignment_fill_not_inferred_id(self):
        support, audit, _ = fuse([atom(0, 100, status='after_only')], (), [(0, 80)])
        self.assertEqual(support['chr'], {'hap1': 0, 'hap2': 80})
        self.assertEqual(audit['chr']['alignment_fill_bp'], 80)

    def test_main_fasta_uses_hybrid_and_retains_minimap_baseline(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'final.fa').write_text('>chr\n' + 'A'*100 + '\n>unplaced\nAA\n')
            agp = 'chr\t1\t100\t1\tW\th1tg000001l\t1\t100\t+\nunplaced\t1\t2\t1\tW\th2tg000002l\t1\t2\t+\n'
            for name in ('before.agp', 'after.agp'):
                (root/name).write_text(agp)
            (root/'h1.paf').write_text('')
            (root/'h2.paf').write_text('q\t100\t0\t100\t+\tchr\t100\t0\t100\t100\t100\t60\n')
            cmd = [sys.executable, '-m', 'hapreturn', '--final-fasta', str(root/'final.fa'),
                   '--hap1-paf', str(root/'h1.paf'), '--hap2-paf', str(root/'h2.paf'),
                   '--before-agp', str(root/'before.agp'), '--after-agp', str(root/'after.agp'),
                   '--chromosome-count', '1', '--min-aln-len', '1', '--no-plots', '--out-dir', str(root/'out')]
            run = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(run.returncode, 2, run.stderr)  # conflict requires review
            def rows(name):
                with (root/'out'/name).open(encoding='utf-8') as f:
                    return list(csv.DictReader(f, delimiter='\t'))
            self.assertEqual(rows('assignment.tsv')[0]['assignment'], 'hap1')
            self.assertEqual(rows('minimap_assignment.tsv')[0]['assignment'], 'hap2')
            self.assertEqual(rows('assignment.tsv')[0]['hap1_covered_bp'], '0')
            self.assertEqual(rows('assignment.tsv')[0]['hap1_support_bp'], '100')
            self.assertIn('>chr', (root/'out/hap1.fa').read_text())
            self.assertEqual((root/'out/hap2.fa').read_text(), '')
            self.assertIn('>unplaced', (root/'out/unplaced.fa').read_text())
            self.assertEqual(rows('hybrid_intervals.tsv')[0]['start'], '1')
            self.assertEqual(rows('evidence_comparison.tsv')[0]['conflict_bp'], '100')
            refresh = subprocess.run([sys.executable, '-m', 'hapreturn.report', '--result-dir', str(root/'out'),
                                      '--out-dir', str(root/'refresh'), '--min-aln-len', '1'], capture_output=True, text=True)
            self.assertEqual(refresh.returncode, 0, refresh.stderr)
            saved = json.loads((root/'refresh/curation/agp_comparison.json').read_text(encoding='utf-8'))
            compared = next(r for r in saved['agp_provenance']['method_comparison'] if r['seq_id']=='chr')
            self.assertEqual(compared['minimap_assignment'], 'hap2')
            self.assertIn('Joint evidence matrix', (root/'refresh/report.html').read_text(encoding='utf-8'))
            self.assertEqual((root/'refresh/evidence_comparison.tsv').read_bytes(), (root/'out/evidence_comparison.tsv').read_bytes())
            bad = subprocess.run(cmd + ['--force', '--assignment-mode', 'balanced', '--ploidy', '2'], capture_output=True, text=True)
            self.assertEqual(bad.returncode, 1)
            self.assertIn('explicit --homolog-groups', bad.stderr)
            fallback = subprocess.run(cmd + ['--force', '--evidence-mode', 'minimap', '--assignment-mode', 'provenance'], capture_output=True, text=True)
            self.assertEqual(fallback.returncode, 0, fallback.stderr)
            self.assertEqual(rows('assignment.tsv')[0]['assignment'], 'hap2')
            self.assertEqual(json.loads((root/'out/evidence_config.json').read_text())['evidence_mode'], 'minimap')


if __name__ == '__main__':
    unittest.main()
