# Hybrid evidence

The main `hapreturn` command combines AGP construction provenance and minimap2
support when both `--before-agp` and `--after-agp` are supplied. Main FASTA outputs
and `assignment.tsv` use the integrated evidence. Without AGPs, the existing
alignment-only workflow remains available.

## Usage

```bash
hapreturn --final-fasta curated.fa \
  --hap1-paf hap1.paf --hap2-paf hap2.paf \
  --before-agp before.agp --after-agp after.agp \
  --chromosome-list chromosomes.txt --out-dir result
```

`--evidence-mode auto` selects hybrid when AGPs are present; `minimap` explicitly
uses alignment-only assignment. Hybrid defaults to independent provenance
assignment. Balanced hybrid assignment requires `--assignment-mode balanced`,
`--homolog-groups groups.tsv` and an even `--ploidy`.

## Evidence rules

- Use the same original component IDs and coordinates in both AGPs. Final AGP
  object IDs and lengths must match the final FASTA.
- Unique, direction-resolved AGP placements with a known source provide direct
  support. Infer sources from `h1tg`/`h2tg` IDs or supply `--agp-source-map` with
  `contig_id` and `source` columns matching the alignment source labels.
- Accepted alignments fill unknown component intervals only when exactly one
  source covers the interval. Gaps never contribute; overlapping evidence is
  counted once. Repeated placements and unknown orientations remain unresolved.
- Opposite-source-only alignment coverage on a known AGP interval is flagged as
  a conflict while retaining the AGP origin.
- Independent assignment requires known support / non-gap length ≥0.25 and
  dominant-source support / total known support ≥0.60 by default. Use
  `--min-covered-fraction` and `--dominance` to change these thresholds.
- A whole-chromosome label preserves its minority-source intervals and sequence.
  Group constraints retain independent labels and flag changed decisions or ties.

AGP records construction history; minimap2 measures sequence similarity. These
are related evidence sources, not independent votes. ID/coordinate validation
does not reconstruct component sequences or establish parental phase continuity.

## Outputs

| File | Contents |
|---|---|
| `assignment.tsv` | Final and independent labels, evidence mode and source support |
| `minimap_assignment.tsv` | Independent alignment-only baseline |
| `evidence_comparison.tsv` | AGP, minimap2 and integrated labels; support and review flags |
| `hybrid_intervals.tsv` | Disjoint interval evidence; 1-based inclusive coordinates |
| `evidence_config.json` | Thresholds, selected records, input checksums and group settings |
| `report.html` | Unified English evidence report; `evidence_report.html` is a compatibility entry |
| `joint_evidence_tracks.*` | Co-registered AGP, alignment and integrated tracks; PDF/SVG/PNG |
| `joint_evidence_matrix.*` | Paired AGP/alignment matrix; PDF/SVG/PNG and source TSV |
| `curation/` | AGP correspondence, composition and curation audit |

In `assignment.tsv`, `*_covered_bp` and `dominant_coverage_fraction` retain
alignment statistics; `*_support_bp` contains the evidence used for assignment.
For hybrid runs, `support_coverage_fraction` uses known support / non-gap length.
The joint matrix uses full target chromosome length for both evidence fractions.

`hapreturn-report` redraws the saved result without changing assignments.
To change assignment thresholds, rerun the main command. Exit status 2 indicates
outputs requiring review; `--allow-unresolved` accepts them.
