# HapReturn

**Recover two chromosome sets from a curated merged assembly using original hifiasm provenance and equal-copy homolog-group constraints.**

[中文说明](README.zh-CN.md) · [Method](docs/method.md)

## Install

Python 3.9+ is required. Install minimap2 for FASTA alignment; existing PAFs can be reused. Matplotlib is installed automatically.

```bash
python -m pip install .
hapreturn --version
```

## Input order

Place chromosome scaffolds first during Juicebox curation and arrange homolog groups consecutively. A tetraploid input has records 1–4 in the first group, 5–8 in the second, and so on. An octoploid input has eight records per group. The equal-copy assumption must match the input design; no fixed subgenome labels are required.

`--chromosome-count N` selects the first N FASTA records without length sorting. Remaining records go to `unplaced.fa`. Alternatively, supply `--homolog-groups groups.tsv` with columns `group_id` and `seq_id`, listing all selected chromosomes exactly once.

## Run

```bash
hapreturn \
  --final-fasta curated.fa \
  --hap1-fasta original.hap1.fa --hap2-fasta original.hap2.fa \
  --chromosome-count 56 --ploidy 4 --threads 24 \
  --out-dir result
```

The default balanced mode assigns two of every four homologs to each source, producing 28 + 28 records. At ploidy eight it assigns 4 + 4 per group. An even `--ploidy` is required.

Use `--hap1-paf` and `--hap2-paf` instead of source FASTAs to reuse alignments. Queries must be original contigs; targets must be the curated assembly. Exactly two custom labels can be supplied with `--source-paf LABEL=PATH` or `--source-fasta LABEL=PATH`.

## Assignment

After filtering at mapQ ≥10 and alignment length ≥10 kb, merge overlapping target intervals within each source. For each chromosome, C1 and C2 are union-covered target bases. Assign the P/2 largest differences C1−C2 to source 1 and the remainder to source 2. This exactly maximizes total support under the group balance constraint.

`assignment.tsv` retains final labels, independent evidence labels, group IDs, and final-source coverage. Independent evidence uses coverage ≥0.25 and dominance ≥0.60 by default. `group_qc.tsv` records total support, the optimality gap, changed evidence, and tied optima. Ties use stable input order and are flagged for review; zero evidence is recorded as such.

Use `--assignment-mode provenance` for independent source tracing with potentially unequal output counts. In that mode, groups and copy quotas provide optional QC only.

## Local unresolved windows

Within each original contiguous gray run, expand an unresolved window's evidence context geometrically up to the full run length. Reuse existing PAF evidence, preserve raw coverage, and do not borrow evidence from resolved neighboring runs. Retain gray when the expanded context does not meet coverage and dominance thresholds.

`window_provenance.tsv` records raw and revised labels, method and evidence span. `window_rescue.tsv` records every context attempt. Local labels are separate from chromosome group assignments.

## Outputs

- `hap1.fa`, `hap2.fa`, `unresolved.fa`, `unplaced.fa`: sequence partitions.
- `assignment.tsv`, `summary.tsv`, `homolog_groups.tsv`, `group_qc.tsv`: assignment evidence and checks.
- `partition.pdf`: source evidence, counts, lengths and chromosome provenance.
- `alignment_matrix.pdf`: original-contig × curated-chromosome support.
- `hap_comparison.pdf`: one continuous, height-adaptive correspondence figure per recovered source.
- `report.html`: all figures, source statistics, window hover and numbered filtering.

The correspondence figures retain accepted links from both original sources. PDFs use editable text; figures are exported as separate PDFs and one portable HTML.

### Refresh an existing report

```bash
hapreturn-report --result-dir result --out-dir refreshed_report
```

This preserves chromosome assignments and redraws figures with local window expansion. To apply balanced assignment to an earlier result, rerun the main command with the original curated FASTA, existing PAFs and ploidy, using a new output directory. Reuse any customized filtering and window settings.

### Standalone correspondence figure

```bash
hapreturn-chromosomes \
  --final-fasta curated.fa --hap1-paf hap1.paf --hap2-paf hap2.paf \
  --assignment result/assignment.tsv --out-pdf original_to_curated.pdf
```

## Synthetic octoploid example

```bash
hapreturn --final-fasta examples/octoploid/final.fa \
  --hap1-paf examples/octoploid/hap1.paf --hap2-paf examples/octoploid/hap2.paf \
  --ploidy 8 --min-aln-len 1 --out-dir example_result
```

Exit status: 0 passed; 1 input/runtime error; 2 outputs available with review flags. `--allow-unresolved` accepts review-required output. MIT license.
