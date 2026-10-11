# HapReturn

**Recover hifiasm hap1/hap2 provenance after merged assembly scaffolding.**

HapReturn aligns the two original hifiasm contig sets to a curated merged
assembly, measures source-specific chromosome coverage, and exports two source
sets with an auditable assignment table. It supports arbitrary chromosome
names and counts and optional ploidy and homolog-group checks.

```text
hifiasm hap1 + hap2 -> merged assembly -> HapHiC / Juicebox
       |                                      |
       +---------- provenance alignment ------+
                                              |
                                 HapReturn -> hap1.fa / hap2.fa
                                              + unresolved.fa / unplaced.fa
```

[中文说明](README.zh-CN.md) · [Method](docs/method.md) ·
[Release instructions](docs/github-upload.md)

## Install

Python 3.9 or newer is required. There are no Python runtime dependencies.
From a cloned or extracted repository:

```bash
python -m pip install .
hapreturn --version
```

When using source FASTAs, install minimap2 separately and ensure `minimap2`
is on PATH. Existing PAFs can be processed without minimap2, including on
Windows. Installation from PyPI is not provided by this release.

## Quick start

Run the bundled synthetic example without minimap2:

```bash
hapreturn \
  --final-fasta examples/final.fa \
  --hap1-paf examples/hap1.paf --hap2-paf examples/hap2.paf \
  --chromosome-list examples/chromosomes.txt \
  --homolog-groups examples/groups.tsv --ploidy 4 --hap1-per-group 2 \
  --min-aln-len 1 --out-dir example_output
```

Expected result: two chromosomes in each source set, one unplaced record and
zero unresolved chromosomes. The 1 bp alignment cutoff is for this small
example only; the default for real data is 10,000 bp.

For a real curated assembly:

```bash
hapreturn \
  --final-fasta out_JBAT.FINAL.fa \
  --hap1-fasta original.hap1.fa --hap2-fasta original.hap2.fa \
  --chromosome-list chromosomes.txt --threads 24 \
  --out-dir source_split
```

`chromosomes.txt` contains one reviewed chromosome FASTA ID per line. Records
outside that list go to `unplaced.fa`. Alternatively, `--chromosome-count N`
selects the first N records in FASTA order. With neither option, all records
are classified. No chromosome pairing is inferred from record order.

Each source must be specified by exactly one of `--hapN-fasta` or `--hapN-paf`.
PAF **targets must be the final assembly**, and queries must be original source
contigs. For manual alignment:

```bash
minimap2 -x asm5 --secondary=no -t 24 final.fa original.hap1.fa > hap1.paf
minimap2 -x asm5 --secondary=no -t 24 final.fa original.hap2.fa > hap2.paf
```

## Assignment and parameters

For each chromosome, let C1 and C2 be the union of accepted target intervals
from the two sources, measured in bases. The dominant source is accepted if
`max(C1,C2) / (C1+C2) >= dominance` and
`max(C1,C2) / chromosome_length >= min_covered_fraction`.

| Option | Default | Meaning |
|---|---:|---|
| `--min-mapq` | 10 | Minimum alignment mapping quality |
| `--min-aln-len` | 10000 | Minimum PAF alignment block length |
| `--dominance` | 0.60 | Minimum dominant-source share |
| `--min-covered-fraction` | 0.25 | Minimum dominant-source chromosome coverage |
| `--preset` | asm5 | minimap2 preset for source provenance alignment |
| `--threads` | 8 | minimap2 threads |
| `--allow-unresolved` | off | Accept partial outputs despite review findings |
| `--force` | off | Allow reuse of a nonempty output directory |

Secondary PAF hits (`tp:A:S`) are excluded. Overlapping intervals within each
source are counted once. Thresholds are screening rules and need validation
on the target data. Run `hapreturn --help` for all options.

## Different ploidies

Assignment always recovers the **two input source labels**. It does not require
equal numbers in those two source sets. The optional `--ploidy P` checks that
the selected chromosome count is divisible by P. When `--homolog-groups` is
provided, each group must have exactly P members.

The group TSV has columns `group_id` and `seq_id`, with one chromosome per row;
every selected chromosome must occur exactly once. `--hap1-per-group K`
additionally checks for K hap1 members in each group. For a balanced tetraploid
group, use P=4 and K=2; for a justified 1:2 triploid source balance, use P=3 and
K=1. These checks never force a chromosome into the other source.

## Outputs

| File | Contents |
|---|---|
| `hap1.fa`, `hap2.fa` | Confidently assigned chromosomes; identifiers preserved |
| `unresolved.fa` | Low-coverage or ambiguous chromosomes |
| `unplaced.fa` | Records outside the chromosome selection |
| `assignment.tsv` | Per-chromosome coverage, dominance, assignment and reason |
| `group_qc.tsv` | Optional group membership and source-balance checks |
| `hap1_to_final.paf`, `hap2_to_final.paf` | Alignments generated in FASTA mode |

Empty categories produce empty FASTAs. Exit codes: **0** accepted run, **1**
input/alignment error, **2** unresolved chromosomes or failed group QC. Exit 2
still produces reviewable outputs. A nonempty output directory requires
`--force`; use a new directory for each analysis.

## Scope and validation

HapReturn labels whole chromosomes by source provenance. In polyploids, a
source set can contain multiple biological haplotypes or subgenomes. Resolving
A/B subgenomes or A1/A2/B1/B2 sets requires further evidence and assignment.
Mixed-source scaffolds and local phase switches may be hidden by a dominant
whole-chromosome label; this release does not detect or split them. The assembly
FASTA is held in memory, so RAM requirements grow with assembly size. FASTA
inputs must be uncompressed.

Version 0.1.0 is an initial release validated on synthetic fixtures, including
triploid source imbalance, tetraploid groups, unresolved evidence and input
validation. Biological performance and the minimap2 execution path have not
been benchmarked on real datasets in this release.

```bash
python -m unittest discover -s tests -v
```

The implementation is derived from a Prunus assembly provenance workflow; see
[method and provenance](docs/method.md). Cite the hifiasm, HapHiC and minimap2
versions/publications relevant to your analysis. HapReturn currently has no
associated publication or archival DOI.

## License

MIT; see [LICENSE](LICENSE).
