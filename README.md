# HapReturn

**Recover source-labelled chromosome sets after merged assembly scaffolding.**

HapReturn maps original source assemblies to a curated merged assembly and
partitions its chromosomes by alignment coverage. Use hifiasm hap1/hap2 inputs
or supply any number of sources with your own labels. Diploid, tetraploid,
octoploid and higher-ploidy assemblies use the same workflow.

[中文说明](README.zh-CN.md) · [Method](docs/method.md)

## Install

Requires Python 3.9+. Install minimap2 for FASTA alignment; existing PAFs can be
processed directly without additional Python dependencies.

```bash
python -m pip install .
hapreturn --version
```

## Return a merged assembly to hap1 and hap2

```bash
hapreturn \
  --final-fasta curated.fa \
  --hap1-fasta original.hap1.fa --hap2-fasta original.hap2.fa \
  --chromosome-list chromosomes.txt \
  --threads 24 --out-dir source_split
```

Use `--hap1-paf` and `--hap2-paf` for existing alignments instead. PAF queries
are original source contigs and targets are the curated assembly:

```bash
minimap2 -x asm5 --secondary=no -t 24 curated.fa original.hap1.fa > hap1.paf
```

`chromosomes.txt` lists one chromosome ID per line. Alternatively,
`--chromosome-count N` selects the first N FASTA records; omit both options to
classify all records. Unselected records are written to `unplaced.fa`.

## Octoploid and higher-ploidy assemblies

For an octoploid assembly with two source sets, add homolog groups and copy
counts to the command above:

```text
--ploidy 8 --homolog-groups groups.tsv --source-copies hap1=4 --source-copies hap2=4
```

`groups.tsv` contains `group_id` and `seq_id` columns, with one chromosome per
row. Each octoploid group has eight members. Use `--ploidy 12` or `--ploidy 16`
for higher ploidies and set the source copy counts accordingly. Unequal source
copy counts are also supported.

For multiple known sources, use a TSV manifest:

```text
source  path         format
copy01  copy01.paf   paf
copy02  copy02.paf   paf
...
copy08  copy08.paf   paf
```

Separate columns with tabs. Paths are relative to the manifest; `format` is
`paf` or `fasta`. Each source gets its own `<source>.fa` output. Run the bundled
eight-source, two-group octoploid example:

```bash
hapreturn \
  --final-fasta examples/octoploid/final.fa \
  --sources examples/octoploid/sources.tsv \
  --homolog-groups examples/octoploid/groups.tsv --ploidy 8 \
  --min-aln-len 1 --out-dir octoploid_output
```

Expected: eight source FASTAs, two chromosomes each. The 1 bp cutoff is for
the small synthetic example. Source labels describe the supplied assemblies;
two input sources produce two source sets, even when the organism is octoploid.

Sources can also be supplied individually:

```bash
hapreturn --final-fasta curated.fa \
  --source-paf copy01=copy01.paf --source-paf copy02=copy02.paf \
  --out-dir result
```

Use repeatable `--source-fasta LABEL=PATH` for FASTA inputs. Ploidy and source
count have no fixed upper limit in the assignment logic.

## Parameters

| Option | Default | Purpose |
|---|---:|---|
| `--min-mapq` | 10 | Minimum PAF mapping quality |
| `--min-aln-len` | 10000 | Minimum alignment block length in bases |
| `--dominance` | 0.60 | Minimum dominant-source share of covered bases |
| `--min-covered-fraction` | 0.25 | Minimum dominant-source chromosome coverage |
| `--ploidy` | unset | Check chromosome count and homolog-group size |
| `--source-copies LABEL=N` | unset | Check source copies per homolog group; repeatable |
| `--preset` | asm5 | minimap2 preset |
| `--threads` | 8 | minimap2 threads |
| `--allow-unresolved` | off | Accept partial results |
| `--force` | off | Reuse a nonempty output directory |

The existing `--hap1-per-group` option remains available. FASTA inputs are
uncompressed; output writing streams sequence lines to support large assemblies.

## Outputs

| File | Contents |
|---|---|
| `<source>.fa` | Assigned chromosomes with their original identifiers |
| `unresolved.fa` | Chromosomes needing further evidence |
| `unplaced.fa` | Unselected assembly records |
| `assignment.tsv` | Source-specific coverage and chromosome assignments |
| `group_qc.tsv` | Homolog-group source counts and copy-count checks |
| `summary.tsv` | Sequence counts and total bases per output |
| `<source>_to_final.paf` | Alignments generated from FASTA inputs |

Exit codes: 0 for success, 1 for an input/alignment error, 2 for unresolved
chromosomes or failed group checks. Outputs are retained for review on exit 2.
HapReturn assigns whole chromosome records; local phase-switch analysis can
be performed downstream.

## Tests

```bash
python -m unittest discover -s tests -v
```

Tests cover two-source and multisource octoploid, 12-ploid and 16-ploid inputs,
source copy counts, interval overlap, ambiguous evidence and input validation.

## License

[MIT](LICENSE).
