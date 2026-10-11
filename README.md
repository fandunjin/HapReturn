# HapReturn

**Return a merged 2n assembly to two n sets using original hifiasm hap1/hap2 provenance.**

HapReturn aligns the two original hifiasm contig sets to a curated merged
assembly and recovers the chromosome sets belonging to hap1 and hap2.
The workflow is identical at any ploidy: **2n → n + n**.

```text
original hap1 + original hap2 -> merge -> scaffold and curate -> final 2n
      |                |                                         |
      +----------------+----------- HapReturn -------------------+
                                      |
                              hap1 (n) + hap2 (n)
```

[中文说明](README.zh-CN.md) · [Method](docs/method.md)

## Install

Requires Python 3.9+. Install minimap2 to align FASTA inputs; existing PAFs can
be processed without minimap2. Matplotlib is installed with HapReturn.

```bash
python -m pip install .
hapreturn --version
```

## Run

Use the same command for diploid, tetraploid, octoploid and higher-ploidy inputs:

```bash
hapreturn \
  --final-fasta curated.fa \
  --hap1-fasta original.hap1.fa --hap2-fasta original.hap2.fa \
  --chromosome-list chromosomes.txt --threads 24 \
  --out-dir result
```

`chromosomes.txt` lists the 2n chromosome identifiers, one per line. Alternatively,
`--chromosome-count N` selects the first N FASTA records; omit both options to
classify all records. Unselected records go to `unplaced.fa`.
The number n is half of the selected chromosome count: for 64 chromosomes,
the expected split is 32 + 32. No ploidy setting is required for assignment.

Use `--hap1-paf` and `--hap2-paf` to supply existing alignments instead of
source FASTAs. PAF queries are original contigs and targets are the final assembly:

```bash
minimap2 -x asm5 --secondary=no -t 24 curated.fa original.hap1.fa > hap1.paf
```

Run the bundled two-source example:

```bash
hapreturn \
  --final-fasta examples/octoploid/final.fa \
  --hap1-paf examples/octoploid/hap1.paf \
  --hap2-paf examples/octoploid/hap2.paf \
  --min-aln-len 1 --out-dir example_result
```

Expected: 16 input chromosomes split into 8 hap1 and 8 hap2 chromosomes.
The 1 bp cutoff is for this small synthetic example.

## Publication figures

Each run writes **three independent PDFs and one self-contained HTML report**:

| File | Contents |
|---|---|
| `partition.pdf` | Two figures: source support with count/length summaries; length-scaled chromosome provenance ideograms |
| `alignment_matrix.pdf` | Original contig × curated chromosome heatmap, with source and recovered-assignment keys |
| `hap_comparison.pdf` | Original hifiasm contigs linked to recovered chromosomes by oriented, coordinate-aware ribbons |
| `report.html` | All figures, window hover details and a searchable assignment table |

The PDFs use editable, embedded TrueType text, 5–7 pt sans-serif labels and
183 mm page width. Large datasets are paginated at readable font sizes; the
complete matrix and links are retained as TSV. PNG/SVG files are not emitted.
Coverage and output lengths use separate axes. Blue/orange always encode the
two original sources; gray denotes unresolved evidence. Chromosome outlines
have rounded sequence ends without an inferred centromere. All ploidies use
the same two-source logic and shared length scales.

![Synthetic chromosome provenance example](docs/images/chromosome_provenance.png)

`--plot-window-size N` sets window width in bp (default: about 100 windows on
the longest chromosome); `--no-plots` skips visualization. Window coverage is
calculated from interval unions, using the same local coverage/dominance
thresholds. Reports are generated even when assignment QC returns exit code 2.
See [figure design and source data](docs/figure-design.md).

Reproduce the synthetic demonstration:

```bash
python examples/visualization/generate.py demo_input
hapreturn --final-fasta demo_input/final.fa --hap1-paf demo_input/hap1.paf \
  --hap2-paf demo_input/hap2.paf --min-aln-len 1 --allow-unresolved --out-dir demo_report
```

### Standalone original-to-curated chromosome plot

Use the two original-to-curated PAFs and an existing HapReturn assignment table:

```bash
hapreturn-chromosomes \
  --final-fasta curated.fa --hap1-paf hap1.paf --hap2-paf hap2.paf \
  --assignment result/assignment.tsv --out-pdf original_to_curated.pdf
```

The same command is available as `python scripts/plot_chromosomes.py` after
installation. `--source-labels FIRST SECOND` supports two custom source names.
Set `--min-mapq` and `--min-aln-len` to the thresholds used for assignment.
The original side shows accepted aligned contigs, with original hap labels;
the curated side shows chromosome assignments. Ribbon endpoints retain PAF
query/target coordinates and strand. Off-source links show sequence similarity,
not a measurement of transferred bases.

## Parameters

| Option | Default | Purpose |
|---|---:|---|
| `--min-mapq` | 10 | Minimum PAF mapping quality |
| `--min-aln-len` | 10000 | Minimum alignment block length in bases |
| `--dominance` | 0.60 | Minimum dominant-source share of covered bases |
| `--min-covered-fraction` | 0.25 | Minimum dominant-source chromosome coverage |
| `--preset` | asm5 | minimap2 preset |
| `--threads` | 8 | minimap2 threads |
| `--allow-unresolved` | off | Accept partial results |
| `--force` | off | Reuse a nonempty output directory |
| `--plot-window-size` | auto | Source window width in bp |
| `--no-plots` | off | Skip PDFs, HTML and visualization source tables |

Optional `--homolog-groups`, `--ploidy` and source-copy counts provide additional
group checks; they do not change the assignment logic. Exactly two named sources
can also be supplied with `--sources` or `--source-paf` / `--source-fasta`.
FASTA inputs are uncompressed and sequence output is streamed.

## Outputs

| File | Contents |
|---|---|
| `hap1.fa`, `hap2.fa` | Recovered source sets with original chromosome identifiers |
| `unresolved.fa` | Chromosomes with insufficient or conflicting coverage |
| `unplaced.fa` | Unselected assembly records |
| `assignment.tsv` | Source coverage and chromosome assignments |
| `summary.tsv` | Sequence counts and total bases per output |
| `group_qc.tsv` | Optional homolog-group checks |
| `window_provenance.tsv` | Local source coverage and labels, zero-based half-open coordinates |
| `alignment_matrix.tsv` | Full target-normalized contig-to-chromosome coverage matrix |
| `alignment_links.tsv` | Original source, query/target coordinates, strand and recovered assignment |
| `hap1_to_final.paf`, `hap2_to_final.paf` | Alignments generated from FASTA inputs |

Exit codes: 0 for success, 1 for an input/alignment error, 2 for unresolved
chromosomes or failed group checks. Outputs are retained on exit 2; source
assignments follow alignment evidence rather than forcing equal counts.

## Tests

```bash
python -m unittest discover -s tests -v
```

Tests verify the same two-source command at 2-, 4-, 8-, 12- and 16-ploid input
sizes, together with overlap handling, ambiguous evidence and input checks.

## License

[MIT](LICENSE).
