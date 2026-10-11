# Figure design and source data

## Scientific question

Which original hifiasm source supports each recovered chromosome, and where is
that evidence consistent or mixed along the curated assembly?

## Figure contract

- `partition.pdf`, figure 1: independent source coverage and recovered output
  counts/lengths establish chromosome assignments and the extent of recovery.
- `partition.pdf`, figure 2: length-scaled chromosome ideograms show local source
  evidence. Rounded ends depict sequence boundaries; no centromere is inferred.
- `alignment_matrix.pdf`: original source contigs by curated chromosomes show
  the target fraction covered by each contig, with source and assignment keys.
  Each recovered source has one continuous matrix with all original-contig rows
  from both sources and all chromosomes assigned to that recovered source.
- `hap_comparison.pdf`: coordinate-aware ribbons link original hifiasm contigs
  to recovered chromosomes; color always denotes the original source and
  reverse-strand links retain their orientation.

These are descriptive assembly figures, not replicated biological comparisons.
Counts are sequence counts; no inferential tests or error bars are applicable.
Whole-chromosome balanced assignments maximize support under group copy quotas;
independent and local window evidence use the reported coverage and dominance
thresholds. Overlapping intervals from the same source (or the
same contig in matrix cells) are counted once. Off-source alignments are evidence
of sequence similarity, not an estimate of bases moved between haplotypes.

## Format and integrity

Python/matplotlib produces vector PDF with embedded TrueType fonts, 5–7 pt
sans-serif text, 8 pt lowercase panel labels, white backgrounds, and one
consistent blue/orange source palette. Partition figures use a width of 183 mm
and page heights up to 170 mm. Alignment matrices increase width and height
with record count; correspondence figures increase height. Each recovered set
remains continuous without reducing label sizes. HTML embeds the same rendered
figures and searchable quantitative evidence without external resources.

Machine-readable source data include `assignment.tsv`, `summary.tsv`,
`window_provenance.tsv`, `window_rescue.tsv`, `window_qc.tsv`,
`alignment_matrix.tsv` and `alignment_links.tsv`. The full matrix is exported
with all selected chromosomes. No contigs are removed to simplify the figure.

Formatting follows the [Nature figure specifications](https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/)
and [panel export guidance](https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/).
Synthetic examples are explicitly labelled as demonstrations.
