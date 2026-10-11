# Method

## Two-source recovery

HapReturn recovers two n chromosome sets from a merged 2n assembly by tracing
chromosomes to the original hifiasm hap1 and hap2 contigs. The assignment method
is the same at every ploidy; n denotes half of the input chromosome count.

1. Supply the curated merged assembly and the two original source assemblies
   or their PAF alignments.
2. Align each source to the final assembly with minimap2, with secondary hits
   disabled.
3. Filter by mapping quality and alignment length, then merge overlapping target
   intervals separately for hap1 and hap2 on each chromosome.
4. Assign each chromosome to its supported original source.
5. Export the two source FASTAs and the chromosome evidence table.

Place all chromosome scaffolds first during Juicebox curation, then select the
first N exported FASTA records with `--chromosome-count N`, where N is the total
chromosome count across both sources. Selection follows FASTA order without
length sorting; unselected records are exported to `unplaced.fa`. A chromosome
identifier list is available when records are interspersed.

## Coverage calculation

PAF target coordinates are zero-based, half-open intervals. For chromosome
length L and interval-union coverages C1 and C2:

```text
dominance = max(C1, C2) / (C1 + C2)
dominant_coverage_fraction = max(C1, C2) / L
```

A unique source is assigned when dominance is at least 0.60 and its chromosome
coverage is at least 0.25 by default. Ties, low coverage and insufficient
dominance are written to `unresolved.fa`. Secondary PAF hits are excluded and
overlapping intervals within each source are counted once.

The expected chromosome counts are n in each recovered set. Labels follow
alignment evidence; unresolved records remain separate rather than being forced
into a source to achieve the expected counts.

## Optional checks and export

User-supplied homolog groups and copy counts can validate a split. The ploidy
parameter checks chromosome counts and group sizes only; it is not required
for source assignment. The output set count is always two.

The assembly is scanned for identifiers and lengths and then streamed to the
output files. Original chromosome identifiers are preserved. Each record goes
to its source set, the unresolved category or the unplaced category.

## Visualization

For each target window, intersect the union of accepted alignment intervals from
each source with the window. Apply the chromosome coverage and dominance
thresholds to these local counts. The final window uses its actual length.
Whole-chromosome labels remain based on whole-chromosome union coverage.
Coverage is displayed independently for each source. Chromosome outlines and
original-to-curated ribbons use shared length scales. Summary plots include all output sequences, including
unresolved and unselected records; coverage plots include selected records only.
Window coordinates in the evidence table are zero-based and half-open.

Matrix cells equal the union of target intervals from one original contig divided
by target chromosome length; contig identifiers are namespaced by original
source. Rows are grouped by source and ordered by the best-supported target,
then identifier. Every retained contig with a selected target is included.
Ribbons map PAF query and target intervals to their respective sequence bodies,
reversing target endpoints for minus-strand alignments. Inter-sequence gaps reserve
label space and are excluded from the shared bp scale. Separate PDFs and one
HTML report are exported; the partition PDF contains source evidence and
chromosome ideograms as two figures, with additional pages for large inputs.
