# Method

## Ordered homolog groups

The selected chromosome records must be ordered in consecutive groups of even size P, the specified ploidy, or supplied in an explicit group table. Each group is assumed to contain P homologous chromosomes with equal representation in the two recovered sets. Unselected records remain unplaced.

## Source evidence

Original hifiasm source contigs are aligned to the curated merged assembly with minimap2, with secondary hits disabled. Accepted PAF records pass mapping-quality and alignment-length filters. Target intervals are merged separately for each source and chromosome to count covered bases once within each source. Overlap between different sources is retained as competing evidence.

For chromosome i, length L_i and coverages C_i1, C_i2, independent source evidence is accepted when max(C_i1,C_i2)/L_i ≥0.25 and max(C_i1,C_i2)/(C_i1+C_i2) ≥0.60, with a unique winner. Otherwise its independent evidence label is unresolved.

## Exact balanced assignment

For each group of P chromosomes, choose P/2 chromosomes for source 1 and the remainder for source 2 to maximize:

```text
S = sum(C_i1 for i assigned to source 1)
  + sum(C_i2 for i assigned to source 2)
```

Start with every chromosome assigned to source 2. Moving chromosome i to source 1 changes S by D_i=C_i1−C_i2. Therefore selecting the P/2 largest D_i is an exact optimum; exhaustive enumeration is unnecessary. Complexity is O(P log P) per group. The score uses covered bases, not length-normalized fractions.

The difference between the P/2-th and the next ranked D_i is the support gap to the next feasible assignment. A zero gap indicates tied optima. Ties are resolved reproducibly in input order and flagged for review. Assignments lacking independent support or differing from the independent label also receive review flags. Equal-copy output is retained, including in these cases; its constrained nature remains visible in the evidence table. Original sequence IDs and sequences are preserved.

Independent tracing is available with `--assignment-mode provenance`; it does not enforce equal counts.

## Expanding local unresolved windows

Partition each selected chromosome into base windows, using the actual length of the final partial window. Record union-covered bases from each source. Apply the independent coverage and dominance thresholds to obtain raw window labels.

Identify runs of consecutive raw unresolved windows. For every window in a run, try contexts spanning 2, 4, 8, ... base windows, clipped and shifted to remain inside that original run. The last context spans the complete run. Aggregate original covered bases over each context and apply the same thresholds. At the first accepted context, record its source as an expanded-context local label; otherwise retain unresolved. Raw window covered bases are never replaced by context totals. The process reuses existing alignments rather than remapping sequence. No evidence from neighboring resolved runs or from group copy quotas enters this local inference.

`window_rescue.tsv` records attempted contexts, their source-covered bases and decisions. `window_provenance.tsv` preserves raw and revised labels, the method, and context coordinates. Thus a context-supported label is distinguishable from direct window evidence. Coordinates are zero-based and half-open.

## Figures

Source-support rows use faint horizontal guides. Count and assembly-length labels use plain numeric formatting. Chromosome outlines share a length scale without inferred centromeres. Alignment-matrix cells are union-covered target bases from one original contig divided by target length. Correspondence ribbons preserve PAF query and target intervals and strand. Every recovered source is drawn in one continuous figure with height adapted to record count; both sides share a base-pair scale and source colors follow the original contigs. HTML retains numbered filtering, source counts, raw evidence and constrained labels.
