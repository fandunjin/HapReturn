# Method

## Workflow

1. Supply the curated merged assembly and two or more original source assemblies
   or their PAF alignments. Source labels are user-defined.
2. Align each FASTA source to the curated assembly using minimap2 with the
   selected preset and secondary alignments disabled.
3. Filter PAF records by mapping quality and alignment length and merge target
   intervals within each source and chromosome.
4. Assign each chromosome to the source with the strongest accepted coverage.
5. Export source FASTAs, chromosome evidence, source summaries and optional
   homolog-group checks.

## Coverage and assignment

PAF target coordinates are zero-based, half-open intervals. Let Ci be the
interval-union coverage from source i, L the chromosome length, and Cmax the
largest source coverage:

```text
dominance = Cmax / sum(Ci)
dominant_coverage_fraction = Cmax / L
```

A unique source is assigned when dominance is at least 0.60 and its chromosome
coverage is at least 0.25 by default. Ties, low coverage and insufficient
dominance produce an unresolved record. Secondary hits marked `tp:A:S` are
excluded, and overlapping intervals within a source contribute bases once.

`assignment.tsv` includes one coverage column per input source, followed by
the dominant coverage fraction, dominance, assigned label and reason.

## Ploidy and homolog groups

The coverage calculation operates over an arbitrary number of sources and
chromosomes. Ploidy P checks that the selected chromosome count is divisible
by P. With user-supplied homolog groups, each group must contain P chromosomes;
optional source-copy quotas check how many are assigned to each source.

For example, an octoploid group may contain four chromosomes from each of two
sources, or one chromosome from each of eight known sources. The same rules
apply to 12-ploid, 16-ploid and higher-ploidy assemblies. Group names and source
labels are independent of biological subgenome naming.

## Sequence export

The assembly is scanned to collect identifiers and lengths, then streamed to
the output FASTAs. Sequence bases from the full assembly are not retained in
memory. Original identifiers are preserved, and each record is exported once
to its assigned source, unresolved category or unplaced category.
