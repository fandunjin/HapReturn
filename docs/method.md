# Method and provenance

## Source workflow

HapReturn was extracted from a local Prunus project:

1. `jobs/haphic.sh` merged cleaned hifiasm hap1 and hap2 contigs, mapped Hi-C
   reads and scaffolded the merged assembly with HapHiC.
2. `jobs/finalize_assembly.sh` / `scripts/finalize_assembly.py` accepted a
   Juicebox-curated assembly, selected chromosomes and generated an adjacent
   homolog-pair table. Those chromosome-order assumptions are project specific.
3. `jobs/final_haplotype_provenance.sh` aligned the original contigs to final
   chromosomes; `scripts/haplotype_provenance.py` used interval-union coverage
   to orient homolog pairs and emit two phase-copy FASTAs.
4. The project's later self-synteny and reference-constrained assignment
   organized chromosomes into A1/A2/B1/B2 sets. Those steps are outside this
   package. The P. spinulosa scripts used a related source-provenance workflow;
   the P. padus partition used an existing public assembly.

The project scripts and research data are not bundled. The standalone package
removes fixed chromosome counts and pair-order assumptions and reports ambiguous
assignments rather than inheriting an original pair order.

## Coverage calculation

PAF uses zero-based, half-open target intervals. After mapQ, block-length and
secondary-hit filtering, intervals are sorted and merged within each source
and chromosome. Let their union lengths be C1 and C2 and chromosome length be L:

```text
dominance = max(C1, C2) / (C1 + C2)
dominant_coverage_fraction = max(C1, C2) / L
```

No coverage or equal source coverage is unresolved. Coverage below the minimum
is `low_coverage`; insufficient dominance is `ambiguous_origin`. PAF target
names and lengths must agree with the supplied final FASTA. Overlap between
the two sources contributes to each source's own union; dominance is therefore
a relative coverage statistic, not a posterior probability of phase.

## Interpretation

Group QC checks user-supplied membership and optional source quotas. It does
not discover homologues or optimize a constrained assignment. Ploidy does not
change the coverage classifier. The same mathematical classifier can process
diploid or polyploid inputs with two source sets, but that does not establish
its biological accuracy across species or ploidies.

Thresholds, genome repetitiveness, structural changes, collapsed sequence and
the source assembler's phase errors all affect the result. Whole-chromosome
labels may conceal local mixed provenance. Validate chromosome structure and
local phase independently before interpreting the FASTAs biologically.
