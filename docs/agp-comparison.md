# AGP curation comparison

`hapreturn-agp` compares two assembly states using original component coordinates.
It does not require FASTA alignment or a Hi-C file.

## Usage

```bash
hapreturn-agp --before-agp before.agp --after-agp after.agp --out-dir agp_result
```

Use a new output directory, or `--force` to replace existing outputs.
`--no-plots` writes JSON/TSV only; `--min-display-bp` controls the overview cutoff,
while full interval tables retain all components.

## Inputs and interpretation

Both AGPs must use the same original sequence IDs and coordinate system. Input
validation checks nine columns, contiguous object coordinates, part numbers,
component spans, gaps and orientations. Component IDs are matched exactly.

The comparison distinguishes renaming, whole-scaffold reversal, scaffold order,
internal order/orientation, cross-scaffold changes, junctions and gaps. Duplicate
placements, unmatched intervals and unknown orientations are retained as QC
rather than forced into a unique correspondence. Scaffold anchors use reciprocal
unique largest overlap; ties remain unanchored.

Tables use **1-based inclusive** coordinates. The comparison describes endpoint
states, not the order of manual edits, and does not verify biological correctness
or component sequence identity.

## Source recovery

```bash
hapreturn-agp --before-agp before.agp --after-agp after.agp \
  --chromosome-list chromosomes.txt --final-fasta curated.fa --out-dir agp_result
```

Sources follow the `h1tg`/`h2tg` naming convention, optionally prefixed by `lcl|`.
An explicit `--source-map` TSV (`contig_id`, `source`) overrides name inference.
Use `--source-labels FIRST SECOND` for custom labels.

Independent source-majority assignment is the default. Non-gap component lengths
are counted by source, retaining mixed origins. `--agp-dominance` defaults to 0.60
and `--agp-min-known-fraction` to 0.25. Balanced assignment requires
`--agp-assignment-mode balanced`, explicit `--homolog-groups` and even `--ploidy`.

`--chromosome-count N` selects the first N AGP objects; a chromosome list selects
explicit IDs. With `--final-fasta`, whole records are exported to `agp_partition/`
as the two source sets, unresolved and unplaced. FASTA IDs and lengths must match
the final AGP.

To compare alignment evidence, supply `--minimap-assignment` with an independent
minimap assignment table, or both `--hap1-paf` and `--hap2-paf`. Default PAF filters
are MAPQ ≥10 and alignment length ≥10 kb. Missing alignment evidence is reported
as unavailable.

## Outputs

- `segment_correspondence.tsv`, `scaffold_correspondence.tsv`: interval and scaffold mappings.
- `curation_events.tsv`, `junction_changes.tsv`, `component_boundary_changes.tsv`: structural differences.
- `curation_qc.tsv`, `curation_summary.tsv`, `agp_comparison.json`: validation and reusable data.
- `agp_chromosome_provenance.tsv`, `agp_chromosome_components.tsv`, `agp_original_scaffold_composition.tsv`: original-source composition.
- `agp_assignment.tsv`, `agp_partition_summary.tsv`, `agp_group_qc.tsv`, `agp_vs_minimap.tsv`: assignments and comparisons.
- `curation_report.html`, correspondence and provenance PDF/SVG files: interactive and static figures.

Exit status: 0 completed; 1 input/runtime error; 2 output with review items.
`--allow-review` accepts review-required output.

For combined AGP and minimap2 assignment through the main command, see
[hybrid evidence](hybrid-evidence.md).
