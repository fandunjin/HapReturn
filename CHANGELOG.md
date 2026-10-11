# Changelog

## 0.3.0

- Automatically export source coverage, assignment, window provenance and output summary figures in PNG, SVG and PDF.
- Add a portable HTML report with window hover details and assignment filtering.
- Export window evidence as TSV; add configurable window size and an option to skip plots.

## 0.2.1

- Use one two-source recovery workflow at every ploidy: 2n to n + n.
- Require exactly two original sources and simplify the examples and methods.
- Verify unchanged commands on 2-, 4-, 8-, 12- and 16-ploid input sizes.

## 0.2.0

- Introduce named FASTA/PAF inputs and TSV source manifests.
- Add source-specific copy counts for arbitrary-ploidy homolog groups.
- Stream assembly sequences to output files and export per-source summaries.
- Add octoploid, 12-ploid and 16-ploid tests.
- Simplify documentation with species-independent methods and source labels.

## 0.1.0

- Recover two-source chromosome provenance using minimap2 or existing PAFs.
- Use union coverage to avoid counting overlapping alignments repeatedly.
- Export hap1, hap2, unresolved and unplaced sequences with original identifiers.
- Validate optional chromosome lists, ploidy and homolog-group source balance.
- Provide an installable CLI, synthetic examples and automated tests.
