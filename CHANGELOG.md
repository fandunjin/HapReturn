# Changelog

## 0.2.0

- Accept any number of named FASTA or PAF sources and TSV source manifests.
- Add source-specific copy counts for arbitrary-ploidy homolog groups.
- Stream assembly sequences to output files and export per-source summaries.
- Add octoploid, 12-ploid and 16-ploid tests and an eight-source example.
- Simplify documentation with species-independent methods and source labels.

## 0.1.0

- Recover two-source chromosome provenance using minimap2 or existing PAFs.
- Use union coverage to avoid counting overlapping alignments repeatedly.
- Export hap1, hap2, unresolved and unplaced sequences with original identifiers.
- Validate optional chromosome lists, ploidy and homolog-group source balance.
- Provide an installable CLI, synthetic examples and automated tests.
