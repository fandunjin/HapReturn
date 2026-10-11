# Changelog

## 0.1.0

- Recover two-source chromosome provenance using minimap2 or existing PAFs.
- Use union coverage to avoid counting overlapping alignments repeatedly.
- Export hap1, hap2, unresolved and unplaced sequences with original identifiers.
- Validate optional chromosome lists, ploidy and homolog-group source balance.
- Provide an installable CLI, synthetic examples and automated tests.

This initial release has synthetic validation. Real assembly benchmarking is
required before making biological accuracy claims.
