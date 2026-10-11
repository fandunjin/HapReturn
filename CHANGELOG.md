# Changelog

## 0.5.2

- Simplify HTML by removing the group-assignment and local-reassessment explanation sections.
- Reserve measured chromosome-label space between provenance rows and adapt page height to prevent overlap.

## 0.5.1

- Draw a complete alignment matrix for each recovered hap, retaining all original-contig rows from both sources and adapting figure dimensions without pagination.
- Add local-window expansion counts, remaining unresolved reasons and report version to HTML; export the same evidence summary as window_qc.tsv.
- Restore the README workflow diagram and reproducible synthetic figure gallery in English and Chinese.

## 0.5.0

- Default to exact equal-copy assignment within homolog groups; require even ploidy and grouped chromosome order or a group table.
- Maximize total source-covered bases and retain independent labels, constrained labels, group scores, and tie diagnostics.
- Expand unresolved local windows up to their original gray-run extent with a complete evidence trace.
- Draw one continuous correspondence figure per recovered source, with height adjusted for label count.
- Add faint evidence-row guides and plain decimal assembly-length labels.
- Preserve independent source tracing with --assignment-mode provenance.

## 0.4.1

- Keep recovered sources on separate chromosome rows and correspondence pages at every source count.
- Add HTML source counts, assembly lengths, numbered rows and live filtered counts.
- Number chromosome and correspondence labels while retaining cross-source alignment evidence.
- Add hapreturn-report to redraw existing results without rerunning minimap2.

## 0.4.0

- Replace the four-panel overview with two publication figures in partition.pdf.
- Export independent alignment-matrix and original-to-curated correspondence PDFs and one complete HTML report.
- Add length-scaled chromosome ideograms, target-normalized contig matrices, oriented ribbons and their source tables.
- Add the hapreturn-chromosomes command and standalone chromosome script.
- Use bounded publication page sizes, embedded editable fonts and automatic pagination.

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
