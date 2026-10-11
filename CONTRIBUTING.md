# Contributing

Install locally with `python -m pip install -e .` and run
`python -m unittest discover -s tests -v`.

For bug reports, include the HapReturn and Python versions, the command, the
error text, and a small synthetic FASTA/PAF reproducer. Avoid uploading private
sequencing data. Contributions should preserve sequence identifiers and report
ambiguous evidence explicitly. Add regression tests for changes to assignment
or coordinate validation.
