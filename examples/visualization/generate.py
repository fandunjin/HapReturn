"""Generate synthetic tetraploid inputs for the README figure gallery.

Usage: python examples/visualization/generate.py OUTPUT_DIRECTORY

Three consecutive groups of four chromosomes illustrate clear support, mixed
local provenance, missing local evidence, and a constrained 2+2 assignment.
Sequences and alignment coordinates are artificial, not biological results.
"""
import sys
from pathlib import Path

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
lengths = [1200000, 1120000, 1160000, 1080000,
           1000000, 940000, 980000, 900000,
           820000, 780000, 800000, 760000]
with (out / "final.fa").open("w") as handle:
    for index, length in enumerate(lengths, 1):
        handle.write(f">chr{index:02d}\n" + "A" * length + "\n")
    handle.write(">unplaced01\n" + "A" * 50000 + "\n")
for source in (1, 2):
    with (out / f"hap{source}.paf").open("w") as handle:
        for index, length in enumerate(lengths, 1):
            major_source = 1 if index in (1, 2, 5, 7, 9) else 2
            fractions = [(0.0, 0.94)] if source == major_source else [(0.35, 0.41)]
            if index == 2:
                fractions = [(0.0, 0.56), (0.70, 0.94)] if source == 1 else [(0.56, 0.70)]
            elif index == 7:
                fractions = [(0.0, 0.36), (0.49, 0.94)] if source == 1 else [(0.15, 0.20)]
            elif index == 10:
                # Third group has 1:3 independent source labels; the optimizer
                # assigns chr10 to hap1 to maximize support with two per source.
                fractions = [(0.0, 0.35)] if source == 1 else [(0.35, 1.0)]
            spans = [(int(length * start), int(length * end)) for start, end in fractions]
            for start, end in spans:
                handle.write(f"ctg{index:02d}\t{length}\t{start}\t{end}\t+\tchr{index:02d}\t{length}\t{start}\t{end}\t{end-start}\t{end-start}\t60\ttp:A:P\n")
