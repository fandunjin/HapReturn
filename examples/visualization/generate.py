"""Generate synthetic inputs illustrating clear, mixed and missing evidence.

Usage: python examples/visualization/generate.py OUTPUT_DIRECTORY
"""
import sys
from pathlib import Path

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
lengths = [100000, 90000, 80000, 95000, 85000, 75000, 100000, 60000]
with (out / "final.fa").open("w") as handle:
    for index, length in enumerate(lengths, 1):
        handle.write(f">chr{index:02d}\n" + "A" * length + "\n")
for source in (1, 2):
    with (out / f"hap{source}.paf").open("w") as handle:
        for index, length in enumerate(lengths, 1):
            spans = []
            if (source == 1 and index <= 3) or (source == 2 and 4 <= index <= 6):
                spans = [(0, int(length * (0.95 - index * 0.015)))]
            elif index == 7:
                spans = [(0, 50000)] if source == 1 else [(50000, length)]
            elif index <= 6:
                spans = [(length // 3, length // 3 + length // 12)]
            for start, end in spans:
                handle.write(f"query{index}\t{end-start}\t0\t{end-start}\t+\tchr{index:02d}\t{length}\t{start}\t{end}\t{end-start}\t{end-start}\t60\ttp:A:P\n")
