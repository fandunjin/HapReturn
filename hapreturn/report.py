"""Regenerate figures from an existing result without rerunning alignments."""
import argparse
import csv
import shutil
import sys
from collections import defaultdict
from pathlib import Path

from .cli import fasta_lengths
from .evidence import read_paf
from .visualize import make_report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--min-mapq", type=int, default=10)
    parser.add_argument("--min-aln-len", type=int, default=10000)
    parser.add_argument("--min-covered-fraction", type=float, default=0.25)
    parser.add_argument("--dominance", type=float, default=0.60)
    parser.add_argument("--plot-window-size", type=int, default=0)
    args = parser.parse_args(argv)
    try:
        if args.min_mapq < 0 or args.min_aln_len < 1 or args.plot_window_size < 0:
            raise ValueError("Invalid alignment filter or window size")
        if not 0 <= args.min_covered_fraction <= 1 or not 0.5 <= args.dominance <= 1:
            raise ValueError("Invalid coverage or dominance threshold")
        with (args.result_dir / "summary.tsv").open(encoding="utf-8-sig", newline="") as handle:
            labels = [row["source"] for row in csv.DictReader(handle, delimiter="\t")
                      if row["source"] not in ("unresolved", "unplaced")]
        if len(labels) != 2 or len(set(labels)) != 2:
            raise ValueError("Expected exactly two source labels in summary.tsv")
        all_records = {}
        for label in labels + ["unresolved", "unplaced"]:
            path = args.result_dir / f"{label}.fa"
            records = fasta_lengths(path) if path.stat().st_size else {}
            if all_records.keys() & records.keys():
                raise ValueError("Duplicate sequence IDs across partition FASTAs")
            all_records.update(records)
        records, assignments = {}, {}
        with (args.result_dir / "assignment.tsv").open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                name, label, length = row["seq_id"], row["assignment"], int(row["length"])
                if name in records or all_records.get(name) != length or label not in labels + ["unresolved"]:
                    raise ValueError("Invalid chromosome assignment or FASTA length")
                records[name], assignments[name] = length, label
        if not records:
            raise ValueError("No selected chromosomes")
        intervals = {label: defaultdict(list) for label in labels}
        hits = []
        for label in labels:
            for hit in read_paf(args.result_dir / f"{label}_to_final.paf", all_records,
                                args.min_mapq, args.min_aln_len, label):
                hits.append(hit)
                intervals[label][hit.target].append((hit.start, hit.end))
        args.out_dir.mkdir(parents=True, exist_ok=True)
        for name in ("assignment.tsv", "summary.tsv", "group_qc.tsv"):
            source, target = args.result_dir / name, args.out_dir / name
            if source.is_file() and source.resolve() != target.resolve():
                shutil.copyfile(source, target)
        make_report(args.out_dir, labels, records, intervals, assignments, all_records,
                    args.min_covered_fraction, args.dominance, args.plot_window_size, hits)
        print(f"Refreshed report for {len(records)} chromosomes: {args.out_dir}")
        return 0
    except (OSError, ValueError, KeyError) as exc:
        print(f"hapreturn-report: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
