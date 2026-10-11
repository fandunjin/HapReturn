#!/usr/bin/env python3
"""Split a curated merged assembly by its two hifiasm source assemblies.

This labels chromosome-scale sequences by assembly provenance. It does not
infer subgenomes or phase additional homologues within either source.
"""

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from . import __version__


def fasta_records(path):
    name, chunks = None, []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                header = line[1:].split()
                if not header:
                    raise ValueError("Empty FASTA identifier")
                name = header[0]
                chunks = []
            elif line.strip():
                if name is None:
                    raise ValueError("Sequence precedes FASTA header")
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def load_fasta(path):
    records = {}
    for name, seq in fasta_records(path):
        if name in records or not seq:
            raise ValueError(f"Duplicate identifier or empty sequence: {name}")
        records[name] = seq
    if not records:
        raise ValueError(f"No FASTA records in {path}")
    return records


def union_length(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    start, end = intervals[0]
    total = 0
    for a, b in intervals[1:]:
        if a > end:
            total += end - start
            start, end = a, b
        else:
            end = max(end, b)
    return total + end - start


def paf_intervals(path, records, min_mapq, min_aln_len):
    result = defaultdict(list)
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF {path}:{number}")
            try:
                target = fields[5]
                target_len, start, end = map(int, fields[6:9])
                aln_len, mapq = int(fields[10]), int(fields[11])
            except ValueError as exc:
                raise ValueError(f"Invalid PAF coordinates {path}:{number}") from exc
            if target not in records:
                raise ValueError(f"PAF target {target} absent from final FASTA")
            if target_len != len(records[target]) or not 0 <= start < end <= target_len:
                raise ValueError(f"PAF target coordinates/length disagree at {path}:{number}")
            # Supplementary chains may cover distinct pieces of a scaffold;
            # secondary hits duplicate loci and must not vote again.
            if "tp:A:S" in fields[12:]:
                continue
            if mapq >= min_mapq and aln_len >= min_aln_len:
                result[target].append((start, end))
    return result


def write_fasta(path, records, names):
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for name in names:
            handle.write(f">{name}\n")
            seq = records[name]
            for start in range(0, len(seq), 80):
                handle.write(seq[start:start + 80] + "\n")


def load_groups(path, records, ploidy):
    groups, seen = defaultdict(list), set()
    if not path:
        return groups
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"group_id", "seq_id"} <= set(reader.fieldnames):
            raise ValueError("Group TSV needs group_id and seq_id columns")
        for row in reader:
            group, name = row["group_id"], row["seq_id"]
            if not group or name not in records or name in seen:
                raise ValueError(f"Invalid/duplicate group member: {group}, {name}")
            groups[group].append(name)
            seen.add(name)
    if seen != set(records):
        raise ValueError("Group TSV must list every final FASTA record exactly once")
    if ploidy is not None:
        for group, names in groups.items():
            if len(names) != ploidy:
                raise ValueError(f"Group {group} has {len(names)} records; expected ploidy {ploidy}")
    return groups


def run_minimap(final, source, output, threads, preset):
    temporary = Path(str(output) + ".partial")
    try:
        with open(temporary, "w", encoding="utf-8") as handle:
            subprocess.run(
                ["minimap2", "-x", preset, "--secondary=no", "-t", str(threads),
                 str(final), str(source)], stdout=handle, check=True,
            )
        temporary.replace(output)
    except (OSError, subprocess.CalledProcessError):
        temporary.unlink(missing_ok=True)
        raise


def parse_args():
    parser = argparse.ArgumentParser(prog="hapreturn", description=__doc__)
    parser.add_argument("--version", action="version", version=f"HapReturn {__version__}")
    parser.add_argument("--final-fasta", type=Path, required=True)
    parser.add_argument("--hap1-fasta", type=Path)
    parser.add_argument("--hap2-fasta", type=Path)
    parser.add_argument("--hap1-paf", type=Path)
    parser.add_argument("--hap2-paf", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--chromosome-count", type=int, help="Select first N FASTA records; required if unplaced records follow")
    parser.add_argument("--chromosome-list", type=Path, help="One chromosome FASTA ID per line; safer when record order is arbitrary")
    parser.add_argument("--ploidy", type=int, help="Optional size check for each homolog group")
    parser.add_argument("--homolog-groups", type=Path, help="TSV: group_id, seq_id; used for QC only")
    parser.add_argument("--hap1-per-group", type=int, help="Expected hap1 records per group; used for QC only")
    parser.add_argument("--min-mapq", type=int, default=10)
    parser.add_argument("--min-aln-len", type=int, default=10000)
    parser.add_argument("--min-covered-fraction", type=float, default=0.25)
    parser.add_argument("--dominance", type=float, default=0.60)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--preset", default="asm5")
    parser.add_argument("--allow-unresolved", action="store_true")
    parser.add_argument("--force", action="store_true", help="Allow reuse of a nonempty output directory")
    return parser.parse_args()


def run():
    args = parse_args()
    if not (0.5 < args.dominance <= 1 and 0 <= args.min_covered_fraction <= 1):
        raise ValueError("dominance must be in (0.5,1]; min-covered-fraction in [0,1]")
    if args.ploidy is not None and args.ploidy < 1:
        raise ValueError("ploidy must be positive")
    if args.hap1_per_group is not None and not args.homolog_groups:
        raise ValueError("--hap1-per-group requires --homolog-groups")
    if args.hap1_per_group is not None and args.ploidy is not None and not 0 <= args.hap1_per_group <= args.ploidy:
        raise ValueError("hap1-per-group must lie between 0 and ploidy")
    if args.threads < 1 or args.min_mapq < 0 or args.min_aln_len < 1:
        raise ValueError("threads/min-mapq/min-aln-len out of range")
    if args.hap1_per_group is not None and args.hap1_per_group < 0:
        raise ValueError("hap1-per-group must be nonnegative")
    for label in ("hap1", "hap2"):
        paf, source = getattr(args, f"{label}_paf"), getattr(args, f"{label}_fasta")
        if (paf is None) == (source is None):
            raise ValueError(f"Provide exactly one of --{label}-paf or --{label}-fasta")
        path = paf if paf is not None else source
        if not path.is_file():
            raise ValueError(f"Input file does not exist: {path}")
    if args.out_dir.exists() and not args.force and any(args.out_dir.iterdir()):
        raise ValueError("Output directory is nonempty; use a new directory or --force")
    all_records = load_fasta(args.final_fasta)
    names = list(all_records)
    if args.chromosome_count is not None and args.chromosome_list is not None:
        raise ValueError("Use chromosome-count or chromosome-list, not both")
    if args.chromosome_list is not None:
        names = [line.strip() for line in args.chromosome_list.read_text(encoding="utf-8").splitlines()
                 if line.strip() and not line.startswith("#")]
        if not names or len(names) != len(set(names)) or set(names) - set(all_records):
            raise ValueError("Chromosome list is empty, duplicated, or names absent from FASTA")
    if args.chromosome_count is not None:
        if not 1 <= args.chromosome_count <= len(names):
            raise ValueError("chromosome-count exceeds FASTA record count")
        names = names[:args.chromosome_count]
    records = {name: all_records[name] for name in names}
    unplaced = [name for name in all_records if name not in records]
    if args.ploidy is not None and len(records) % args.ploidy:
        raise ValueError("Selected chromosome count is not divisible by ploidy")
    groups = load_groups(args.homolog_groups, records, args.ploidy)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    paf_paths = {}
    for label in ("hap1", "hap2"):
        paf = getattr(args, f"{label}_paf")
        source = getattr(args, f"{label}_fasta")
        if paf is None:
            if source is None:
                raise ValueError(f"Provide --{label}-paf or --{label}-fasta")
            paf = args.out_dir / f"{label}_to_final.paf"
            run_minimap(args.final_fasta, source, paf, args.threads, args.preset)
        paf_paths[label] = paf
    intervals = {label: paf_intervals(path, all_records, args.min_mapq, args.min_aln_len)
                 for label, path in paf_paths.items()}
    rows, assigned = [], {"hap1": [], "hap2": [], "unresolved": []}
    for name in names:
        length = len(records[name])
        a, b = (union_length(intervals[label][name]) for label in ("hap1", "hap2"))
        total = a + b
        dominant = max(a, b) / total if total else 0.0
        covered = max(a, b) / length
        label = ("hap1" if a > b else "hap2") if a != b else "unresolved"
        reason = "assigned"
        if covered < args.min_covered_fraction:
            reason, label = "low_coverage", "unresolved"
        elif dominant < args.dominance or a == b:
            reason, label = "ambiguous_origin", "unresolved"
        assigned[label].append(name)
        rows.append((name, length, a, b, f"{covered:.6f}", f"{dominant:.6f}", label, reason))
    with open(args.out_dir / "assignment.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("seq_id", "length", "hap1_covered_bp", "hap2_covered_bp",
                         "dominant_coverage_fraction", "dominance", "assignment", "reason"))
        writer.writerows(rows)
    qc_errors = []
    with open(args.out_dir / "group_qc.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("group_id", "n_records", "hap1", "hap2", "unresolved", "status"))
        for group, members in groups.items():
            counts = {label: len(set(members) & set(assigned[label])) for label in assigned}
            ok = counts["unresolved"] == 0 and (args.hap1_per_group is None or counts["hap1"] == args.hap1_per_group)
            writer.writerow((group, len(members), counts["hap1"], counts["hap2"], counts["unresolved"], "pass" if ok else "review"))
            if not ok:
                qc_errors.append(group)
    for label in assigned:
        write_fasta(args.out_dir / f"{label}.fa", records, assigned[label])
    write_fasta(args.out_dir / "unplaced.fa", all_records, unplaced)
    print(f"hap1={len(assigned['hap1'])} hap2={len(assigned['hap2'])} unresolved={len(assigned['unresolved'])} groups_to_review={len(qc_errors)}")
    if (assigned["unresolved"] or qc_errors) and not args.allow_unresolved:
        print("Review assignment.tsv and group_qc.tsv; use --allow-unresolved to accept partial output", file=sys.stderr)
        return 2
    return 0


def main():
    try:
        return run()
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
