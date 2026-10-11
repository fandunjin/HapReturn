#!/usr/bin/env python3
"""Return a merged 2n assembly to two n source sets using hap1/hap2 provenance."""

import argparse
import csv
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from contextlib import ExitStack
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


def fasta_lengths(path):
    """Read lengths without retaining chromosome sequences in memory."""
    lengths = {}
    name = None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                fields = line[1:].split()
                if not fields:
                    raise ValueError("Empty FASTA identifier")
                name = fields[0]
                if name in lengths:
                    raise ValueError(f"Duplicate FASTA identifier: {name}")
                lengths[name] = 0
            elif line.strip():
                if name is None:
                    raise ValueError("Sequence precedes FASTA header")
                lengths[name] += len(line.strip())
    if not lengths or any(length == 0 for length in lengths.values()):
        raise ValueError("FASTA has no records or contains an empty sequence")
    return lengths


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
    from .evidence import read_paf
    result = defaultdict(list)
    for hit in read_paf(path, records, min_mapq, min_aln_len):
        result[hit.target].append((hit.start, hit.end))
    return result


def write_fasta(path, records, names):
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for name in names:
            handle.write(f">{name}\n")
            seq = records[name]
            for start in range(0, len(seq), 80):
                handle.write(seq[start:start + 80] + "\n")


def write_partitioned_fasta(path, out_dir, labels, assignments):
    """Stream FASTA lines directly to the selected source output."""
    with ExitStack() as stack:
        handles = {label: stack.enter_context(open(out_dir / f"{label}.fa", "w", encoding="utf-8", newline="\n"))
                   for label in list(labels) + ["unresolved", "unplaced"]}
        current = None
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    current = handles[assignments.get(name, "unplaced")]
                    current.write(f">{name}\n")
                elif line.strip():
                    current.write(line.strip() + "\n")


def source_config(args):
    sources = {}

    def add(label, path, kind):
        reserved = {"unresolved", "unplaced", "con", "prn", "aux", "nul", "group_id", "n_records", "status"}
        reserved.update(f"{prefix}{number}" for prefix in ("com", "lpt") for number in range(1, 10))
        if (not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*", label)
                or label.casefold() in reserved
                or label.casefold() in {name.casefold() for name in sources}):
            raise ValueError(f"Invalid, reserved or duplicate source label: {label}")
        if kind not in ("paf", "fasta") or not path.is_file():
            raise ValueError(f"Invalid source format or missing input: {label}, {path}, {kind}")
        sources[label] = (kind, path)

    legacy = any(getattr(args, f"hap{i}_{kind}") is not None for i in (1, 2) for kind in ("paf", "fasta"))
    if legacy and (args.sources or args.source_paf or args.source_fasta):
        raise ValueError("Use hap1/hap2 options or named sources, not both")
    if legacy:
        for label in ("hap1", "hap2"):
            paf, fasta = getattr(args, f"{label}_paf"), getattr(args, f"{label}_fasta")
            if (paf is None) == (fasta is None):
                raise ValueError(f"Provide exactly one of --{label}-paf or --{label}-fasta")
            add(label, paf if paf is not None else fasta, "paf" if paf is not None else "fasta")
    if args.sources:
        with open(args.sources, encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or not {"source", "path", "format"} <= set(reader.fieldnames):
                raise ValueError("Sources TSV requires source, path and format columns")
            for row in reader:
                path = Path(row["path"])
                if not path.is_absolute():
                    path = args.sources.parent / path
                add(row["source"], path, row["format"])
    for kind in ("paf", "fasta"):
        for spec in getattr(args, f"source_{kind}"):
            label, separator, path = spec.partition("=")
            if not separator or not path:
                raise ValueError(f"Expected LABEL=PATH for --source-{kind}")
            add(label, Path(path), kind)
    if len(sources) != 2:
        raise ValueError("Exactly two original sources are required for the 2n-to-n+n split")
    return sources


def source_quotas(args, sources):
    quotas = {}
    for spec in args.source_copies:
        label, separator, count = spec.partition("=")
        if not separator or label not in sources or label in quotas:
            raise ValueError(f"Invalid or duplicate source quota: {spec}")
        quotas[label] = int(count)
    if args.hap1_per_group is not None:
        if "hap1" not in sources or "hap1" in quotas:
            raise ValueError("--hap1-per-group requires source hap1 without a duplicate quota")
        quotas["hap1"] = args.hap1_per_group
    if quotas and not args.homolog_groups and args.assignment_mode == "provenance":
        raise ValueError("Source copy quotas require --homolog-groups")
    if any(count < 0 or (args.ploidy is not None and count > args.ploidy) for count in quotas.values()):
        raise ValueError("Source copy quota must be between zero and ploidy")
    if args.ploidy is not None and (sum(quotas.values()) > args.ploidy or
            (len(quotas) == len(sources) and sum(quotas.values()) != args.ploidy)):
        raise ValueError("Source copy quotas disagree with ploidy")
    return quotas


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
    parser.add_argument("--sources", type=Path, help="TSV for exactly two original sources: source, path, format (paf/fasta)")
    parser.add_argument("--source-paf", action="append", default=[], metavar="LABEL=PATH", help="Named PAF input; supply exactly two sources")
    parser.add_argument("--source-fasta", action="append", default=[], metavar="LABEL=PATH", help="Named FASTA input; supply exactly two sources")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--chromosome-count", type=int, help="Select first N FASTA records; required if unplaced records follow")
    parser.add_argument("--chromosome-list", type=Path, help="One chromosome FASTA ID per line; safer when record order is arbitrary")
    parser.add_argument("--assignment-mode", choices=("balanced", "provenance"), default="balanced",
                        help="balanced: equal-copy group optimization (default); provenance: independent source tracing")
    parser.add_argument("--ploidy", type=int, help="Even number of chromosomes per homolog group; required in balanced mode")
    parser.add_argument("--homolog-groups", type=Path, help="TSV: group_id, seq_id; otherwise consecutive FASTA records form groups")
    parser.add_argument("--hap1-per-group", type=int, help="Expected hap1 records per group; used for QC only")
    parser.add_argument("--source-copies", action="append", default=[], metavar="LABEL=N", help="Expected chromosome copies per homolog group for a source; repeat as needed")
    parser.add_argument("--min-mapq", type=int, default=10)
    parser.add_argument("--min-aln-len", type=int, default=10000)
    parser.add_argument("--min-covered-fraction", type=float, default=0.25)
    parser.add_argument("--dominance", type=float, default=0.60)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--preset", default="asm5")
    parser.add_argument("--allow-unresolved", action="store_true")
    parser.add_argument("--force", action="store_true", help="Allow reuse of a nonempty output directory")
    parser.add_argument("--no-plots", action="store_true", help="Skip visualization outputs")
    parser.add_argument("--plot-window-size", type=int, default=0, help="Window size in bp; 0 selects about 100 windows on the longest chromosome")
    return parser.parse_args()


def run():
    args = parse_args()
    if not (0.5 < args.dominance <= 1 and 0 <= args.min_covered_fraction <= 1):
        raise ValueError("dominance must be in (0.5,1]; min-covered-fraction in [0,1]")
    if args.ploidy is not None and args.ploidy < 1:
        raise ValueError("ploidy must be positive")
    if args.threads < 1 or args.min_mapq < 0 or args.min_aln_len < 1:
        raise ValueError("threads/min-mapq/min-aln-len out of range")
    if args.plot_window_size < 0:
        raise ValueError("plot-window-size must be nonnegative")
    sources = source_config(args)
    quotas = source_quotas(args, sources)
    if args.assignment_mode == "balanced":
        if args.ploidy is None or args.ploidy < 2 or args.ploidy % 2:
            raise ValueError("Balanced assignment requires an even --ploidy >=2; use --assignment-mode provenance for independent tracing")
        if any(count != args.ploidy // 2 for count in quotas.values()):
            raise ValueError("Balanced assignment requires equal source copy counts")
        quotas = {label: args.ploidy // 2 for label in sources}
    if args.out_dir.exists() and not args.force and any(args.out_dir.iterdir()):
        raise ValueError("Output directory is nonempty; use a new directory or --force")
    all_records = fasta_lengths(args.final_fasta)
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
    if args.ploidy is not None and len(records) % args.ploidy:
        raise ValueError("Selected chromosome count is not divisible by ploidy")
    groups = load_groups(args.homolog_groups, records, args.ploidy)
    if args.assignment_mode == "balanced" and not groups:
        from .assignment import consecutive_groups
        groups = consecutive_groups(records, args.ploidy)
    input_paths = {args.final_fasta.resolve(), *(path.resolve() for kind, path in sources.values())}
    input_paths.update(path.resolve() for path in (args.sources, args.chromosome_list, args.homolog_groups) if path is not None)
    output_names = [f"{label}.fa" for label in list(sources) + ["unresolved", "unplaced"]]
    output_names += [f"{label}_to_final.paf" for label in sources]
    output_names += ["assignment.tsv", "group_qc.tsv", "summary.tsv", "homolog_groups.tsv"]
    if not args.no_plots:
        output_names += ["report.html", "window_provenance.tsv", "window_rescue.tsv", "window_qc.tsv", "alignment_matrix.tsv", "alignment_links.tsv",
                         "partition.pdf", "alignment_matrix.pdf", "hap_comparison.pdf"]
        size = args.plot_window_size or max(1, (max(records.values()) + 99) // 100)
        if sum((length + size - 1) // size for length in records.values()) > 200000:
            raise ValueError("Too many plotting windows; increase --plot-window-size")
    if any((args.out_dir / name).resolve() in input_paths for name in output_names):
        raise ValueError("An output would overwrite an input; choose a separate output directory")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    paf_paths = {}
    for label, (kind, path) in sources.items():
        paf = path
        if kind == "fasta":
            paf = args.out_dir / f"{label}_to_final.paf"
            run_minimap(args.final_fasta, path, paf, args.threads, args.preset)
        else:
            paf = args.out_dir / f"{label}_to_final.paf"
            shutil.copyfile(path, paf)
        paf_paths[label] = paf
    from .evidence import read_paf
    alignments = []
    intervals = {label: defaultdict(list) for label in sources}
    for label, path in paf_paths.items():
        for hit in read_paf(path, all_records, args.min_mapq, args.min_aln_len, label):
            intervals[hit.source][hit.target].append((hit.start, hit.end))
            if not args.no_plots:
                alignments.append(hit)
    rows, support_by_name, evidence, evidence_reason = [], {}, {}, {}
    for name in names:
        length = records[name]
        support = {label: union_length(intervals[label][name]) for label in sources}
        maximum = max(support.values())
        total = sum(support.values())
        dominant = maximum / total if total else 0.0
        covered = maximum / length
        winners = [label for label, bp in support.items() if bp == maximum]
        label = winners[0] if len(winners) == 1 else "unresolved"
        reason = "assigned"
        if covered < args.min_covered_fraction:
            reason, label = "low_coverage", "unresolved"
        elif dominant < args.dominance or len(winners) != 1:
            reason, label = "ambiguous_origin", "unresolved"
        support_by_name[name] = support
        evidence[name], evidence_reason[name] = label, reason
        rows.append((name, length, *support.values(), f"{covered:.6f}", f"{dominant:.6f}", label, reason))
    assignments, group_audit = dict(evidence), {}
    if args.assignment_mode == "balanced":
        from .assignment import balanced_assignment
        assignments, group_audit = balanced_assignment(groups, support_by_name, list(sources), args.ploidy)
    group_of = {name: group for group, members in groups.items() for name in members}
    assigned = {label: [] for label in list(sources) + ["unresolved"]}
    for name in names:
        assigned[assignments[name]].append(name)
    with open(args.out_dir / "assignment.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("seq_id", "length", *(f"{label}_covered_bp" for label in sources),
                         "dominant_coverage_fraction", "dominance", "assignment", "reason",
                         "evidence_assignment", "evidence_reason", "assignment_mode", "group_id",
                         "assigned_coverage_fraction", "assigned_source_share"))
        for row in rows:
            name, label = row[0], assignments[row[0]]
            selected = support_by_name[name].get(label, 0)
            total = sum(support_by_name[name].values())
            reason = ("group_constrained" if label != evidence[name] else evidence_reason[name])
            writer.writerow((*row[:6], label, reason, evidence[name], evidence_reason[name],
                             args.assignment_mode, group_of.get(name, ""),
                             f"{selected / records[name]:.6f}", f"{selected / total if total else 0:.6f}"))
    with open(args.out_dir / "homolog_groups.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("group_id", "seq_id"))
        writer.writerows((group, name) for group, members in groups.items() for name in members)
    qc_errors = []
    with open(args.out_dir / "group_qc.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("group_id", "n_records", *sources, "unresolved", "status", "review_reason",
                         "changed_records", "total_support_bp", "optimality_gap_bp"))
        assigned_sets = {label: set(names) for label, names in assigned.items()}
        for group, members in groups.items():
            counts = {label: len(set(members) & assigned_sets[label]) for label in assigned}
            ok = counts["unresolved"] == 0 and all(counts[label] == count for label, count in quotas.items())
            review = []
            audit = group_audit.get(group, {})
            changed = sum(assignments[n] != evidence[n] for n in members)
            if audit.get("tied_optimum"):
                review.append("tied_optimum")
            if args.assignment_mode == "balanced" and any(evidence[n] != assignments[n] for n in members):
                review.append("constraint_changed_evidence")
            if not ok:
                review.append("copy_count_or_unresolved")
            ok = ok and not review
            writer.writerow((group, len(members), *counts.values(), "pass" if ok else "review", ";".join(review),
                             changed, audit.get("total_support_bp", ""), audit.get("optimality_gap_bp", "")))
            if not ok:
                qc_errors.append(group)
    write_partitioned_fasta(args.final_fasta, args.out_dir, sources, assignments)
    with open(args.out_dir / "summary.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("source", "n_sequences", "total_bp"))
        for label in list(assigned) + ["unplaced"]:
            members = assigned[label] if label in assigned else [name for name in all_records if name not in records]
            writer.writerow((label, len(members), sum(all_records[name] for name in members)))
    if not args.no_plots:
        from .visualize import make_report
        make_report(args.out_dir, sources, records, intervals, assignments, all_records,
                    args.min_covered_fraction, args.dominance, args.plot_window_size, alignments)
        if args.force:
            for stem in ("coverage", "assignment", "tracks", "summary", "overview"):
                for extension in ("png", "svg", "pdf"):
                    legacy = args.out_dir / f"{stem}.{extension}"
                    if legacy.is_file() and legacy.resolve() not in input_paths:
                        legacy.unlink()
    print(" ".join(f"{label}={len(members)}" for label, members in assigned.items()) + f" groups_to_review={len(qc_errors)}")
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
