"""Draw coordinate-aware original hifiasm to recovered hap correspondence."""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


def write_links(path, alignments, assignments):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["original_source", "original_contig", "query_length", "query_start", "query_end",
                         "strand", "curated_chromosome", "target_length", "target_start", "target_end", "recovered_source"])
        for hit in alignments:
            if hit.target in assignments:
                writer.writerow([hit.source, hit.query, hit.query_length, hit.query_start, hit.query_end,
                                 hit.strand, hit.target, hit.target_length, hit.start, hit.end, assignments[hit.target]])


def stacked_positions(lengths, scale, minimum=0.022, gap=0.005):
    """Reserve label space separately from each length-scaled sequence body."""
    cursor, result = 1.0, {}
    for key, length in lengths.items():
        body = length * scale
        span = max(body, minimum)
        result[key] = (cursor - span / 2 - body / 2, body)
        cursor -= span + gap
    return result


def common_scale(left, right):
    lo, hi = 0, 1 / max(list(left.values()) + list(right.values()))
    for _ in range(60):
        middle = (lo + hi) / 2
        extent = max(sum(max(length * middle, 0.022) + 0.005 for length in side.values()) for side in (left, right))
        if extent > 0.95:
            hi = middle
        else:
            lo = middle
    return lo


def ribbon_vertices(hit, query_position, target_position):
    """Return oriented ribbon endpoints; minus-strand target ends are swapped."""
    qb, qh = query_position
    tb, th = target_position
    q0 = qb + qh * hit.query_start / hit.query_length
    q1 = qb + qh * hit.query_end / hit.query_length
    t0 = tb + th * hit.start / hit.target_length
    t1 = tb + th * hit.end / hit.target_length
    return (q0, q1, t0, t1) if hit.strand == "+" else (q0, q1, t1, t0)


def comparison_pages(data, alignments, pdf):
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import PathPatch, Patch
    from .plotting import WIDTH, palette, capsule, save_page
    colors, pages = palette(data.labels), []
    names = sorted(data.records, key=lambda n: (list(data.labels + ["unresolved"]).index(data.assignments[n]), list(data.records).index(n)))
    for offset in range(0, len(names), 16):
        targets = names[offset:offset + 16]
        hits = [h for h in alignments if h.target in targets]
        queries = {(h.source, h.query): h.query_length for h in hits}
        queries = dict(sorted(queries.items(), key=lambda item: (data.labels.index(item[0][0]), item[0][1])))
        keys = list(queries)
        for qoffset in range(0, max(1, len(keys)), 32):
            subset = keys[qoffset:qoffset + 32]
            left = {key: queries[key] for key in subset}
            right = {name: data.records[name] for name in targets}
            scale = common_scale(left, right)
            lpos, rpos = stacked_positions(left, scale), stacked_positions(right, scale)
            fig, ax = plt.subplots(figsize=(WIDTH, 160 / 25.4))
            fig.subplots_adjust(left=0.025, right=0.975, top=0.84, bottom=0.14)
            for hit in hits:
                key = (hit.source, hit.query)
                if key not in lpos:
                    continue
                q0, q1, t0, t1 = ribbon_vertices(hit, lpos[key], rpos[hit.target])
                verts = [(0.33, q0), (0.46, q0), (0.54, t0), (0.67, t0),
                         (0.67, t1), (0.54, t1), (0.46, q1), (0.33, q1), (0.33, q0)]
                codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
                         MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.CLOSEPOLY]
                ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=colors[hit.source], edgecolor="none", alpha=0.20, zorder=1))
            for (source, query), (bottom, height) in lpos.items():
                capsule(ax, 0.32, bottom, height, 0.02, 0.009, facecolor=colors[source])
                ax.text(0.295, bottom + height / 2, f"{source} | {query}", va="center", ha="right", fontsize=5)
            for target, (bottom, height) in rpos.items():
                capsule(ax, 0.68, bottom, height, 0.02, 0.009, facecolor=colors[data.assignments[target]])
                ax.text(0.705, bottom + height / 2, f"{target} | {data.assignments[target]}", va="center", fontsize=5)
            ax.text(0.32, 1.04, "Original hifiasm\naligned contigs", ha="center", fontsize=7)
            ax.text(0.68, 1.04, "Curated assembly\nrecovered chromosomes", ha="center", fontsize=7)
            if not hits:
                ax.text(0.5, 0.55, "No accepted alignments", ha="center", fontsize=7)
            ax.set(xlim=(0, 1), ylim=(-0.03, 1.1))
            ax.axis("off")
            # One shared bp scale applies to both sides; gaps only reserve label space.
            largest = max(list(left.values()) + list(right.values()))
            scale_bp = 10 ** max(0, len(str(largest)) - 1)
            ax.plot([0.05, 0.05], [0.04, 0.04 + scale_bp * scale], color="#333333", lw=1)
            ax.text(0.065, 0.04 + scale_bp * scale / 2, f"{scale_bp / 1e6:g} Mb", va="center", fontsize=6)
            fig.text(0.04, 0.97, "Original hap provenance and curated hap correspondence", fontsize=7, weight="bold", va="top")
            fig.legend(handles=[Patch(facecolor=colors[l], label=l) for l in data.labels + ["unresolved"]],
                       loc="upper left", bbox_to_anchor=(0.04, 0.94), ncol=3)
            fig.text(0.04, 0.055, "Ribbons retain query/target coordinates and strand; color denotes original source. Both sides share a bp scale.\n"
                     "Twisted ribbons denote reverse-strand links. Only accepted, aligned original contigs are shown.\n"
                     f"Target sequences {offset + 1}–{offset + len(targets)} of {len(names)}; original contigs {qoffset + 1 if keys else 0}–{qoffset + len(subset)} of {len(keys)} in this target block.", fontsize=6)
            pages.append(save_page(fig, pdf, f"comparison_{offset}_{qoffset}"))
    return pages


def main(argv=None):
    parser = argparse.ArgumentParser(description="Plot original hifiasm hap provenance against recovered chromosomes")
    parser.add_argument("--final-fasta", type=Path, required=True)
    parser.add_argument("--hap1-paf", type=Path, required=True)
    parser.add_argument("--hap2-paf", type=Path, required=True)
    parser.add_argument("--assignment", type=Path, required=True, help="HapReturn assignment.tsv")
    parser.add_argument("--source-labels", nargs=2, default=["hap1", "hap2"], metavar=("FIRST", "SECOND"))
    parser.add_argument("--min-mapq", type=int, default=10)
    parser.add_argument("--min-aln-len", type=int, default=10000)
    parser.add_argument("--out-pdf", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        from .cli import fasta_lengths
        from .evidence import read_paf
        from .visualize import prepare_data
        from .plotting import STYLE
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        if len(set(args.source_labels)) != 2 or set(args.source_labels) & {"unresolved", "unplaced"}:
            raise ValueError("Provide two distinct source labels")
        if args.min_mapq < 0 or args.min_aln_len < 1 or args.out_pdf.suffix.lower() != ".pdf":
            raise ValueError("Check mapping thresholds and the .pdf output suffix")
        inputs = [args.final_fasta, args.hap1_paf, args.hap2_paf, args.assignment]
        if args.out_pdf.resolve() in {p.resolve() for p in inputs}:
            raise ValueError("Output would overwrite an input")
        all_records = fasta_lengths(args.final_fasta)
        assignments = {}
        with args.assignment.open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                name, label = row["seq_id"], row["assignment"]
                if name in assignments or name not in all_records or label not in args.source_labels + ["unresolved"]:
                    raise ValueError("Invalid or duplicate chromosome assignment")
                assignments[name] = label
        if not assignments:
            raise ValueError("No chromosome assignments")
        records = {name: all_records[name] for name in assignments}
        alignments = [hit for label, path in zip(args.source_labels, (args.hap1_paf, args.hap2_paf))
                      for hit in read_paf(path, all_records, args.min_mapq, args.min_aln_len, label)]
        intervals = {l: defaultdict(list) for l in args.source_labels}
        for hit in alignments:
            intervals[hit.source][hit.target].append((hit.start, hit.end))
        data = prepare_data(args.source_labels, records, intervals, assignments, all_records)
        args.out_pdf.parent.mkdir(parents=True, exist_ok=True)
        with plt.rc_context(STYLE), PdfPages(args.out_pdf) as pdf:
            comparison_pages(data, alignments, pdf)
        print(f"Wrote {args.out_pdf.name}")
        return 0
    except (ValueError, OSError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
