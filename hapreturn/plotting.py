"""Publication figures with shared scales, editable text and continuous sets."""
import io
import xml.etree.ElementTree as ET
from itertools import groupby

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle, Patch
from matplotlib.ticker import MaxNLocator, StrMethodFormatter

WIDTH = 183 / 25.4
STYLE = {"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
         "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
         "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
         "axes.spines.top": False, "axes.spines.right": False,
         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
         "pdf.fonttype": 42, "svg.fonttype": "none", "legend.frameon": False,
         "savefig.facecolor": "white"}


def palette(labels):
    return dict(zip(list(labels) + ["unresolved", "unplaced"],
                    ["#3274A1", "#D88A3D", "#A5ABB3", "#DCE0E5"]))


def panel(ax, letter, title):
    ax.set_title(title, loc="left", pad=8)
    ax.text(-0.10, 1.03, letter, transform=ax.transAxes, fontsize=8, weight="bold", va="bottom")


def save_page(fig, pdf, identifier):
    """Save vector PDF and return the same page as an embedded SVG."""
    try:
        pdf.savefig(fig)
        buffer = io.StringIO()
        fig.savefig(buffer, format="svg")
        root = ET.fromstring(buffer.getvalue())
        ids = {item.attrib["id"]: identifier + "_" + item.attrib["id"]
               for item in root.iter() if "id" in item.attrib}
        for item in root.iter():
            for key, value in list(item.attrib.items()):
                if key == "id":
                    item.attrib[key] = ids[value]
                else:
                    if value.startswith("#") and value[1:] in ids:
                        value = "#" + ids[value[1:]]
                    for old, new in ids.items() if "url(#" in value else []:
                        value = value.replace(f"url(#{old})", f"url(#{new})")
                    item.attrib[key] = value
        ET.register_namespace("", "http://www.w3.org/2000/svg")
        ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
        return ET.tostring(root, encoding="unicode")
    finally:
        plt.close(fig)


def capsule(ax, x, bottom, length, width, cap, facecolor="none", zorder=3):
    """Round-ended sequence outline; no unmeasured centromere is invented."""
    radius = min(cap, length * 0.2)
    left, right, top = x - width / 2, x + width / 2, bottom + length
    vertices = [(x, bottom), (left, bottom), (left, bottom), (left, bottom + radius),
                (left, top - radius), (left, top), (left, top), (x, top),
                (right, top), (right, top), (right, top - radius), (right, bottom + radius),
                (right, bottom), (right, bottom), (x, bottom), (x, bottom)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    patch = PathPatch(Path(vertices, codes), facecolor=facecolor, edgecolor="#3D4650", lw=0.45, zorder=zorder)
    ax.add_patch(patch)
    return patch


def source_blocks(data, size=16):
    """Chunk within recovered sources so a row or target page never mixes them."""
    blocks = []
    for label in data.labels + ["unresolved"]:
        members = [name for name in data.records if data.assignments[name] == label]
        for offset in range(0, len(members), size):
            blocks.append((label, members[offset:offset + size]))
    return blocks


def partition_pages(data, pdf):
    labels, names, colors = data.labels, list(data.records), palette(data.labels)
    categories = labels + ["unresolved", "unplaced"]
    pages = []
    # Figure 1 combines nonredundant evidence and output metrics; no dual axes.
    for offset in range(0, len(names), 64):
        selected = names[offset:offset + 64]
        fig = plt.figure(figsize=(WIDTH, min(170, max(110, 35 + len(selected) * 2.7)) / 25.4))
        grid = fig.add_gridspec(2, 2, width_ratios=[2.0, 1], hspace=0.65, wspace=0.6,
                               left=0.15, right=0.97, bottom=0.22, top=0.88)
        ax = fig.add_subplot(grid[:, 0])
        for index, (label, marker) in enumerate(zip(labels, ("o", "s"))):
            ax.scatter([data.coverage[n][index] for n in selected],
                       [i + (index - 0.5) * 0.22 for i in range(len(selected))],
                       s=13, marker=marker, color=colors[label], label=label, linewidths=0)
        ax.axvline(data.min_coverage, lw=0.6, ls="--", color="#767C85")
        ax.set_yticks(range(len(selected)), [f"{offset + i + 1}. {n}" for i, n in enumerate(selected)])
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#E3E7EC", linewidth=0.45)
        if len(selected) > 48:
            ax.tick_params(axis="y", labelsize=5)
        ax.set(xlim=(-0.02, 1.03), ylim=(len(selected) - 0.5, -0.6), xlabel="Chromosome coverage (fraction)")
        ax.tick_params(axis="y", length=0)
        panel(ax, "a", "Source evidence")
        ax.legend(loc="upper left", bbox_to_anchor=(0, -0.14), ncol=2, borderaxespad=0)
        for row, (metric, title, unit) in enumerate(((data.counts, "Recovered sequence counts", "Sequences (count)"),
                                                    ([x / 1e6 for x in data.bases], "Recovered assembly lengths", "Sequence length (Mb)"))):
            ax = fig.add_subplot(grid[row, 1])
            ax.barh(categories, metric, height=0.52, color=[colors[l] for l in categories])
            ax.invert_yaxis()
            ax.set_xlabel(unit)
            ax.xaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
            ax.tick_params(axis="y", length=0)
            if row == 0:
                ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=4))
            for i, value in enumerate(metric):
                ax.annotate(f"{value:,.0f}" if row == 0 else f"{value:,.2f}", (value, i),
                            xytext=(3, 0), textcoords="offset points", va="center", fontsize=6)
            ax.set_xlim(0, max(metric or [0]) * 1.3 or 1)
            panel(ax, "bc"[row], title)
        fig.text(0.06, 0.97, "Source support and recovered chromosome sets", fontsize=7, weight="bold", va="top")
        fig.text(0.06, 0.045, f"Selected sequences: {len(names)}. Coverage threshold: {data.min_coverage:g}; dominance: {data.dominance:g}.\n"
                 "Coverage is calculated independently for each source; output statistics include unselected records.", fontsize=6)
        pages.append(save_page(fig, pdf, f"partition_evidence_{offset}"))

    # Figure 2: source windows inside chromosome outlines, all at one bp scale.
    blocks = source_blocks(data)
    maximum = max(data.records.values()) / 1e6
    for offset in range(0, len(blocks), 4):
        selected_blocks = blocks[offset:offset + 4]
        nrows = len(selected_blocks)
        fig = plt.figure(figsize=(WIDTH, min(170, max(95, 40 + 31 * nrows)) / 25.4))
        grid = fig.add_gridspec(nrows, 1, left=0.10, right=0.98, top=0.82, bottom=0.22, hspace=0.72)
        for row in range(nrows):
            ax = fig.add_subplot(grid[row, 0])
            assignment, subset = selected_blocks[row]
            for x, name in enumerate(subset):
                length = data.records[name] / 1e6
                shape = capsule(ax, x, 0, length, 0.48, maximum * 0.045)
                runs = []
                for _, start, end, first, second, origin in data.by_chromosome[name]:
                    if runs and runs[-1][2] == origin and runs[-1][1] == start:
                        runs[-1][1] = end
                    else:
                        runs.append([start, end, origin])
                for start, end, origin in runs:
                    band = Rectangle((x - 0.24, start / 1e6), 0.48, (end - start) / 1e6,
                                     facecolor=colors[origin], edgecolor="none", antialiased=False, zorder=2)
                    band.set_clip_path(shape)
                    ax.add_patch(band)
            for assignment, group in groupby(enumerate(subset), key=lambda pair: data.assignments[pair[1]]):
                indices = [pair[0] for pair in group]
                start, end = indices[0] - 0.25, indices[-1] + 0.25
                ax.plot([start, end], [maximum * 1.08] * 2, lw=1.8, color=colors[assignment], solid_capstyle="butt")
                count = sum(data.assignments[n] == assignment for n in names)
                ax.text((start + end) / 2, maximum * 1.12, f"{assignment} (n = {count})", ha="center", va="bottom", fontsize=6)
            members = [n for n in names if data.assignments[n] == assignment]
            ax.set_xticks(range(len(subset)), [f"{members.index(n) + 1}. {n}" for n in subset], rotation=60, ha="right")
            ax.set(xlim=(-0.7, max(7.5, max(len(b[1]) for b in selected_blocks) - 0.3)), ylim=(0, maximum * 1.28), ylabel="Position (Mb)")
            ax.tick_params(axis="x", length=0, pad=2)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
            ax.spines["bottom"].set_visible(False)
        fig.text(0.05, 0.97, "Chromosome-scale source provenance", fontsize=7, weight="bold", va="top")
        fig.legend(handles=[Patch(facecolor=colors[l], label=l) for l in categories[:3]],
                   loc="upper left", bbox_to_anchor=(0.08, 0.93), ncol=3)
        fig.text(0.05, 0.025, f"Window: {data.window_size:,} bp. Group labels above chromosomes show whole-chromosome assignments.\n"
                 "Band colors include expanded-context labels; evidence spans are recorded in window_rescue.tsv.", fontsize=6)
        pages.append(save_page(fig, pdf, f"partition_chromosomes_{offset}"))
    return pages


def matrix_pages(data, matrix, pdf):
    import numpy as np
    labels, colors, names = data.labels, palette(data.labels), list(data.records)
    rowkeys = list(matrix)
    pages = []
    for offset, (recovered_source, columns) in enumerate(source_blocks(data, size=len(names))):
        # Keep every original-contig row from both sources, including zero cells.
        # Increase physical dimensions to retain readable labels without tiling.
        width_mm = max(183, 76 + 4 * len(columns))
        height_mm = max(100, 64 + 2.8 * len(rowkeys))
        values = np.array([[matrix[key].get(n, 0.0) for n in columns] for key in rowkeys])
        fig, ax = plt.subplots(figsize=(width_mm / 25.4, height_mm / 25.4))
        fig.subplots_adjust(left=54 / width_mm, right=1 - 22 / width_mm,
                            bottom=30 / height_mm, top=1 - 34 / height_mm)
        mesh = ax.pcolormesh(values, cmap="cividis", vmin=0, vmax=1, shading="flat", rasterized=False,
                            edgecolors="face", linewidth=0, antialiased=False)
        ax.set_xticks([i + 0.5 for i in range(len(columns))], columns, rotation=65, ha="right")
        ax.set_yticks([i + 0.5 for i in range(len(rowkeys))], [f"{s} | {q}" for s, q in rowkeys])
        ax.invert_yaxis()
        ax.tick_params(length=0, pad=3)
        ax.tick_params(axis="y", pad=16)
        ax.set(xlabel="Curated chromosome", ylabel="Original hifiasm source | contig")
        for i, (source, query) in enumerate(rowkeys):
            ax.add_patch(Rectangle((-len(columns) * 0.035, i), len(columns) * 0.018, 1,
                                   facecolor=colors[source], edgecolor="none", clip_on=False))
        for i, name in enumerate(columns):
            ax.add_patch(Rectangle((i, -0.35), 1, 0.20, facecolor=colors[data.assignments[name]], edgecolor="none", clip_on=False))
        colorbar_mm = min(50, (height_mm - 64) * 0.6)
        cax = fig.add_axes([1 - 15 / width_mm, 0.5 - colorbar_mm / height_mm / 2,
                            2 / width_mm, colorbar_mm / height_mm])
        bar = fig.colorbar(mesh, cax=cax, ticks=[0, 0.5, 1])
        bar.set_label("Target fraction covered")
        fig.text(8 / width_mm, 1 - 7 / height_mm,
                 f"Original-contig alignment matrix → recovered {recovered_source} (n = {len(columns)})",
                 fontsize=7, weight="bold", va="top")
        fig.legend(handles=[Patch(facecolor=colors[l], label=l) for l in labels + ["unresolved"]],
                   loc="upper left", bbox_to_anchor=(54 / width_mm, 1 - 16 / height_mm), ncol=3)
        fig.text(8 / width_mm, 4 / height_mm,
                 "Cells: union of target bases aligned from one original contig / target chromosome length.\n"
                 f"All {len(rowkeys)} original-contig rows from both sources and all {len(columns)} recovered {recovered_source} chromosomes are shown continuously.",
                 fontsize=6)
        pages.append(save_page(fig, pdf, f"matrix_{offset}"))
    return pages
