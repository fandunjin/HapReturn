"""Source evidence plots and a portable, interactive HTML report."""
import csv
import html
import math
from pathlib import Path


def window_support(length, intervals, size):
    """Count union-covered target bases in half-open windows."""
    result = [0] * math.ceil(length / size)
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(end, merged[-1][1])
        else:
            merged.append([start, end])
    for start, end in merged:
        for index in range(start // size, (end - 1) // size + 1):
            result[index] += min(end, (index + 1) * size) - max(start, index * size)
    return result


def make_report(out_dir, labels, records, intervals, assignments, all_records,
                min_coverage=0.25, dominance=0.6, window_size=0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.ticker import MaxNLocator

    labels = list(labels)
    colors = {labels[0]: "#2878B5", labels[1]: "#E58B35",
              "unresolved": "#949BA5", "unplaced": "#D3D7DD"}
    size = window_size or max(1, math.ceil(max(records.values()) / 100))
    if sum(math.ceil(length / size) for length in records.values()) > 200000:
        raise ValueError("Too many plotting windows; increase --plot-window-size")
    windows, coverage = [], {}
    for name, length in records.items():
        support = [window_support(length, intervals[label].get(name, []), size) for label in labels]
        coverage[name] = [sum(values) / length for values in support]
        for index, (first, second) in enumerate(zip(*support)):
            start, end = index * size, min(length, (index + 1) * size)
            total, maximum = first + second, max(first, second)
            share = maximum / total if total else 0
            origin = labels[0] if first > second else labels[1]
            if not total or first == second or maximum / (end - start) < min_coverage or share < dominance:
                origin = "unresolved"
            windows.append((name, start, end, first, second, origin))
    with open(out_dir / "window_provenance.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["seq_id", "start", "end", *(f"{label}_covered_bp" for label in labels), "window_assignment"])
        writer.writerows(windows)

    names = list(records)
    categories = labels + ["unresolved", "unplaced"]
    members = {label: [name for name in all_records if assignments.get(name, "unplaced") == label] for label in categories}
    counts = [len(members[label]) for label in categories]
    bases = [sum(all_records[name] for name in members[label]) for label in categories]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none", "pdf.fonttype": 42})

    def draw(ax, kind):
        if kind == "coverage":
            for offset, label in zip([-0.18, 0.18], labels):
                index = labels.index(label)
                ax.barh([i + offset for i in range(len(names))], [coverage[n][index] for n in names],
                        height=0.32, color=colors[label], label=label)
            ax.set_yticks(range(len(names)), names)
            ax.invert_yaxis()
            ax.set_xlim(0, 1.03)
            ax.set_xlabel("Fraction of chromosome covered")
            ax.legend(frameon=False)
            ax.set_title("Chromosome source coverage", loc="left", weight="bold")
        elif kind == "assignment":
            for label in categories[:3]:
                selected = [n for n in names if assignments[n] == label]
                ax.scatter([coverage[n][0] for n in selected], [coverage[n][1] for n in selected],
                           s=[25 + 130 * records[n] / max(records.values()) for n in selected],
                           color=colors[label], alpha=0.8, edgecolor="white", label=label)
            ax.plot([0, 1], [0, 1], color="#ADB5BD", ls="--", lw=1)
            ratio = (1 - dominance) / dominance
            ax.plot([0, 1], [0, ratio], color="#ADB5BD", ls=":", lw=1)
            ax.plot([0, ratio], [0, 1], color="#ADB5BD", ls=":", lw=1)
            ax.axvline(min_coverage, color="#D3D7DD", lw=0.7)
            ax.axhline(min_coverage, color="#D3D7DD", lw=0.7)
            ax.set(xlim=(-0.03, 1.05), ylim=(-0.03, 1.05), xlabel=f"{labels[0]} coverage fraction", ylabel=f"{labels[1]} coverage fraction")
            ax.legend(frameon=False, fontsize=8)
            ax.set_title("Source assignment · size = length", loc="left", weight="bold")
        elif kind == "tracks":
            positions = {name: i for i, name in enumerate(names)}
            for name, start, end, first, second, origin in windows:
                ax.broken_barh([(start / 1e6, (end - start) / 1e6)], (positions[name] - 0.32, 0.64), facecolors=colors[origin], linewidth=0)
            ax.set_yticks(range(len(names)), names)
            ax.invert_yaxis()
            ax.set_xlabel("Chromosome position (Mb)")
            ax.set_xlim(0, max(records.values()) / 1e6)
            ax.legend(handles=[Patch(color=colors[l], label=l) for l in categories[:3]], frameon=False, fontsize=8)
            ax.set_title("Window source evidence", loc="left", weight="bold")
        else:
            x = list(range(4))
            ax.bar([i - 0.2 for i in x], counts, width=0.36, color=[colors[l] for l in categories])
            twin = ax.twinx()
            twin.spines["top"].set_visible(False)
            twin.bar([i + 0.2 for i in x], [b / 1e6 for b in bases], width=0.36,
                     color=[colors[l] for l in categories], alpha=0.4, hatch="//")
            ax.set_xticks(x, categories, rotation=15)
            ax.set_ylabel("Sequence count · solid")
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            twin.set_ylabel("Total Mb · hatched")
            ax.set_title("Recovered output sets", loc="left", weight="bold")

    height = max(4, len(names) * 0.24 + 1.8)
    figures = []
    for kind in ("coverage", "assignment", "tracks", "summary"):
        fig, ax = plt.subplots(figsize=(9, height if kind in ("coverage", "tracks") else 5), layout="constrained")
        draw(ax, kind)
        figures.append((kind, fig))
    fig, axes = plt.subplots(2, 2, figsize=(16, max(10, height * 2)), layout="constrained")
    for ax, kind in zip(axes.flat, ("coverage", "assignment", "tracks", "summary")):
        draw(ax, kind)
    fig.suptitle("HapReturn | two-source recovery", fontsize=18, weight="bold")
    figures.append(("overview", fig))
    for name, fig in figures:
        for extension in ("png", "svg", "pdf"):
            fig.savefig(out_dir / f"{name}.{extension}", dpi=160)
        plt.close(fig)

    # SVG is embedded directly: no external resources or JavaScript dependencies.
    panels = []
    for name in ("coverage", "assignment", "tracks", "summary"):
        svg = (out_dir / f"{name}.svg").read_text(encoding="utf-8")
        panels.append(svg[svg.index("<svg"):])
    escape = html.escape
    table = "".join(f'<tr data-source="{escape(assignments[n], quote=True)}"><td>{escape(n)}</td><td>{records[n]:,}</td><td>{coverage[n][0]:.3f}</td><td>{coverage[n][1]:.3f}</td><td>{escape(assignments[n])}</td></tr>' for n in names)
    tracks = []
    for n in names:
        rectangles = []
        for name, start, end, first, second, origin in windows:
            if name == n:
                tip = escape(f"{n} [{start:,}, {end:,}) | {labels[0]}: {first:,} bp | {labels[1]}: {second:,} bp | {origin}")
                rectangles.append(f'<rect x="{start / records[n] * 1000}" y="0" width="{(end-start) / records[n] * 1000}" height="24" fill="{colors[origin]}"><title>{tip}</title></rect>')
        tracks.append(f'<div class="track"><span>{escape(n)}</span><svg viewBox="0 0 1000 24" role="img" aria-label="{escape(n, quote=True)} source windows">{"".join(rectangles)}</svg><small>{records[n]:,} bp</small></div>')
    options = ''.join(f'<option value="{escape(l, quote=True)}">{escape(l)}</option>' for l in categories[:3])
    document = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>HapReturn report</title>
<style>body{font:15px system-ui,sans-serif;color:#253247;background:#f4f6f9;margin:0}main{max-width:1400px;margin:auto;padding:32px}h1{margin-bottom:8px}section{background:white;padding:24px;margin:20px 0;border-radius:12px}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(400px,1fr));gap:16px}.grid svg{width:100%;height:auto}.track{display:grid;grid-template-columns:150px 1fr 110px;align-items:center;gap:10px;margin:8px 0}.track svg{width:100%;height:24px}table{border-collapse:collapse;width:100%}td,th{text-align:left;padding:9px;border-bottom:1px solid #e4e8ee}tr:hover{background:#eef4fc}input,select{padding:8px;margin:8px}small{color:#637083}@media(max-width:600px){main{padding:12px}.grid{grid-template-columns:1fr}.track{grid-template-columns:90px 1fr}.track small{display:none}}</style>
<main><h1>HapReturn</h1><p>Two-source recovery · alignment evidence and output summary</p><section><div class="grid">''' + ''.join(panels) + '''</div></section><section><h2>Explore chromosome windows</h2><p>Hover over a window for coordinates and union-covered bases. Each row spans its own chromosome length. Gray indicates unresolved evidence.</p>''' + ''.join(tracks) + '''</section><section><h2>Chromosome assignments</h2><input id="search" placeholder="Search chromosome" aria-label="Search chromosome"><select id="source" aria-label="Filter source"><option value="">All sources</option>''' + options + f'''</select><table><thead><tr><th>Chromosome</th><th>Length (bp)</th><th>{escape(labels[0])} coverage</th><th>{escape(labels[1])} coverage</th><th>Assignment</th></tr></thead><tbody>{table}</tbody></table></section><p>Window size: {size:,} bp. Coordinates are zero-based, half-open. Window labels describe local source evidence; chromosome assignments use whole-chromosome coverage.</p></main>''' + '''<script>const search=document.getElementById('search'),source=document.getElementById('source');function filter(){document.querySelectorAll('tbody tr').forEach(row=>{row.hidden=!(row.cells[0].textContent.toLowerCase().includes(search.value.toLowerCase())&&(!source.value||row.dataset.source===source.value));});}search.addEventListener('input',filter);source.addEventListener('change',filter);</script></html>'''
    (out_dir / "report.html").write_text(document, encoding="utf-8")
