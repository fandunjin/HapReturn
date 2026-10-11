"""Source data, publication PDF figures and a portable HTML report."""
import csv
import html
import math
from collections import defaultdict
from dataclasses import dataclass


def window_support(length, intervals, size):
    """Count union-covered target bases in half-open windows."""
    if size < 1 or length < 1:
        raise ValueError("Window size and sequence length must be positive")
    result, merged = [0] * math.ceil(length / size), []
    for start, end in sorted(intervals):
        if not 0 <= start < end <= length:
            raise ValueError("Window interval outside sequence")
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(end, merged[-1][1])
        else:
            merged.append([start, end])
    for start, end in merged:
        for index in range(start // size, (end - 1) // size + 1):
            result[index] += min(end, (index + 1) * size) - max(start, index * size)
    return result


@dataclass
class ReportData:
    labels: list
    records: dict
    assignments: dict
    coverage: dict
    windows: list
    by_chromosome: dict
    counts: list
    bases: list
    min_coverage: float
    dominance: float
    window_size: int


def prepare_data(labels, records, intervals, assignments, all_records,
                 min_coverage=0.25, dominance=0.6, window_size=0):
    labels = list(labels)
    size = window_size or max(1, math.ceil(max(records.values()) / 100))
    if size < 1 or sum(math.ceil(length / size) for length in records.values()) > 200000:
        raise ValueError("Too many plotting windows; increase --plot-window-size")
    windows, coverage, by_chromosome = [], {}, defaultdict(list)
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
            row = (name, start, end, first, second, origin)
            windows.append(row)
            by_chromosome[name].append(row)
    categories = labels + ["unresolved", "unplaced"]
    members = {label: [n for n in all_records if assignments.get(n, "unplaced") == label] for label in categories}
    return ReportData(labels, records, assignments, coverage, windows, by_chromosome,
                      [len(members[l]) for l in categories],
                      [sum(all_records[n] for n in members[l]) for l in categories],
                      min_coverage, dominance, size)


def alignment_matrix(data, alignments, intervals):
    from .cli import union_length
    cells = defaultdict(lambda: defaultdict(list))
    for hit in alignments:
        if hit.target in data.records:
            cells[(hit.source, hit.query)][hit.target].append((hit.start, hit.end))
    if not cells:
        for label in data.labels:
            for name in data.records:
                cells[(label, "all accepted contigs")][name] = intervals[label].get(name, [])
    # Group sources; order contigs by their best supported chromosome, then ID.
    fractions = {key: {name: union_length(spans) / data.records[name] for name, spans in values.items()}
                 for key, values in cells.items()}
    names = list(data.records)
    def order(key):
        values = fractions[key]
        best = max(values, key=values.get) if values else names[0]
        return data.labels.index(key[0]), names.index(best), key[1]
    return {key: fractions[key] for key in sorted(fractions, key=order)}


def write_html(path, data, sections):
    from .plotting import palette
    escape, colors = html.escape, palette(data.labels)
    figures = "".join(f'<section><h2>{escape(title)}</h2>{"".join(pages)}</section>' for title, pages in sections)
    table = "".join(f'<tr data-source="{escape(data.assignments[n], quote=True)}"><td>{escape(n)}</td><td>{data.records[n]:,}</td><td>{data.coverage[n][0]:.3f}</td><td>{data.coverage[n][1]:.3f}</td><td>{escape(data.assignments[n])}</td></tr>' for n in data.records)
    maximum = max(data.records.values())
    tracks = []
    for index, (name, length) in enumerate(data.records.items()):
        width = length / maximum * 1000
        clip = f"window_track_{index}"
        pieces = [f'<defs><clipPath id="{clip}"><rect width="{width}" height="24" rx="12"/></clipPath></defs>']
        for _, start, end, first, second, origin in data.by_chromosome[name]:
            tip = escape(f"{name} [{start:,}, {end:,}) | {data.labels[0]}: {first:,} bp | {data.labels[1]}: {second:,} bp | {origin}")
            pieces.append(f'<rect clip-path="url(#{clip})" x="{start / maximum * 1000}" y="0" width="{(end-start) / maximum * 1000}" height="24" fill="{colors[origin]}"><title>{tip}</title></rect>')
        tracks.append(f'<div class="track"><span>{escape(name)}</span><svg viewBox="0 0 1000 24" role="img" aria-label="{escape(name, quote=True)} source windows">{"".join(pieces)}</svg><small>{length:,} bp</small></div>')
    options = ''.join(f'<option value="{escape(l, quote=True)}">{escape(l)}</option>' for l in data.labels + ["unresolved"])
    document = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>HapReturn report</title>
<style>body{font:15px Arial,sans-serif;color:#253247;background:#f5f6f8;margin:0}main{max-width:1100px;margin:auto;padding:28px}section{background:white;padding:24px;margin:20px 0}section>svg{display:block;width:100%;height:auto;max-height:960px;margin:24px auto}.track{display:grid;grid-template-columns:150px 1fr 110px;align-items:center;gap:10px;margin:8px 0}.track svg{width:100%;height:24px}table{border-collapse:collapse;width:100%}td,th{text-align:left;padding:9px;border-bottom:1px solid #e4e8ee}tr:hover{background:#eef4fc}input,select{padding:8px;margin:8px}small{color:#637083}@media(max-width:600px){main{padding:10px}.track{grid-template-columns:90px 1fr}.track small{display:none}}</style>
<main><h1>HapReturn</h1><p>Two-source chromosome recovery</p>''' + figures + '''<section><h2>Explore chromosome windows</h2><p>Hover for coordinates and union-covered bases. All rows share the same length scale; gray indicates unresolved evidence.</p>''' + ''.join(tracks) + '''</section><section><h2>Chromosome assignments</h2><input id="search" placeholder="Search chromosome" aria-label="Search chromosome"><select id="source" aria-label="Filter source"><option value="">All sources</option>''' + options + f'''</select><table><thead><tr><th>Chromosome</th><th>Length (bp)</th><th>{escape(data.labels[0])} coverage</th><th>{escape(data.labels[1])} coverage</th><th>Assignment</th></tr></thead><tbody>{table}</tbody></table></section><p>Window: {data.window_size:,} bp. Coordinates: zero-based, half-open. Local windows and whole-chromosome assignments use their respective union coverage.</p></main>''' + '''<script>const search=document.getElementById('search'),source=document.getElementById('source');function filter(){document.querySelectorAll('tbody tr').forEach(row=>{row.hidden=!(row.cells[0].textContent.toLowerCase().includes(search.value.toLowerCase())&&(!source.value||row.dataset.source===source.value));});}search.addEventListener('input',filter);source.addEventListener('change',filter);</script></html>'''
    path.write_text(document, encoding="utf-8")


def make_report(out_dir, labels, records, intervals, assignments, all_records,
                min_coverage=0.25, dominance=0.6, window_size=0, alignments=None):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from .plotting import STYLE, partition_pages, matrix_pages
    from .comparison import comparison_pages, write_links
    data = prepare_data(labels, records, intervals, assignments, all_records, min_coverage, dominance, window_size)
    alignments = list(alignments or [])
    with open(out_dir / "window_provenance.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["seq_id", "start", "end", *(f"{label}_covered_bp" for label in data.labels), "window_assignment"])
        writer.writerows(data.windows)
    matrix = alignment_matrix(data, alignments, intervals)
    with open(out_dir / "alignment_matrix.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["source", "original_contig", *records])
        for key, values in matrix.items():
            writer.writerow([*key, *(f"{values.get(n, 0):.8f}" for n in records)])
    write_links(out_dir / "alignment_links.tsv", alignments, assignments)
    sections = []
    with plt.rc_context(STYLE):
        with PdfPages(out_dir / "partition.pdf") as pdf:
            sections.append(("1. Source evidence and chromosome provenance", partition_pages(data, pdf)))
        with PdfPages(out_dir / "alignment_matrix.pdf") as pdf:
            sections.append(("2. Original-contig alignment matrix", matrix_pages(data, matrix, pdf)))
        with PdfPages(out_dir / "hap_comparison.pdf") as pdf:
            sections.append(("3. Original hifiasm to recovered hap correspondence", comparison_pages(data, alignments, pdf)))
    write_html(out_dir / "report.html", data, sections)
