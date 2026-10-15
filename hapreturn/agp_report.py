"""Portable interactive SVG report and vector PDF for AGP comparisons."""
import html
import json
from collections import defaultdict
from pathlib import Path

from .agp import OUTPUT_NAMES


COLORS = {"hap1": "#3274A1", "hap2": "#D88A3D", "unresolved": "#A5ABB3"}
CHANGE_COLOR = "#8a3fab"
LABELS = {
    "equivalent_scaffold_rename": "等价 scaffold 重命名",
    "whole_scaffold_reverse": "整条 scaffold 反向",
    "scaffold_order_changed": "scaffold 排列顺序改变",
    "cross_scaffold": "跨 scaffold 区间转移",
    "orientation_changed": "区间相对方向改变",
    "relative_order_changed": "内部相对顺序改变",
    "scaffold_split": "scaffold 拆分",
    "scaffold_merge": "scaffold 合并",
    "complex_scaffold_reorganization": "复杂 scaffold 重组",
    "junction_added": "新增连接", "junction_removed": "消失连接", "gap_changed": "Gap 改变",
    "component_boundary_added": "新增 component 边界", "component_boundary_removed": "移除 component 边界",
    "before_only": "仅调整前找到", "after_only": "仅调整后找到",
    "ambiguous": "重复区间／对应歧义", "orientation_unknown": "方向未知",
    "unanchored_scaffold": "缺少唯一 scaffold 锚点",
}


def palette(result):
    labels = result.get("agp_provenance", {}).get("labels", result.get("source_labels", []))
    colors = dict(zip(list(labels)+["unresolved"], [COLORS["hap1"], COLORS["hap2"], COLORS["unresolved"]])) if labels else dict(COLORS)
    extras = ["#249a83", "#b45a7c", "#7b75bd", "#9a7732"]
    for label in sorted({x["source_hap"] for x in result["atoms"]} - colors.keys()):
        colors[label] = extras[len(colors) % len(extras)]
    return colors


def default_selection(result, minimum):
    before = {n for n, size in result["before_lengths"].items() if size >= minimum}
    after = {n for n, size in result["after_lengths"].items() if size >= minimum}
    # Always include the counterpart of a visible interval, even if very short.
    # Iterate to preserve the complete connected neighborhood for split/merge cases.
    changed = True
    while changed:
        changed = False
        for atom in result["atoms"]:
            bs = {p["object"] for p in atom["before"]}
            aas = {p["object"] for p in atom["after"]}
            if bs & before or aas & after:
                if not bs <= before or not aas <= after:
                    before.update(bs)
                    after.update(aas)
                    changed = True
    if not before and not after:
        before = {max(result["before_lengths"], key=result["before_lengths"].get)}
        after = {p["object"] for a in result["atoms"] for p in a["after"]
                 if any(b["object"] in before for b in a["before"])}
    return [n for n in result["before_lengths"] if n in before], [n for n in result["after_lengths"] if n in after]


def geometry(result, before, after):
    lengths = result["before_lengths"], result["after_lengths"]
    names = before, after
    count = max(len(before), len(after), 1)
    desired = max(400, count * 28)
    lo, hi = 0.0, desired / max(list(lengths[0].values()) + list(lengths[1].values()))
    for _ in range(60):
        scale = (lo + hi) / 2
        extent = max(sum(max(12, lengths[i][name] * scale) + 14 for name in names[i]) for i in (0, 1))
        if extent > desired:
            hi = scale
        else:
            lo = scale
    positions = [{}, {}]
    for i in (0, 1):
        cursor = 98.0
        for name in names[i]:
            body = lengths[i][name] * lo
            reserve = max(12, body)
            positions[i][name] = (cursor + (reserve - body) / 2, body)
            cursor += reserve + 14
    return {"positions": positions, "scale": lo, "width": 1240, "height": desired + 160,
            "left_x": 290, "right_x": 950}


def visible_links(result, geom):
    for atom in result["atoms"]:
        if atom["mapping_status"] != "unique":
            continue
        bp, ap = atom["before"][0], atom["after"][0]
        if bp["object"] not in geom["positions"][0] or ap["object"] not in geom["positions"][1]:
            continue
        ys = []
        for version, p in enumerate((bp, ap)):
            bottom, _ = geom["positions"][version][p["object"]]
            a = bottom + (p["start"] - 1) * geom["scale"]
            b = bottom + p["end"] * geom["scale"]
            ys.append((b, a) if p["strand"] == "-" else (a, b))
        yield atom, ys[0], ys[1]


def static_svg(result, before, after, title="AGP curation correspondence", left_title="Before", right_title="After"):
    esc, colors = html.escape, palette(result)
    geom = geometry(result, before, after)
    width, height, scale = geom["width"], geom["height"], geom["scale"]
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="AGP before and after source interval correspondence">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#253247}.link:hover{fill-opacity:.8!important;stroke:#253247;stroke-width:1}</style>',
             f'<text x="24" y="28" font-size="18" font-weight="bold">{esc(title)}</text>',
             f'<text x="24" y="49" font-size="12">Before: {len(before)}/{len(result["before_lengths"])} tracks · After: {len(after)}/{len(result["after_lengths"])} tracks · shared bp scale</text>',
             f'<text x="290" y="79" font-size="13" text-anchor="middle">{esc(left_title)}</text><text x="950" y="79" font-size="13" text-anchor="middle">{esc(right_title)}</text>']
    for atom, (b0, b1), (a0, a1) in visible_links(result, geom):
        flags = ", ".join(atom["flags"]) or "No internal structural difference"
        tip = esc(f'{atom["source_contig"]}:{atom["source_start"]}-{atom["source_end"]} | {atom["length_bp"]:,} bp | {flags}')
        path = f'M 298 {b0} C 520 {b0},720 {a0},942 {a0} L 942 {a1} C 720 {a1},520 {b1},298 {b1} Z'
        parts.append(f'<path class="link" data-segment="{atom["segment_id"]}" d="{path}" fill="{colors[atom["source_hap"]]}" fill-opacity=".24"><title>{tip}</title></path>')
    for index, version in enumerate(("before", "after")):
        x = geom["left_x"] if index == 0 else geom["right_x"]
        for name, (top, body) in geom["positions"][index].items():
            parts.append(f'<rect x="{x-8}" y="{top}" width="16" height="{body}" fill="#eef1f5" stroke="#566274" stroke-width=".6"/>')
            label_x, anchor = (x - 18, "end") if index == 0 else (x + 18, "start")
            label = esc(f'{name} | {result[version + "_lengths"][name] / 1e6:.3f} Mb')
            parts.append(f'<text x="{label_x}" y="{top + body / 2 + 3.5}" font-size="11" text-anchor="{anchor}">{label}</text>')
        for atom in result["atoms"]:
            for p in atom[version]:
                if p["object"] not in geom["positions"][index]:
                    continue
                top = geom["positions"][index][p["object"]][0] + (p["start"] - 1) * scale
                body = atom["length_bp"] * scale
                tip = esc(f'{atom["source_contig"]}:{atom["source_start"]}-{atom["source_end"]} → {p["object"]}:{p["start"]}-{p["end"]} ({p["strand"]}) | {atom["mapping_status"]}')
                stroke = CHANGE_COLOR if atom["flags"] else "#566274"
                parts.append(f'<rect data-segment="{atom["segment_id"]}" x="{x-7}" y="{top}" width="14" height="{body}" fill="{colors[atom["source_hap"]]}" stroke="{stroke}" stroke-width=".5"><title>{tip}</title></rect>')
                if body > 14 and p["strand"] in "+-":
                    y = top + body / 2
                    direction = 1 if p["strand"] == "+" else -1
                    parts.append(f'<path d="M {x-3} {y-3*direction} L {x} {y+3*direction} L {x+3} {y-3*direction}" fill="none" stroke="white" stroke-width="1.1"/>')
    legend_y = height - 34
    for i, (label, color) in enumerate(colors.items()):
        parts.append(f'<rect x="{24+i*145}" y="{legend_y-10}" width="12" height="12" fill="{color}"/><text x="{42+i*145}" y="{legend_y}" font-size="12">{esc(label)}</text>')
    parts.append(f'<text x="24" y="{height-12}" font-size="11">Colors: source labels (name convention unless a source map was supplied). Only unique interval links are drawn.</text></svg>')
    return "".join(parts)


def plot_pdf(result, out_dir, minimum, filename="curation_comparison.pdf", title="AGP curation correspondence",
             left_title="Before", right_title="After"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import Rectangle, PathPatch
    from .plotting import STYLE
    colors = palette(result)
    before, after = default_selection(result, minimum)
    # Append every remaining track in batches, closing over source correspondence.
    # No source record is silently dropped from the PDF.
    selections = [(before, after, "Overview")]
    done_b, done_a = set(before), set(after)
    while len(done_b) < len(result["before_lengths"]) or len(done_a) < len(result["after_lengths"]):
        bs = set([n for n in result["before_lengths"] if n not in done_b][:18])
        aas = set([n for n in result["after_lengths"] if n not in done_a][:18])
        changed = True
        while changed:
            changed = False
            for atom in result["atoms"]:
                bnames = {p["object"] for p in atom["before"]}
                anames = {p["object"] for p in atom["after"]}
                if bnames & bs or anames & aas:
                    if not bnames <= bs or not anames <= aas:
                        bs.update(bnames)
                        aas.update(anames)
                        changed = True
        bnames = [n for n in result["before_lengths"] if n in bs]
        anames = [n for n in result["after_lengths"] if n in aas]
        selections.append((bnames, anames, "Additional tracks"))
        done_b.update(bs)
        done_a.update(aas)
    with plt.rc_context(STYLE), PdfPages(Path(out_dir) / filename) as pdf:
        for page, (before, after, page_title) in enumerate(selections, 1):
            geom = geometry(result, before, after)
            w, h = geom["width"], geom["height"]
            fig = plt.figure(figsize=(12.4, h / 100))
            ax = fig.add_axes([0, 0, 1, 1])
            ax.set(xlim=(0, w), ylim=(h, 0))
            ax.axis("off")
            ax.text(24, 28, f"{title} | {page_title} | page {page}/{len(selections)}", fontsize=12, weight="bold")
            ax.text(24, 49, f"Before {len(before)}/{len(result['before_lengths'])} tracks; after {len(after)}/{len(result['after_lengths'])} tracks. Shared bp scale within this page.", fontsize=8)
            for x, label in ((290, left_title), (950, right_title)):
                ax.text(x, 79, label, fontsize=9, ha="center")
            for atom, (b0, b1), (a0, a1) in visible_links(result, geom):
                verts = [(298, b0), (520, b0), (720, a0), (942, a0), (942, a1), (720, a1), (520, b1), (298, b1), (298, b0)]
                codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.CLOSEPOLY]
                ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=colors[atom["source_hap"]], edgecolor="none", alpha=.24))
            for index, version in enumerate(("before", "after")):
                x = 290 if index == 0 else 950
                for name, (top, body) in geom["positions"][index].items():
                    ax.add_patch(Rectangle((x - 8, top), 16, body, facecolor="#eef1f5", edgecolor="#566274", lw=.4))
                    ax.text(x - 18 if index == 0 else x + 18, top + body / 2,
                            f"{name} | {result[version + '_lengths'][name] / 1e6:.3f} Mb", fontsize=7,
                            ha="right" if index == 0 else "left", va="center")
                for atom in result["atoms"]:
                    for p in atom[version]:
                        if p["object"] not in geom["positions"][index]:
                            continue
                        top = geom["positions"][index][p["object"]][0] + (p["start"] - 1) * geom["scale"]
                        body = atom["length_bp"] * geom["scale"]
                        ax.add_patch(Rectangle((x - 7, top), 14, body, facecolor=colors[atom["source_hap"]],
                                               edgecolor=CHANGE_COLOR if atom["flags"] else "#566274", lw=.35))
                        if body > 14 and p["strand"] in "+-":
                            y, direction = top + body / 2, 1 if p["strand"] == "+" else -1
                            ax.plot([x-3, x, x+3], [y-3*direction, y+3*direction, y-3*direction], color="white", lw=.8)
            for i, (label, color) in enumerate(colors.items()):
                ax.add_patch(Rectangle((24 + i*145, h-44), 12, 12, facecolor=color, edgecolor="none"))
                ax.text(42 + i*145, h-34, label, fontsize=8)
            ax.text(24, h-12, "Source labels use contig names unless a source map was supplied. All tracks appear across these pages. Not an edit history or Hi-C quality assessment.", fontsize=7)
            pdf.savefig(fig)
            plt.close(fig)


def write_html(result, out_dir, minimum):
    out_dir = Path(out_dir)
    before, after = default_selection(result, minimum)
    svg = static_svg(result, before, after)
    (out_dir / "curation_comparison.svg").write_text(svg, encoding="utf-8")
    from .agp_provenance import source_view, PROVENANCE_OUTPUTS
    view = source_view(result)
    origin_svgs = {}
    if view["atoms"]:
        args = dict(title="Original hifiasm components to current chromosomes", left_title="Original contigs", right_title="Current chromosomes")
        source_svg = static_svg(view, list(view["before_lengths"]), list(view["after_lengths"]), **args)
        (out_dir / "agp_source_correspondence.svg").write_text(source_svg, encoding="utf-8")
        origin_svgs["all"] = source_svg
        for name in result["agp_provenance"]["selected"]:
            source_names = {a["source_contig"] for a in view["atoms"] if any(p["object"] == name for p in a["after"])}
            origin_svgs[name] = static_svg(view, [n for n in view["before_lengths"] if n in source_names], [name], **args)
    payload = dict(result, colors=palette(result), labels=LABELS, default_before=before, default_after=after, origin_svgs=origin_svgs)
    # Escape script delimiters and JS line separators in untrusted AGP identifiers.
    data = json.dumps(payload, ensure_ascii=False).replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    files = OUTPUT_NAMES + PROVENANCE_OUTPUTS
    files += result.get("sequence_export", {}).get("files", [])
    downloads = " ".join(f'<a href="{html.escape(name, quote=True)}" download>{html.escape(name)}</a>' for name in files if name != "curation_report.html")
    template = Path(__file__).with_name("agp_report_template.html").read_text(encoding="utf-8")
    document = template.replace("<!--INITIAL_SVG-->", svg).replace("<!--DOWNLOADS-->", downloads).replace("/*DATA_PAYLOAD*/null", data)
    (out_dir / "curation_report.html").write_text(document, encoding="utf-8")


def provenance_pdf(result, out_dir):
    """Length-scaled current chromosomes, grouped by AGP assignment, preserving mixtures."""
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import Rectangle
    from .plotting import STYLE
    p, colors = result["agp_provenance"], palette(result)
    rows = [r for r in p["chromosomes"] if r["selected"]]
    if not rows:
        return
    max_length = max(r["length"] for r in rows)
    parts = defaultdict(list)
    for part in p["components"]:
        parts[part["seq_id"]].append(part)
    categories = p["labels"] + ["unresolved"]
    with plt.rc_context(STYLE), PdfPages(Path(out_dir)/"agp_provenance.pdf") as pdf:
        for source in categories:
            members = [r for r in rows if r["assignment"] == source]
            if not members:
                continue
            height = max(3.2, 1.65+len(members)*.42)
            fig = plt.figure(figsize=(11.5, height))
            ax = fig.add_axes([.02,.05,.96,.90])
            ax.set(xlim=(0,1100),ylim=(len(members)+2,-1))
            ax.axis("off")
            ax.text(12,-.65,f"AGP chromosome provenance | recovered {source} | n={len(members)}",fontsize=12,weight="bold")
            ax.text(12,-.15,"Color = actual original component source; whole-chromosome assignment does not erase mixed segments.",fontsize=8)
            for i,row in enumerate(members):
                x,width = 170,row["length"]/max_length*670
                ax.text(155,i+.65,row["seq_id"],ha="right",va="center",fontsize=8)
                ax.add_patch(Rectangle((x,i+.45),width,.38,facecolor="#eef1f5",edgecolor="#566274",lw=.35))
                for part in parts[row["seq_id"]]:
                    ax.add_patch(Rectangle((x+(part["start"]-1)/max_length*670,i+.45),part["length_bp"]/max_length*670,.38,
                                           facecolor=colors.get(part["source_hap"],COLORS["unresolved"]),edgecolor="none"))
                fractions=" | ".join(f"{label}: {row[label+'_fraction']*100:.2f}%" for label in p["labels"])
                ax.text(860,i+.65,fractions,fontsize=7,va="center")
            y=len(members)+.5
            for i,(label,color) in enumerate(colors.items()):
                ax.add_patch(Rectangle((12+i*160,y),14,.30,facecolor=color,edgecolor="none"))
                ax.text(32+i*160,y+.15,label,fontsize=8,va="center")
            ax.text(12,len(members)+1.25,f"Assignment mode: {p['assignment_mode']}; source dominance >= {p['dominance']:.2f}; gap bases excluded from fractions. All bars share a bp scale across pages.",fontsize=7)
            pdf.savefig(fig)
            plt.close(fig)
