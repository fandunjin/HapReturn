"""Coordinate-wise provenance-first fusion; alignment support is never added twice."""
import csv
import html
import json
from collections import Counter, defaultdict


def fuse_evidence(result, intervals, names, labels):
    """Return disjoint support and auditable 0-based half-open intervals.

    Input comparison placements are 1-based inclusive; convert to PAF coordinates.
    AGP is construction provenance, not verified parental phase. Only unique
    pre/post traceable placements with known orientation supply direct support.
    Ambiguous placements are quarantined. Alignment fills unknown sequence
    components only when exactly one source covers the target interval.
    Gaps never supply evidence. Opposite-only alignment on known AGP is flagged
    but does not erase direct provenance.
    """
    components = defaultdict(list)
    for atom in result['atoms']:
        for loc in atom['after']:
            components[loc['object']].append((atom, loc))
    support, audits, segments = {}, {}, []
    for name in names:
        events = defaultdict(list)
        parts = components[name]
        for i, (_, loc) in enumerate(parts):
            events[loc['start'] - 1].append(('part', i, 1))
            events[loc['end']].append(('part', i, -1))
        for label in labels:
            for start, end in intervals[label].get(name, []):
                events[start].append(('aln', label, 1))
                events[end].append(('aln', label, -1))
        active, align = Counter(), Counter()
        bp = {label: 0 for label in labels}
        stats = Counter()
        points = sorted(events)
        for index, start in enumerate(points[:-1]):
            for kind, key, delta in events[start]:
                (active if kind == 'part' else align)[key] += delta
            end = points[index + 1]
            present = [i for i, count in active.items() if count]
            if not present:
                continue  # AGP gap, even if a PAF spans it.
            if len(present) != 1:
                raise ValueError('Overlapping AGP target components')
            atom, loc = parts[present[0]]
            hits = [label for label in labels if align[label] > 0]
            direct = atom['source_hap'] if atom['source_hap'] in labels else ''
            unsafe = (atom['mapping_status'] == 'ambiguous' or loc['strand'] not in ('+', '-')
                      or any(x['strand'] not in ('+', '-') for x in atom['before']))
            traceable = atom['mapping_status'] == 'unique' and not unsafe
            conflict = False
            if unsafe:
                chosen, method = 'unresolved', 'unsafe_agp'
            elif direct and traceable:
                chosen, method = direct, 'agp_direct'
                conflict = len(hits) == 1 and hits[0] != direct
            elif len(hits) == 1:
                chosen, method = hits[0], 'alignment_fill'
            else:
                chosen = 'unresolved'
                method = 'alignment_ambiguous' if len(hits) > 1 else 'unknown'
            size = end - start
            stats['sequence_bp'] += size
            stats[method + '_bp'] += size
            stats['conflict_bp'] += size if conflict else 0
            if chosen in bp:
                bp[chosen] += size
            segments.append({'seq_id': name, 'start': start, 'end': end, 'length_bp': size,
                             'source_contig': atom['source_contig'], 'agp_source': direct or 'unresolved',
                             'alignment_sources': ','.join(hits), 'source': chosen, 'method': method,
                             'conflict': conflict})
        expected = sum(loc['end'] - loc['start'] + 1 for _, loc in parts)
        if stats['sequence_bp'] != expected:
            raise ValueError(f"Hybrid coordinate accounting mismatch: {name}")
        support[name], audits[name] = bp, dict(stats)
    return support, audits, segments


def write_integration(out_dir, labels, names, support, audits, segments, minimap,
                      agp, evidence, assignments, config):
    from .agp import write_tsv
    rows = []
    for name in names:
        stats = audits[name]
        flags = [key for key in ('conflict_bp', 'unsafe_agp_bp', 'alignment_ambiguous_bp') if stats.get(key)]
        if evidence[name] == 'unresolved':
            flags.append('insufficient_or_mixed_evidence')
        if assignments[name] != evidence[name]:
            flags.append('constraint_changed_evidence')
        rows.append({'seq_id': name, 'minimap_assignment': minimap[name],
                     'agp_assignment': agp[name], 'hybrid_independent_assignment': evidence[name],
                     'assignment': assignments[name],
                     **{label + '_support_bp': support[name][label] for label in labels},
                     **{key: stats.get(key, 0) for key in ('sequence_bp', 'agp_direct_bp', 'alignment_fill_bp',
                        'unknown_bp', 'unsafe_agp_bp', 'alignment_ambiguous_bp', 'conflict_bp')},
                     'changed_from_minimap': assignments[name] != minimap[name],
                     'review_flags': ';'.join(flags)})
    write_tsv(out_dir / 'evidence_comparison.tsv', rows, list(rows[0]))
    write_tsv(out_dir / 'hybrid_intervals.tsv',
              [dict(r, start=r['start'] + 1) for r in segments], list(segments[0]) if segments else
              ['seq_id', 'start', 'end', 'length_bp', 'source_contig', 'agp_source', 'alignment_sources', 'source', 'method', 'conflict'])
    (out_dir / 'evidence_config.json').write_text(json.dumps(config, indent=2), encoding='utf-8')
    render_integration_report(out_dir, generate_figures=False)
    return rows


def render_integration_report(out_dir, sections=None, generate_figures=True):
    """Render from saved evidence, so report refresh never changes decisions."""
    from pathlib import Path
    import re
    out_dir = Path(out_dir)
    def read_table(name):
        with (out_dir / name).open(encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    config = json.loads((out_dir / 'evidence_config.json').read_text(encoding='utf-8'))
    if config.get('evidence_mode') != 'hybrid':
        return
    rows = read_table('evidence_comparison.tsv')
    intervals = read_table('hybrid_intervals.tsv')
    assignments = read_table('assignment.tsv')
    lengths = {r['seq_id']: int(r['length']) for r in assignments}
    labels = config.get('labels') or [k[:-11] for k in rows[0] if k.endswith('_support_bp')]
    for row in rows:
        row['length'] = lengths[row['seq_id']]
        for key in row:
            if key.endswith('_bp'):
                row[key] = int(row[key])
    for row in intervals:
        for key in ('start', 'end', 'length_bp'):
            row[key] = int(row[key])
        row['conflict'] = str(row['conflict']).lower() == 'true'
    from .joint_figures import make_joint_figures
    if generate_figures:
        make_joint_figures(out_dir, labels, rows, intervals)
    figures = {}
    for key, stem in [('1', 'joint_evidence_tracks'), ('2', 'joint_evidence_matrix')]:
        path = out_dir / (stem + '.svg')
        if generate_figures and path.exists():
            content = path.read_text(encoding='utf-8')
            figures[key] = content[content.index('<svg'):]
    payload = json.dumps({'labels': labels, 'rows': rows, 'segments': intervals}, ensure_ascii=False).replace('<', '\u003c')
    # Escape script terminators in user-provided sequence IDs.
    payload = payload.replace(chr(60), chr(92) + 'u003c').replace(chr(0x2028), chr(92)+'u2028').replace(chr(0x2029), chr(92)+'u2029')
    template = Path(__file__).with_name('hybrid_report_template.html').read_text(encoding='utf-8')
    for key in ('1', '2'):
        template = template.replace('<!--FIGURE_'+key+'-->', figures.get(key, '<p class="empty">Joint evidence is not available for this view. Regenerate the report from saved alignment files.</p>'))
    template = template.replace('/*PAYLOAD*/null', payload)
    (out_dir / 'evidence_report.html').write_text(template, encoding='utf-8')
