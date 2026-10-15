"""Coordinate-exact AGP state comparison. Internal coordinates are 0-based half-open.

Only unique source intervals participate in structural inference. The module
compares states, not edit history, sequence identity, phasing or Hi-C quality.
"""
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class Component:
    object: str
    start: int
    end: int
    part: int
    contig: str
    source_start: int
    source_end: int
    strand: str

    def project(self, start, end):
        if self.strand == "-":
            return self.start + self.source_end - end, self.start + self.source_end - start
        return self.start + start - self.source_start, self.start + end - self.source_start


@dataclass
class Assembly:
    path: Path
    lengths: dict = field(default_factory=dict)
    components: list = field(default_factory=list)
    gaps: dict = field(default_factory=lambda: defaultdict(list))


def read_agp(path):
    """Validate tiling, lengths, part order and standard AGP component/gap fields."""
    result = Assembly(Path(path))
    previous, closed, part = None, set(), 0
    with result.path.open(encoding="utf-8-sig") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            try:
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 9 or any(not x for x in fields):
                    raise ValueError("expected 9 nonempty tab-separated columns")
                obj, a, b, p, kind, name, c, d, strand = fields
                start, end, current_part = int(a) - 1, int(b), int(p)
                if obj != previous:
                    if obj in closed:
                        raise ValueError("object appears in noncontiguous blocks")
                    if previous is not None:
                        closed.add(previous)
                    previous, part = obj, 0
                if start != result.lengths.get(obj, 0) or end <= start or current_part != part + 1:
                    raise ValueError("object coordinates must tile from 1; part numbers must be consecutive")
                if kind in ("N", "U"):
                    if int(name) != end - start or (kind == "U" and int(name) != 100):
                        raise ValueError("gap length disagrees with span or U gap is not 100 bp")
                    if c not in {"scaffold", "contig", "centromere", "short_arm", "heterochromatin", "telomere", "repeat", "contamination"} or d not in {"yes", "no"}:
                        raise ValueError("invalid gap type or linkage")
                    result.gaps[obj].append((start, end, kind, name, c, d, strand))
                elif kind in {"A", "D", "F", "G", "O", "P", "W"}:
                    source_start, source_end = int(c) - 1, int(d)
                    if source_start < 0 or source_end <= source_start or source_end - source_start != end - start:
                        raise ValueError("component and object spans disagree")
                    if strand not in {"+", "-", "?", "0", "na"}:
                        raise ValueError("invalid orientation")
                    result.components.append(Component(obj, start, end, current_part, name, source_start, source_end, strand))
                else:
                    raise ValueError("invalid component type")
                result.lengths[obj], part = end, current_part
            except ValueError as exc:
                raise ValueError(f"{path}:{number}: {exc}") from exc
    if not result.components:
        raise ValueError(f"{path}: no sequence components")
    return result


def read_source_map(path):
    if path is None:
        return None
    result = {}
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"contig_id", "source"} <= set(reader.fieldnames or []):
            raise ValueError("Source map requires contig_id and source columns")
        for row in reader:
            if not row["contig_id"] or not row["source"] or row["contig_id"] in result:
                raise ValueError("Empty or duplicate contig in source map")
            result[row["contig_id"]] = row["source"]
    return result


def source_label(contig, mapping):
    if mapping is not None:
        return mapping.get(contig, "unresolved"), "source_map" if contig in mapping else "unresolved"
    match = re.fullmatch(r"(?:lcl\|)?h([12])tg\d+[A-Za-z]*", contig)
    return (f"hap{match[1]}", "contig_name_convention") if match else ("unresolved", "unresolved")


def atomic_intervals(before, after, mapping=None):
    """Sweep source endpoints; never expand duplicate placements into Cartesian links."""
    events = defaultdict(lambda: defaultdict(list))
    for version, assembly in (("before", before), ("after", after)):
        for index, component in enumerate(assembly.components):
            events[component.contig][component.source_start].append((version, index, True))
            events[component.contig][component.source_end].append((version, index, False))
    atoms = []
    for contig in sorted(events):
        active = {"before": set(), "after": set()}
        points = sorted(events[contig])
        for start, end in zip(points, points[1:]):
            for version, index, add in events[contig][start]:
                if add:
                    active[version].add(index)
                else:
                    active[version].remove(index)
            if not any(active.values()):
                continue
            identity = json.dumps([contig, start, end], ensure_ascii=False).encode("utf-8")
            label, method = source_label(contig, mapping)
            atom = {"segment_id": "S" + hashlib.sha256(identity).hexdigest()[:16],
                    "source_contig": contig, "source_start": start + 1, "source_end": end,
                    "length_bp": end - start, "source_hap": label, "source_method": method,
                    "flags": [], "events": []}
            for version, assembly in (("before", before), ("after", after)):
                placements = []
                for index in sorted(active[version]):
                    component = assembly.components[index]
                    a, b = component.project(start, end)
                    placements.append({"object": component.object, "start": a + 1, "end": b,
                                       "strand": component.strand, "part": component.part,
                                       "coordinate_status": "exact" if component.strand in "+-" else "orientation_assumed_plus"})
                atom[version] = placements
            atom["mapping_status"] = ("ambiguous" if max(len(atom["before"]), len(atom["after"])) > 1 else
                                      "unique" if atom["before"] and atom["after"] else
                                      "before_only" if atom["before"] else "after_only")
            atoms.append(atom)
    return atoms


def object_atoms(atoms, version):
    result = defaultdict(list)
    for atom in atoms:
        for placement in atom[version]:
            result[placement["object"]].append((atom, placement))
    for items in result.values():
        items.sort(key=lambda item: item[1]["start"])
    return result


def reciprocal_anchors(atoms):
    overlaps = Counter()
    for atom in atoms:
        if atom["mapping_status"] == "unique":
            overlaps[atom["before"][0]["object"], atom["after"][0]["object"]] += atom["length_bp"]
    def best(axis):
        grouped = defaultdict(list)
        for (b, a), bp in overlaps.items():
            grouped[(b, a)[axis]].append(((a, b)[axis], bp))
        result = {}
        for name, items in grouped.items():
            ranked = sorted(items, key=lambda item: -item[1])
            if len(ranked) == 1 or ranked[0][1] > ranked[1][1]:
                result[name] = ranked[0][0]
        return result
    bbest, abest = best(0), best(1)
    anchors = {b: a for b, a in bbest.items() if abest.get(a) == b}
    return overlaps, anchors


def sequence_equivalence(before_items, after_items):
    """Return forward/reverse only when complete, uniquely mapped tracks agree."""
    if len(before_items) != len(after_items) or any(a["mapping_status"] != "unique" for a, _ in before_items + after_items):
        return "changed"
    first = [(a["segment_id"], p["strand"]) for a, p in before_items]
    second = [(a["segment_id"], p["strand"]) for a, p in after_items]
    if any(strand not in "+-" for _, strand in first + second):
        return "orientation_unknown"
    if first == second:
        return "forward"
    if [(identity, "+" if strand == "-" else "-") for identity, strand in reversed(first)] == second:
        return "reverse"
    return "changed"


def order_changed_ids(first, second):
    """Mark inversion participants, not downstream coordinate shifts or insertions."""
    ranks = {identity: index for index, identity in enumerate(first)}
    common = [identity for identity in second if identity in ranks]
    values = [ranks[x] for x in common]
    prefix, largest = [], -1
    for value in values:
        prefix.append(largest)
        largest = max(largest, value)
    suffix, smallest = [], len(first)
    for value in reversed(values):
        suffix.append(smallest)
        smallest = min(smallest, value)
    return {identity for identity, value, left, right in zip(common, values, prefix, reversed(suffix))
            if value < left or value > right}


def junctions(assembly, by_object):
    """Canonical source-coordinate ends make junctions invariant to whole reversal."""
    result = defaultdict(list)
    skipped = 0
    for obj, items in by_object.items():
        for (left, lp), (right, rp) in zip(items, items[1:]):
            if (left["mapping_status"] == "ambiguous" or right["mapping_status"] == "ambiguous"
                    or lp["strand"] not in "+-" or rp["strand"] not in "+-"):
                skipped += 1
                continue
            le = (left["segment_id"], "end" if lp["strand"] == "+" else "start")
            re_ = (right["segment_id"], "start" if rp["strand"] == "+" else "end")
            key = tuple(sorted((le, re_)))
            gaps = [g for g in assembly.gaps.get(obj, []) if g[0] >= lp["end"] and g[1] <= rp["start"] - 1]
            # Order-independent gap metadata: track reversal must not invent changes.
            signature = sorted(tuple(g[2:]) for g in gaps)
            result[key].append({"object": obj, "left_end": lp["end"], "right_start": rp["start"],
                                "gap_bp": rp["start"] - 1 - lp["end"], "gap_signature": signature,
                                "segments": [left["segment_id"], right["segment_id"]]})
    return result, skipped


def compare_agp(before, after, source_map=None, assignments=None):
    atoms = atomic_intervals(before, after, source_map)
    by_before, by_after = object_atoms(atoms, "before"), object_atoms(atoms, "after")
    overlaps, anchors = reciprocal_anchors(atoms)
    equivalences = {b: sequence_equivalence(by_before[b], by_after[a]) for b, a in anchors.items()}
    events, qc = [], []
    atom_by_id = {a["segment_id"]: a for a in atoms}
    def event(kind, ids, before_objects=(), after_objects=(), note="Observed state difference"):
        ids = sorted(set(ids))
        row = {"event_id": f"E{len(events) + 1:04d}", "type": kind,
               "segment_ids": ids, "before_objects": list(before_objects), "after_objects": list(after_objects),
               "length_bp": sum(atom_by_id[i]["length_bp"] for i in ids), "interpretation": note}
        events.append(row)
        for identity in ids:
            atom_by_id[identity]["events"].append(row["event_id"])
        return row

    for b, a in anchors.items():
        ids = [atom["segment_id"] for atom, _ in by_before[b] if atom["mapping_status"] == "unique" and atom["after"][0]["object"] == a]
        eq = equivalences[b]
        if eq == "reverse":
            event("whole_scaffold_reverse", ids, [b], [a], "Complete source order and known orientations reverse; internal topology equivalent")
        if eq in {"forward", "reverse"} and b != a:
            event("equivalent_scaffold_rename", ids, [b], [a], "Sequence layout equivalent; gap changes are reported separately")

    before_rank = {name: index for index, name in enumerate(before.lengths, 1)}
    after_rank = {name: index for index, name in enumerate(after.lengths, 1)}
    ordered_before = sorted(anchors, key=before_rank.get)
    ordered_after = sorted(anchors, key=lambda b: after_rank[anchors[b]])
    reordered_scaffolds = order_changed_ids(ordered_before, ordered_after)
    for b in ordered_before:
        if b in reordered_scaffolds:
            a = anchors[b]
            ids = [x["segment_id"] for x, _ in by_before[b] if x["mapping_status"] == "unique" and x["after"][0]["object"] == a]
            event("scaffold_order_changed", ids, [b], [a],
                  f"Relative FASTA/AGP record order changed (rank {before_rank[b]} to {after_rank[a]}); not an internal rearrangement")

    reordered = set()
    for b, a in anchors.items():
        first = [atom["segment_id"] for atom, _ in by_before[b] if atom["mapping_status"] == "unique" and atom["after"][0]["object"] == a]
        second = [atom["segment_id"] for atom, _ in by_after[a] if atom["mapping_status"] == "unique" and atom["before"][0]["object"] == b]
        if equivalences[b] == "reverse":
            second.reverse()
        reordered.update(order_changed_ids(first, second))

    for atom in atoms:
        status = atom["mapping_status"]
        if status != "unique":
            atom["flags"].append(status)
            if status == "ambiguous":
                qc.append({"severity": "review", "code": "duplicate_source_interval", "segment_id": atom["segment_id"],
                           "detail": "Multiple placements; structural inference excluded, all occurrences retained"})
            else:
                qc.append({"severity": "review", "code": "unmatched_source_interval", "segment_id": atom["segment_id"],
                           "detail": f"{status}; absent from the other AGP, not automatically a sequence insertion/deletion"})
        else:
            bp, ap = atom["before"][0], atom["after"][0]
            b, a = bp["object"], ap["object"]
            if b not in anchors:
                atom["flags"].append("unanchored_scaffold")
            elif anchors[b] != a:
                atom["flags"].append("cross_scaffold")
            if bp["strand"] not in "+-" or ap["strand"] not in "+-":
                atom["flags"].append("orientation_unknown")
                qc.append({"severity": "review", "code": "unknown_orientation", "segment_id": atom["segment_id"],
                           "detail": "Within-component coordinates assume plus for display; direction/junction inference unavailable"})
            elif bp["strand"] != ap["strand"] and not (anchors.get(b) == a and equivalences[b] == "reverse"):
                atom["flags"].append("orientation_changed")
            if atom["segment_id"] in reordered:
                atom["flags"].append("relative_order_changed")
        atom["after_assignment"] = (assignments or {}).get(atom["after"][0]["object"], "") if len(atom["after"]) == 1 else ""

    # Group consecutive atomic intervals sharing exactly the same difference flags
    # and object pair, without claiming a unique sequence of manual operations.
    previous, group = None, []
    def flush():
        if group:
            first = group[0]
            event("+".join(first["flags"]), [x["segment_id"] for x in group],
                  sorted({p["object"] for x in group for p in x["before"]}),
                  sorted({p["object"] for x in group for p in x["after"]}),
                  "Consecutive source intervals with the same observed differences; not an edit-history reconstruction")
    for atom in atoms:
        signature = (atom["source_contig"], tuple(atom["flags"]),
                     tuple(p["object"] for p in atom["before"]), tuple(p["object"] for p in atom["after"]))
        if previous is None or signature != previous[0] or atom["source_start"] != previous[1] + 1:
            flush()
            group = []
        if atom["flags"]:
            group.append(atom)
        previous = signature, atom["source_end"]
    flush()

    # Connectivity groups describe one-to-many splits and many-to-one merges.
    neighbors = defaultdict(set)
    for b, a in overlaps:
        neighbors[("before", b)].add(("after", a))
        neighbors[("after", a)].add(("before", b))
    visited = set()
    for node in sorted(neighbors):
        if node in visited:
            continue
        todo, cluster = [node], set()
        while todo:
            n = todo.pop()
            if n not in cluster:
                cluster.add(n)
                todo.extend(neighbors[n] - cluster)
        visited.update(cluster)
        bs = sorted(n for v, n in cluster if v == "before")
        aas = sorted(n for v, n in cluster if v == "after")
        kind = "scaffold_split" if len(bs) == 1 and len(aas) > 1 else "scaffold_merge" if len(bs) > 1 and len(aas) == 1 else "complex_scaffold_reorganization" if len(bs) > 1 and len(aas) > 1 else None
        if kind:
            ids = [x["segment_id"] for x in atoms if x["mapping_status"] == "unique" and x["before"][0]["object"] in bs and x["after"][0]["object"] in aas]
            event(kind, ids, bs, aas, "Connected group of shared unique source intervals; does not imply a single operation")
    # Different component boundaries can occur even when sequence continuity is retained.
    boundaries = []
    for version, assembly in (("before", before), ("after", after)):
        boundaries.append({(c.contig, x) for c in assembly.components for x in (c.source_start, c.source_end)})
    boundary_changes = []
    for direction, points in (("added", boundaries[1] - boundaries[0]), ("removed", boundaries[0] - boundaries[1])):
        for contig, point in sorted(points):
            adjacent = [x for x in atoms if x["source_contig"] == contig and (x["source_end"] == point or x["source_start"] == point + 1)]
            if len(adjacent) == 2 and all(x["mapping_status"] == "unique" for x in adjacent):
                boundary_changes.append({"change": direction, "source_contig": contig, "boundary_after_bp": point,
                                         "segment_ids": [x["segment_id"] for x in adjacent]})
                event("component_boundary_" + direction, [x["segment_id"] for x in adjacent],
                      sorted({p["object"] for x in adjacent for p in x["before"]}),
                      sorted({p["object"] for x in adjacent for p in x["after"]}),
                      f"AGP component boundary {direction} after source base {point}; not necessarily a sequence break")

    bj, bskip = junctions(before, by_before)
    aj, askip = junctions(after, by_after)
    junction_changes = []
    for key in sorted(bj.keys() | aj.keys()):
        bitems, aitems = bj.get(key, []), aj.get(key, [])
        if max(len(bitems), len(aitems)) > 1:
            qc.append({"severity": "review", "code": "ambiguous_junction", "segment_id": "", "detail": str(key)})
            continue
        kind = "junction_added" if not bitems else "junction_removed" if not aitems else (
            "gap_changed" if (bitems[0]["gap_bp"], bitems[0]["gap_signature"]) != (aitems[0]["gap_bp"], aitems[0]["gap_signature"]) else None)
        if kind:
            row = event(kind, [key[0][0], key[1][0]], [x["object"] for x in bitems], [x["object"] for x in aitems],
                        "Canonical source-end adjacency difference; gap metadata compared separately")
            junction_changes.append({"event_id": row["event_id"], "type": kind, "endpoints": key,
                                     "before": bitems, "after": aitems})
    if bskip or askip:
        qc.append({"severity": "review", "code": "junctions_not_auditable", "segment_id": "",
                   "detail": f"before={bskip}, after={askip}; duplicate source or unknown orientation"})
    for obj in before.lengths:
        if obj not in anchors:
            qc.append({"severity": "review", "code": "no_reciprocal_anchor", "segment_id": "", "detail": obj})

    summary = {"before_objects": len(before.lengths), "after_objects": len(after.lengths),
               "before_component_bp": sum(c.end - c.start for c in before.components),
               "after_component_bp": sum(c.end - c.start for c in after.components),
               "before_gap_bp": sum(g[1] - g[0] for gs in before.gaps.values() for g in gs),
               "after_gap_bp": sum(g[1] - g[0] for gs in after.gaps.values() for g in gs),
               "source_union_bp": sum(x["length_bp"] for x in atoms),
               "atomic_intervals": len(atoms), "review_items": sum(x["severity"] == "review" for x in qc)}
    for status in ("unique", "ambiguous", "before_only", "after_only"):
        summary[status + "_bp"] = sum(x["length_bp"] for x in atoms if x["mapping_status"] == status)
    for flag in ("cross_scaffold", "orientation_changed", "relative_order_changed", "unanchored_scaffold"):
        summary[flag + "_bp"] = sum(x["length_bp"] for x in atoms if flag in x["flags"])
    summary["changed_unique_bp"] = sum(x["length_bp"] for x in atoms if x["mapping_status"] == "unique" and set(x["flags"]) & {"cross_scaffold", "orientation_changed", "relative_order_changed", "unanchored_scaffold"})
    for kind in ("junction_added", "junction_removed", "gap_changed"):
        summary[kind] = sum(x["type"] == kind for x in junction_changes)
    for kind in ("whole_scaffold_reverse", "equivalent_scaffold_rename", "scaffold_order_changed", "scaffold_split", "scaffold_merge", "complex_scaffold_reorganization"):
        summary[kind] = sum(x["type"] == kind for x in events)
    scaffold_map = [{"before_object": b, "after_object": a, "shared_unique_bp": bp,
                     "before_index": before_rank[b], "after_index": after_rank[a],
                     "relative_order_changed": b in reordered_scaffolds if anchors.get(b) == a else False,
                     "reciprocal_anchor": anchors.get(b) == a,
                     "sequence_equivalence": equivalences.get(b, "unanchored") if anchors.get(b) == a else "non_anchor"}
                    for (b, a), bp in sorted(overlaps.items())]
    return {"schema_version": 1, "coordinates": "1-based inclusive; junction left_end/right_start are flanking bases",
            "before_path": str(before.path.resolve()), "after_path": str(after.path.resolve()),
            "before_lengths": before.lengths, "after_lengths": after.lengths,
            "before_gaps": dict(before.gaps), "after_gaps": dict(after.gaps),
            "atoms": atoms, "events": events, "junction_changes": junction_changes,
            "boundary_changes": boundary_changes, "scaffold_map": scaffold_map, "qc": qc, "summary": summary}


def write_tsv(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, ensure_ascii=False, separators=(",", ":")) if isinstance(value, (list, dict, tuple)) else value
                             for key, value in row.items() if key in fields})


OUTPUT_NAMES = ["agp_comparison.json", "segment_correspondence.tsv", "junction_changes.tsv", "curation_events.tsv",
                "component_boundary_changes.tsv", "scaffold_correspondence.tsv", "curation_qc.tsv", "curation_summary.tsv",
                "curation_comparison.pdf", "curation_comparison.svg", "curation_report.html"]


def write_tables(result, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "agp_comparison.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    rows = []
    for atom in result["atoms"]:
        row = dict(atom)
        for version in ("before", "after"):
            row[version + "_occurrences"] = atom[version]
            row[version + "_count"] = len(atom[version])
            for field_name in ("object", "start", "end", "strand"):
                row[version + "_" + field_name] = atom[version][0][field_name] if len(atom[version]) == 1 else ""
        rows.append(row)
    fields = ["segment_id", "source_contig", "source_start", "source_end", "length_bp", "source_hap", "source_method", "mapping_status"]
    fields += [version + "_" + key for version in ("before", "after") for key in ("count", "object", "start", "end", "strand", "occurrences")]
    write_tsv(out_dir / "segment_correspondence.tsv", rows, fields + ["flags", "events", "after_assignment"])
    write_tsv(out_dir / "curation_events.tsv", result["events"], ["event_id", "type", "segment_ids", "before_objects", "after_objects", "length_bp", "interpretation"])
    write_tsv(out_dir / "junction_changes.tsv", result["junction_changes"], ["event_id", "type", "endpoints", "before", "after"])
    write_tsv(out_dir / "component_boundary_changes.tsv", result["boundary_changes"], ["change", "source_contig", "boundary_after_bp", "segment_ids"])
    write_tsv(out_dir / "scaffold_correspondence.tsv", result["scaffold_map"], ["before_object", "after_object", "before_index", "after_index", "relative_order_changed", "shared_unique_bp", "reciprocal_anchor", "sequence_equivalence"])
    write_tsv(out_dir / "curation_qc.tsv", result["qc"], ["severity", "code", "segment_id", "detail"])
    write_tsv(out_dir / "curation_summary.tsv", [{"metric": k, "value": v} for k, v in result["summary"].items()], ["metric", "value"])
