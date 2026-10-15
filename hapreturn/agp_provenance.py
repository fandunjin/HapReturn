"""Trace actual AGP component provenance independently of alignment similarity."""
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

from .agp import write_tsv


PROVENANCE_OUTPUTS = ["agp_chromosome_provenance.tsv", "agp_chromosome_components.tsv",
                      "agp_original_scaffold_composition.tsv", "agp_assignment.tsv", "agp_partition_summary.tsv",
                      "agp_group_qc.tsv", "agp_vs_minimap.tsv", "agp_source_correspondence.pdf",
                      "agp_source_correspondence.svg", "agp_provenance.pdf"]


def read_alignment_assignments(path, lengths):
    if path is None:
        return {}
    result = {}
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"seq_id", "assignment"} <= set(reader.fieldnames or []):
            raise ValueError("Minimap assignment requires seq_id and assignment columns")
        for row in reader:
            name = row["seq_id"]
            if name not in lengths or name in result or not row["assignment"]:
                raise ValueError(f"Invalid, missing or duplicate alignment-assigned object: {name}")
            if "length" in row and int(row["length"]) != lengths[name]:
                raise ValueError(f"Alignment assignment/AGP length mismatch: {name}")
            result[name] = dict(row)
    return result


def assignments_from_paf(paths, lengths, labels=("hap1", "hap2"), min_mapq=10, min_aln_len=10000,
                         min_coverage=0.25, dominance=0.60):
    """Use the same independent union-coverage rule as the minimap2 pipeline."""
    from .evidence import read_paf
    from .cli import union_length
    if min_mapq < 0 or min_aln_len < 1 or not 0 <= min_coverage <= 1 or not .5 <= dominance <= 1:
        raise ValueError("Invalid PAF comparison filters")
    spans = {label: defaultdict(list) for label in labels}
    for label, path in zip(labels, paths):
        for hit in read_paf(path, lengths, min_mapq, min_aln_len, label):
            spans[label][hit.target].append((hit.start, hit.end))
    result = {}
    first, second = labels
    for name, length in lengths.items():
        bp = {label: union_length(spans[label][name]) for label in labels}
        total, maximum = sum(bp.values()), max(bp.values())
        share, coverage = maximum/total if total else 0, maximum/length
        label = first if bp[first] > bp[second] else second
        if bp[first] == bp[second] or coverage < min_coverage or share < dominance:
            label = "unresolved"
        result[name] = {"seq_id": name, "length": length, **{l+"_covered_bp":bp[l] for l in labels},
                        "dominance": share, "dominant_coverage_fraction": coverage,
                        "assignment": label, "evidence_assignment": label, "assignment_mode":"provenance"}
    return result


def attach_provenance(result, labels=("hap1", "hap2"), selected=None, dominance=0.60,
                      min_known_fraction=0.25, mode="provenance", groups=None, ploidy=None,
                      alignment_rows=None):
    """Assign whole records while preserving mixed local sources and unknown bases.

    AGP support counts target sequence bases actually built from each source;
    gaps are excluded. Repeated source placements still occupy real target bases,
    but are retained as QC rather than mistaken for competing alignment support.
    """
    labels = list(labels)
    reserved = {"unresolved", "unplaced", "con", "prn", "aux", "nul"} | {f"{p}{i}" for p in ("com", "lpt") for i in range(1,10)}
    if len(labels) != 2 or len({x.casefold() for x in labels}) != 2 or any(not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*",x) or x.casefold() in reserved for x in labels):
        raise ValueError("AGP provenance requires two distinct nonreserved source labels")
    if not 0.5 <= dominance <= 1 or not 0 <= min_known_fraction <= 1:
        raise ValueError("AGP dominance/known-fraction threshold out of range")
    if mode not in {"provenance", "balanced"}:
        raise ValueError("Unknown AGP assignment mode")
    names = list(result["after_lengths"]) if selected is None else list(selected)
    if not names or len(names) != len(set(names)) or set(names) - result["after_lengths"].keys():
        raise ValueError("Empty, duplicate or absent selected AGP chromosomes")
    selected_set = set(names)
    groups = groups or {}
    if groups:
        members = [n for ns in groups.values() for n in ns]
        if len(members) != len(set(members)) or set(members) != selected_set:
            raise ValueError("AGP homolog groups must list every selected chromosome exactly once")
    if mode == "balanced" and (not groups or ploidy is None or ploidy < 2 or ploidy % 2):
        raise ValueError("Balanced AGP assignment requires explicit homolog groups and even ploidy")
    group_of = {n: g for g, ns in groups.items() for n in ns}
    components, by_chromosome = [], defaultdict(list)
    source_axes = defaultdict(int)
    for atom in result["atoms"]:
        source_axes[atom["source_contig"]] = max(source_axes[atom["source_contig"]], atom["source_end"])
        for occurrence, after in enumerate(atom["after"], 1):
            before = atom["before"][0] if len(atom["before"]) == 1 else None
            row = {"seq_id": after["object"], "start": after["start"], "end": after["end"],
                   "strand": after["strand"], "length_bp": atom["length_bp"], "source_hap": atom["source_hap"],
                   "source_method": atom["source_method"], "source_contig": atom["source_contig"],
                   "source_start": atom["source_start"], "source_end": atom["source_end"],
                   "original_scaffold": before["object"] if before else "",
                   "before_start": before["start"] if before else "", "before_end": before["end"] if before else "",
                   "before_strand": before["strand"] if before else "", "before_occurrences": atom["before"],
                   "mapping_status": atom["mapping_status"], "segment_id": atom["segment_id"], "after_occurrence": occurrence}
            components.append(row)
            by_chromosome[after["object"]].append(row)
    for rows in by_chromosome.values():
        rows.sort(key=lambda row: row["start"])
    object_rank = {name: i for i,name in enumerate(result["after_lengths"])}
    components.sort(key=lambda row: (object_rank[row["seq_id"]],row["start"]))
    rows, support, independent = [], {}, {}
    first, second = labels
    for name, length in result["after_lengths"].items():
        parts = by_chromosome[name]
        bp = Counter()
        for p in parts:
            bp[p["source_hap"] if p["source_hap"] in labels else "unresolved"] += p["length_bp"]
        seq_bp = sum(bp.values())
        known = bp[first] + bp[second]
        maximum = max(bp[first], bp[second])
        share = maximum / known if known else 0
        fraction = known / seq_bp if seq_bp else 0
        winner = first if bp[first] > bp[second] else second
        reason = "source_majority"
        if not known or fraction < min_known_fraction:
            winner, reason = "unresolved", "insufficient_known_source"
        elif bp[first] == bp[second] or share < dominance:
            winner, reason = "unresolved", "mixed_source_below_dominance"
        if not parts:
            winner, reason = "unresolved", "no_sequence_component"
        independent[name] = winner
        support[name] = {first: bp[first], second: bp[second]}
        status = "mixed" if bp[first] and bp[second] else "pure" if known and not bp["unresolved"] else "partly_unknown" if known else "unknown"
        related_qc = [q["code"] for q in result["qc"] if q["segment_id"] in {p["segment_id"] for p in parts}]
        rows.append({"seq_id": name, "length": length, "sequence_bp": seq_bp, "gap_bp": length-seq_bp,
                     first + "_bp": bp[first], second + "_bp": bp[second], "unknown_source_bp": bp["unresolved"],
                     first + "_fraction": bp[first]/seq_bp if seq_bp else 0,
                     second + "_fraction": bp[second]/seq_bp if seq_bp else 0,
                     "known_source_fraction": fraction, "source_dominance": share,
                     "composition_status": status, "source_contig_count": len({p["source_contig"] for p in parts}),
                     "original_scaffolds": sorted({p["original_scaffold"] for p in parts if p["original_scaffold"]}),
                     "independent_assignment": winner, "independent_reason": reason,
                     "assignment": winner if name in selected_set else "unplaced", "selected": name in selected_set,
                     "assignment_mode": mode, "group_id": group_of.get(name, ""), "review_flags": sorted(set(related_qc))})
    audit, constrained = {}, {}
    if mode == "balanced":
        from .assignment import balanced_assignment
        constrained, audit = balanced_assignment(groups, support, labels, ploidy)
    for row in rows:
        name = row["seq_id"]
        if name in constrained:
            row["assignment"] = constrained[name]
        if row["selected"] and row["independent_assignment"] == "unresolved":
            row["review_flags"].append(row["independent_reason"])
        if row["selected"] and row["assignment"] != row["independent_assignment"]:
            row["review_flags"].append("constraint_changed_source_majority")
        if audit.get(row["group_id"], {}).get("tied_optimum"):
            row["review_flags"].append("tied_group_optimum")
    group_rows = []
    for group, ns in groups.items():
        chosen = [row for row in rows if row["seq_id"] in ns]
        flags = sorted({f for row in chosen for f in row["review_flags"]})
        group_rows.append({"group_id": group, "n_records": len(ns),
                           first: sum(r["assignment"] == first for r in chosen),
                           second: sum(r["assignment"] == second for r in chosen),
                           "review_flags": flags, "status": "review" if flags else "pass", **audit.get(group, {})})
    composition = defaultdict(lambda: {"length_bp": 0, "segments": []})
    for row in components:
        origin = row["original_scaffold"] or ("[ambiguous]" if len(row["before_occurrences"]) > 1 else "[unmatched]")
        key = row["seq_id"], origin, row["source_hap"]
        composition[key]["length_bp"] += row["length_bp"]
        composition[key]["segments"].append(row["segment_id"])
    origin_rows = [{"seq_id": key[0], "original_scaffold": key[1], "source_hap": key[2], **value}
                   for key, value in composition.items()]
    partition = [{"source": label, "n_sequences": sum(r["assignment"] == label for r in rows),
                  "total_bp": sum(r["length"] for r in rows if r["assignment"] == label)}
                 for label in labels + ["unresolved", "unplaced"]]
    alignment_rows = alignment_rows or {}
    comparisons = []
    for row in rows:
        name = row["seq_id"]
        alignment = alignment_rows.get(name)
        mm = alignment.get("assignment", "") if alignment else ""
        mm_independent = alignment.get("evidence_assignment", mm) if alignment else ""
        comparable = bool(alignment) and row["selected"]
        status = ("not_selected" if not row["selected"] else "not_available" if not alignment else
                  "label_mismatch" if mm not in labels + ["unresolved"] else
                  "agree" if mm == row["assignment"] else
                  "one_unresolved" if "unresolved" in (mm, row["assignment"]) else "different_source")
        comparisons.append({"seq_id": name, "selected": row["selected"],
                            "agp_assignment": row["assignment"], "agp_independent_assignment": row["independent_assignment"],
                            "minimap_assignment": mm, "minimap_independent_assignment": mm_independent,
                            "comparison_status": status,
                            "independent_labels_agree": row["independent_assignment"] == mm_independent if comparable else "",
                            "agp_"+first+"_bp": row[first+"_bp"], "agp_"+second+"_bp": row[second+"_bp"],
                            "minimap_"+first+"_covered_bp": alignment.get(first+"_covered_bp", "") if alignment else "",
                            "minimap_"+second+"_covered_bp": alignment.get(second+"_covered_bp", "") if alignment else "",
                            "minimap_assignment_mode": alignment.get("assignment_mode", "") if alignment else ""})
    provenance = {"labels": labels, "selected": names, "dominance": dominance, "min_known_fraction": min_known_fraction,
                  "assignment_mode": mode, "groups": groups, "ploidy": ploidy,
                  "chromosomes": rows, "components": components, "original_scaffold_composition": origin_rows,
                  "partition_summary": partition, "group_qc": group_rows, "method_comparison": comparisons,
                  "alignment_rows": alignment_rows, "source_axes": dict(source_axes),
                  "source_axis_contract": "Max original coordinate represented in either AGP; not necessarily full original contig length",
                  "review_chromosomes": sum(bool(r["review_flags"]) for r in rows if r["selected"])}
    result["agp_provenance"] = provenance
    assignment_by_name = {r["seq_id"]: r["assignment"] for r in rows}
    for atom in result["atoms"]:
        atom["agp_assignment"] = assignment_by_name[atom["after"][0]["object"]] if len(atom["after"]) == 1 else ""
    return provenance


def source_view(result):
    """A graphical source-coordinate view, independent of pre-curation scaffold names."""
    p = result["agp_provenance"]
    selected = set(p["selected"])
    atoms, sources = [], set()
    for atom in result["atoms"]:
        after = [x for x in atom["after"] if x["object"] in selected]
        if not after:
            continue
        sources.add(atom["source_contig"])
        atoms.append(dict(atom, before=[{"object": atom["source_contig"], "start": atom["source_start"],
                                        "end": atom["source_end"], "strand": "+", "part": 1}],
                          after=after, mapping_status="unique" if len(after) == 1 else "ambiguous"))
    return {"before_lengths": {name: p["source_axes"][name] for name in sorted(sources)},
            "after_lengths": {name: result["after_lengths"][name] for name in p["selected"]},
            "atoms": atoms, "view_kind": "hifiasm_source", "source_labels": p["labels"]}


def write_provenance(result, out_dir):
    p = result["agp_provenance"]
    labels = p["labels"]
    fields = ["seq_id", "length", "sequence_bp", "gap_bp", *(x+"_bp" for x in labels), "unknown_source_bp",
              *(x+"_fraction" for x in labels), "known_source_fraction", "source_dominance", "composition_status",
              "source_contig_count", "original_scaffolds", "independent_assignment", "independent_reason", "assignment",
              "selected", "assignment_mode", "group_id", "review_flags"]
    write_tsv(Path(out_dir)/"agp_chromosome_provenance.tsv", p["chromosomes"], fields)
    write_tsv(Path(out_dir)/"agp_assignment.tsv", p["chromosomes"], fields)
    write_tsv(Path(out_dir)/"agp_chromosome_components.tsv", p["components"],
              ["seq_id", "start", "end", "strand", "source_hap", "source_method", "source_contig", "source_start", "source_end",
               "length_bp", "original_scaffold", "before_start", "before_end", "before_strand", "before_occurrences",
               "mapping_status", "segment_id", "after_occurrence"])
    write_tsv(Path(out_dir)/"agp_original_scaffold_composition.tsv", p["original_scaffold_composition"],
              ["seq_id", "original_scaffold", "source_hap", "length_bp", "segments"])
    write_tsv(Path(out_dir)/"agp_partition_summary.tsv", p["partition_summary"], ["source", "n_sequences", "total_bp"])
    write_tsv(Path(out_dir)/"agp_group_qc.tsv", p["group_qc"],
              ["group_id", "n_records", *labels, "status", "review_flags", "total_support_bp", "optimality_gap_bp", "tied_optimum"])
    compare_fields = ["seq_id", "selected", "agp_assignment", "agp_independent_assignment", "minimap_assignment",
                      "minimap_independent_assignment", "comparison_status", "independent_labels_agree",
                      *("agp_"+x+"_bp" for x in labels), *("minimap_"+x+"_covered_bp" for x in labels), "minimap_assignment_mode"]
    write_tsv(Path(out_dir)/"agp_vs_minimap.tsv", p["method_comparison"], compare_fields)


def export_partition(result, final_fasta, out_dir):
    from .cli import fasta_lengths, write_partitioned_fasta
    lengths = fasta_lengths(final_fasta)
    if lengths != result["after_lengths"]:
        raise ValueError("Final FASTA IDs/lengths must match post-curation AGP for sequence export")
    p = result["agp_provenance"]
    assignments = {r["seq_id"]: r["assignment"] for r in p["chromosomes"]}
    partition_dir = Path(out_dir)/"agp_partition"
    partition_dir.mkdir(parents=True, exist_ok=True)
    write_partitioned_fasta(final_fasta, partition_dir, p["labels"], assignments)
    return partition_dir
