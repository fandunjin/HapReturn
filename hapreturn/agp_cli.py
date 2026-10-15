"""Compare pre/post-curation AGPs without FASTA alignment or Hi-C recomputation."""
import argparse
import csv
import hashlib
import sys
from pathlib import Path

from .agp import OUTPUT_NAMES, compare_agp, read_agp, read_source_map, write_tables
from .agp_provenance import (PROVENANCE_OUTPUTS, attach_provenance, export_partition,
                             read_alignment_assignments, write_provenance, assignments_from_paf)


def read_assignments(path, after):
    if path is None:
        return {}
    result = {}
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"seq_id", "assignment"} <= set(reader.fieldnames or []):
            raise ValueError("Assignment TSV requires seq_id and assignment columns")
        for row in reader:
            name = row["seq_id"]
            if name not in after.lengths or name in result or not row["assignment"]:
                raise ValueError(f"Invalid, absent or duplicate assignment: {name}")
            if "length" in row and int(row["length"]) != after.lengths[name]:
                raise ValueError(f"Assignment/AGP length mismatch: {name}")
            result[name] = row["assignment"]
    return result


def prepare_comparison(before_path, after_path, source_map=None, assignment=None):
    before, after = read_agp(before_path), read_agp(after_path)
    result = compare_agp(before, after, read_source_map(source_map), read_assignments(assignment, after))
    result["input_sha256"] = {"before": hashlib.sha256(Path(before_path).read_bytes()).hexdigest(),
                              "after": hashlib.sha256(Path(after_path).read_bytes()).hexdigest()}
    result["source_label_method"] = "explicit_source_map" if source_map else "h1tg_h2tg_name_convention_else_unresolved"
    attach_provenance(result, alignment_rows=read_alignment_assignments(assignment, after.lengths))
    return result


def generate_comparison(result, out_dir, min_display_bp=1000000, no_plots=False):
    out_dir = Path(out_dir)
    if no_plots and any((out_dir / name).exists() for name in ("curation_report.html", "curation_comparison.pdf", "curation_comparison.svg")):
        raise ValueError("--no-plots would leave stale figures in the output directory; choose a new directory")
    if "agp_provenance" not in result:
        attach_provenance(result)
    write_tables(result, out_dir)
    write_provenance(result, out_dir)
    if not no_plots:
        from .agp_report import plot_pdf, write_html, provenance_pdf
        from .agp_provenance import source_view
        plot_pdf(result, out_dir, min_display_bp)
        view = source_view(result)
        if view["atoms"]:
            plot_pdf(view, out_dir, 0, filename="agp_source_correspondence.pdf", title="Original hifiasm components to current chromosomes",
                     left_title="Original contigs (observed source axes)", right_title="Current chromosomes")
        provenance_pdf(result, out_dir)
        write_html(result, out_dir, min_display_bp)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before-agp", type=Path, required=True, help="Pre-curation AGP using original component IDs")
    parser.add_argument("--after-agp", type=Path, required=True, help="Post-curation AGP using the same source coordinates")
    parser.add_argument("--source-map", type=Path, help="Optional TSV: contig_id, source; overrides h1tg/h2tg name inference")
    parser.add_argument("--assignment", "--minimap-assignment", dest="assignment", type=Path,
                        help="Optional minimap2-based HapReturn assignment.tsv; kept separate from AGP assignments")
    parser.add_argument("--final-fasta", type=Path, help="Export AGP-based chromosome partitions from this post-curation FASTA")
    parser.add_argument("--chromosome-count", type=int, help="Select first N post-curation AGP records; all others become unplaced")
    parser.add_argument("--chromosome-list", type=Path, help="Selected post-curation chromosome IDs, one per line")
    parser.add_argument("--source-labels", nargs=2, default=["hap1", "hap2"])
    parser.add_argument("--agp-dominance", type=float, default=0.60)
    parser.add_argument("--agp-min-known-fraction", type=float, default=0.25)
    parser.add_argument("--agp-assignment-mode", choices=["provenance", "balanced"], default="provenance")
    parser.add_argument("--homolog-groups", type=Path, help="Explicit group_id/seq_id TSV; required for balanced AGP assignment")
    parser.add_argument("--ploidy", type=int, help="Even group size for balanced AGP assignment")
    parser.add_argument("--hap1-paf", type=Path, help="Optional original hap1-to-current PAF for independent comparison")
    parser.add_argument("--hap2-paf", type=Path, help="Optional original hap2-to-current PAF for independent comparison")
    parser.add_argument("--min-mapq", type=int, default=10)
    parser.add_argument("--min-aln-len", type=int, default=10000)
    parser.add_argument("--min-covered-fraction", type=float, default=0.25)
    parser.add_argument("--minimap-dominance", type=float, default=0.60)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--min-display-bp", type=int, default=1000000, help="Overview track cutoff only; all intervals remain in tables, HTML selector and PDF appendix")
    parser.add_argument("--no-plots", action="store_true", help="Write JSON/TSV only")
    parser.add_argument("--force", action="store_true", help="Reuse a nonempty output directory")
    parser.add_argument("--allow-review", action="store_true", help="Exit 0 instead of 2 when ambiguity/unknown orientation requires review")
    args = parser.parse_args(argv)
    try:
        if args.min_display_bp < 0:
            raise ValueError("min-display-bp must be nonnegative")
        if args.chromosome_count is not None and args.chromosome_list:
            raise ValueError("Use chromosome-count or chromosome-list, not both")
        if bool(args.hap1_paf) != bool(args.hap2_paf) or (args.assignment and args.hap1_paf):
            raise ValueError("Provide both PAFs, or one minimap assignment TSV; do not mix them")
        inputs = {p.resolve() for p in (args.before_agp, args.after_agp, args.source_map, args.assignment,
                                       args.final_fasta, args.chromosome_list, args.homolog_groups, args.hap1_paf, args.hap2_paf) if p is not None}
        outputs = OUTPUT_NAMES + PROVENANCE_OUTPUTS + ["agp_partition/"+label+".fa" for label in args.source_labels+["unresolved","unplaced"]]
        if any((args.out_dir / name).resolve() in inputs for name in outputs):
            raise ValueError("Output would overwrite an input")
        if args.out_dir.exists() and (not args.out_dir.is_dir() or (any(args.out_dir.iterdir()) and not args.force)):
            raise ValueError("Output directory is nonempty; choose a new directory or use --force")
        result = prepare_comparison(args.before_agp, args.after_agp, args.source_map, args.assignment)
        selected = list(result["after_lengths"])
        if args.chromosome_count is not None:
            if not 1 <= args.chromosome_count <= len(selected):
                raise ValueError("chromosome-count must be between 1 and post-curation record count")
            selected = selected[:args.chromosome_count]
        if args.chromosome_list:
            selected = [line.strip() for line in args.chromosome_list.read_text(encoding="utf-8-sig").splitlines()
                        if line.strip() and not line.startswith("#")]
        if not selected or len(set(selected)) != len(selected) or set(selected) - result["after_lengths"].keys():
            raise ValueError("Selected chromosome IDs are empty, duplicated or absent from AGP")
        from .cli import load_groups, fasta_lengths
        groups = load_groups(args.homolog_groups, {n: result["after_lengths"][n] for n in selected}, args.ploidy)
        alignment_rows = (assignments_from_paf((args.hap1_paf,args.hap2_paf), result["after_lengths"], args.source_labels,
                                              args.min_mapq,args.min_aln_len,args.min_covered_fraction,args.minimap_dominance)
                          if args.hap1_paf else read_alignment_assignments(args.assignment, result["after_lengths"]))
        result["alignment_comparison"] = {
            "method": "paf_union_coverage" if args.hap1_paf else "imported_assignment" if args.assignment else "not_available",
            "assignment_path": str(args.assignment.resolve()) if args.assignment else "",
            "paf_paths": [str(p.resolve()) for p in (args.hap1_paf,args.hap2_paf)] if args.hap1_paf else [],
            "paf_filters": {"min_mapq":args.min_mapq,"min_aln_len":args.min_aln_len,
                            "min_coverage":args.min_covered_fraction,"dominance":args.minimap_dominance} if args.hap1_paf else {}}
        attach_provenance(result, args.source_labels, selected, args.agp_dominance, args.agp_min_known_fraction,
                          args.agp_assignment_mode, groups, args.ploidy,
                          alignment_rows)
        if args.final_fasta and fasta_lengths(args.final_fasta) != result["after_lengths"]:
            raise ValueError("Final FASTA IDs/lengths disagree with post-curation AGP")
        if args.final_fasta:
            result["sequence_export"] = {"input": str(args.final_fasta.resolve()),
                                         "files": ["agp_partition/"+l+".fa" for l in args.source_labels+["unresolved","unplaced"]]}
        generate_comparison(result, args.out_dir, args.min_display_bp, args.no_plots)
        if args.final_fasta:
            export_partition(result, args.final_fasta, args.out_dir)
        s = result["summary"]
        print(f"AGP comparison: unique_bp={s['unique_bp']} whole_reverse={s['whole_scaffold_reverse']} "
              f"scaffold_order_changed={s['scaffold_order_changed']} internal_or_membership_changed_bp={s['changed_unique_bp']} "
              f"review_items={s['review_items']} agp_review_chromosomes={result['agp_provenance']['review_chromosomes']}\nOutputs: {args.out_dir}")
        return 2 if (s["review_items"] or result["agp_provenance"]["review_chromosomes"]) and not args.allow_review else 0
    except (OSError, ValueError, KeyError) as exc:
        print(f"hapreturn-agp: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
