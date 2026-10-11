"""Exact two-source optimization within equal-copy homolog groups."""


def consecutive_groups(records, ploidy):
    if ploidy is None or ploidy < 2 or ploidy % 2:
        raise ValueError("Balanced assignment requires an even --ploidy >=2")
    names = list(records)
    if len(names) % ploidy:
        raise ValueError("Selected chromosome count must be divisible by ploidy")
    return {f"group_{i // ploidy + 1:03d}": names[i:i + ploidy]
            for i in range(0, len(names), ploidy)}


def balanced_assignment(groups, support, labels, ploidy):
    """Maximize sum C_i,source with exactly ploidy/2 records per source.

    Relative to assigning all records to source 2, selecting a record for source
    1 changes the objective by C_i,1-C_i,2. The largest k changes are optimal.
    Equal changes at the cutoff yield multiple optima; use input order and audit.
    """
    if ploidy is None or ploidy < 2 or ploidy % 2:
        raise ValueError("Balanced assignment requires an even --ploidy >=2")
    first, second = labels
    assignments, audit = {}, {}
    for group, members in groups.items():
        if len(members) != ploidy:
            raise ValueError(f"Group {group} must contain exactly {ploidy} chromosomes")
        ordered = sorted(members, key=lambda name: support[name][first] - support[name][second], reverse=True)
        k = ploidy // 2
        chosen = set(ordered[:k])
        for name in members:
            assignments[name] = first if name in chosen else second
        gap = ((support[ordered[k - 1]][first] - support[ordered[k - 1]][second])
               - (support[ordered[k]][first] - support[ordered[k]][second]))
        audit[group] = {"total_support_bp": sum(support[n][assignments[n]] for n in members),
                        "optimality_gap_bp": gap, "tied_optimum": gap == 0}
    return assignments, audit
