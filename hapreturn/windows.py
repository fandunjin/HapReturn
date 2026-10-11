"""Expand local unresolved windows within their original contiguous gray run."""


def source_label(first, second, length, labels, min_coverage, dominance):
    total, maximum = first + second, max(first, second)
    if (not total or first == second or maximum / length < min_coverage
            or maximum / total < dominance):
        return "unresolved"
    return labels[0] if first > second else labels[1]


def expand_unresolved(rows, labels, min_coverage, dominance):
    """Aggregate original PAF evidence; never borrow from a resolved neighbor.

    Context widths double in whole base-window units, stopping at the full gray
    run. Raw per-window coverage is unchanged. Every attempted context is traced.
    """
    result, annotations, attempts = list(rows), {}, []
    prefix = [[0], [0]]
    for row in rows:
        for index in (0, 1):
            prefix[index].append(prefix[index][-1] + row[3 + index])
        annotations[(row[0], row[1])] = (row[5], "window_evidence" if row[5] != "unresolved" else "unresolved", row[1], row[2])
    start = 0
    while start < len(rows):
        if rows[start][5] != "unresolved":
            start += 1
            continue
        end = start + 1
        while end < len(rows) and rows[end][5] == "unresolved":
            end += 1
        for index in range(start, end):
            width = min(2, end - start)
            while True:
                lo = max(start, min(index - width // 2, end - width))
                hi = lo + width
                first, second = (values[hi] - values[lo] for values in prefix)
                a, b = rows[lo][1], rows[hi - 1][2]
                label = source_label(first, second, b - a, labels, min_coverage, dominance)
                attempts.append((rows[index][0], rows[index][1], rows[index][2], a, b, b - a, first, second, label))
                if label != "unresolved":
                    result[index] = (*rows[index][:5], label)
                    annotations[(rows[index][0], rows[index][1])] = ("unresolved", "expanded_window", a, b)
                    break
                if width == end - start:
                    break
                width = min(width * 2, end - start)
        start = end
    return result, annotations, attempts
