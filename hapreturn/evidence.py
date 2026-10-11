"""Validated PAF records shared by assignment and provenance figures."""
from dataclasses import dataclass


@dataclass(frozen=True)
class Alignment:
    source: str
    query: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target: str
    target_length: int
    start: int
    end: int


def read_paf(path, records, min_mapq=10, min_aln_len=10000, source=""):
    query_lengths = {}
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF {path}:{number}")
            try:
                query_length, query_start, query_end = map(int, fields[1:4])
                target_length, start, end = map(int, fields[6:9])
                block_length, mapq = int(fields[10]), int(fields[11])
            except ValueError as exc:
                raise ValueError(f"Invalid PAF coordinates {path}:{number}") from exc
            query, target, strand = fields[0], fields[5], fields[4]
            if target not in records:
                raise ValueError(f"PAF target {target} absent from final FASTA")
            expected = records[target] if isinstance(records[target], int) else len(records[target])
            if target_length != expected or not 0 <= start < end <= target_length:
                raise ValueError(f"PAF target coordinates/length disagree at {path}:{number}")
            if not query or not 0 <= query_start < query_end <= query_length or strand not in ("+", "-"):
                raise ValueError(f"PAF query coordinates/strand disagree at {path}:{number}")
            if query in query_lengths and query_lengths[query] != query_length:
                raise ValueError(f"Conflicting PAF query lengths at {path}:{number}")
            query_lengths[query] = query_length
            if "tp:A:S" in fields[12:] or mapq < min_mapq or block_length < min_aln_len:
                continue
            yield Alignment(source, query, query_length, query_start, query_end,
                            strand, target, target_length, start, end)
