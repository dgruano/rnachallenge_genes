import gzip
import re
from typing import Callable, Dict, Iterable, Optional

import pandas as pd

RESOLVED_COLS = [
    "transcript_id",
    "db_source",
    "gene_id",
    "gene_symbol",
    "organism",
    "assembly_name",       # human-readable id (R64-5-1, GRCh38, NC_000001.11…)
    "assembly_accession",  # GCF_/GCA_ if known at Stage 1; pd.NA otherwise
    "fasta_url",           # HTTPS URL if in config; pd.NA otherwise
    "gtf_url",             # HTTPS URL if in config; pd.NA otherwise
    "gtf_format",          # "gtf" | "gff3" | pd.NA
    "chrom",
    "start",
    "end",
    "strand",
    "is_ambiguous",
]

UNRESOLVED_COLS = ["transcript_id", "raw_header", "source_file", "db_source", "reason"]

_GTF_ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def open_annotation(path: str):
    opener = gzip.open if str(path).endswith(".gz") else open
    return opener(path, "rt", encoding="utf-8", errors="replace")


def parse_annotation_attrs(attr_str: str) -> Dict[str, str]:
    if "=" in attr_str and '"' not in attr_str:
        attrs: Dict[str, str] = {}
        for field in attr_str.split(";"):
            field = field.strip()
            if "=" in field:
                key, _, val = field.partition("=")
                attrs[key.strip()] = val.strip()
        return attrs
    return dict(_GTF_ATTR_RE.findall(attr_str))


def normalize_strand(value) -> str:
    if value in (1, "1", "+", "+1"):
        return "+"
    if value in (-1, "-1", "-", "-1"):
        return "-"
    return "."


def generic_candidates(tid: str) -> list[str]:
    seeds = [tid]
    for sep in ("|", ":"):
        if sep in tid:
            seeds.append(tid.split(sep)[-1])

    candidates = []
    for seed in seeds:
        candidates.append(seed)
        if seed.endswith("_cdna"):
            candidates.append(seed[: -len("_cdna")])
        candidates.append(re.sub(r"_\d+_cdna$", "", seed))
        candidates.append(re.sub(r"_\d+$", "", seed))
        if "_" in seed:
            candidates.append(seed.split("_")[0])
        if re.search(r"\.\d+(?:\.\d+)?$", seed):
            parts = seed.split(".")
            while len(parts) > 1:
                parts = parts[:-1]
                candidates.append(".".join(parts))
    return [c for c in dict.fromkeys(candidates) if c]


def wormbase_candidates(tid: str) -> list[str]:
    candidates = generic_candidates(tid)
    candidates.append(re.sub(r"\.[a-z]?\d+$", "", tid))
    candidates.append(re.sub(r"\.[a-z]\.\d+$", "", tid))
    return [c for c in dict.fromkeys(candidates) if c]


_SGD_SOURCE_RE = re.compile(r"^Source:SGD;Acc:(S\d+)$", re.IGNORECASE)


def yeast_candidates(tid: str) -> list[str]:
    candidates = [tid]

    # ROI #6: canonicalize Source:SGD;Acc:S000028522 → the bare SGDID and its
    # dbxref form (SGD:S000028522), the keys the GFF3 index exposes.
    m = _SGD_SOURCE_RE.match(tid)
    if m:
        sgdid = m.group(1)
        candidates.append(sgdid)
        candidates.append(f"SGD:{sgdid}")

    stripped = re.sub(r"_[AB]$", "", tid)
    if stripped != tid:
        candidates.append(stripped)

    if "_" in tid:
        dash_var = tid.replace("_", "-")
        candidates.append(dash_var)
        stripped_dash = re.sub(r"-[AB]$", "", dash_var)
        if stripped_dash != dash_var:
            candidates.append(stripped_dash)

    no_ver = re.sub(r"_\d+$", "", tid)
    if no_ver != tid:
        candidates.append(no_ver)

    lower = tid.lower()
    for suffix in ("_cdna", "_mrna"):
        if lower.endswith(suffix):
            candidates.append(tid[: -len(suffix)])
            break

    return [c for c in dict.fromkeys(candidates) if c]


def build_annotation_index(
    annotation_path: str,
    *,
    feature_types: Iterable[str],
    transcript_fields: Iterable[str],
    gene_id_fields: Iterable[str],
    gene_symbol_fields: Iterable[str] = (),
    alias_fields: Iterable[str] = (),
    log=None,
) -> Dict[str, dict]:
    index: Dict[str, dict] = {}
    total_features = 0
    feature_types = set(feature_types)

    with open_annotation(annotation_path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith(">") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in feature_types:
                continue

            attrs = parse_annotation_attrs(parts[8])
            gene_id = next((attrs.get(field, "").strip() for field in gene_id_fields if attrs.get(field, "").strip()), "")
            gene_symbol = next(
                (attrs.get(field, "").strip() for field in gene_symbol_fields if attrs.get(field, "").strip()),
                gene_id,
            ) or gene_id

            record = {
                "gene_id": gene_id,
                "gene_symbol": gene_symbol,
                "chrom": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": normalize_strand(parts[6]),
            }

            keys = []
            for field in tuple(transcript_fields) + tuple(gene_id_fields) + tuple(gene_symbol_fields):
                value = attrs.get(field, "").strip()
                if value:
                    keys.append(value)

            for field in alias_fields:
                value = attrs.get(field, "").strip()
                if value:
                    keys.extend(v.strip() for v in value.split(",") if v.strip())

            for key in dict.fromkeys(keys):
                if key not in index:
                    index[key] = record
            total_features += 1

    if log:
        log.info(
            f"Parsed {total_features} features from {annotation_path}; index contains {len(index)} keys"
        )
    return index


def build_resolved_row(
    transcript_id: str,
    *,
    db_source: str,
    organism: str,
    assembly_name: str,
    hit: dict,
    assembly_accession=None,
    fasta_url=None,
    gtf_url=None,
    gtf_format=None,
) -> dict:
    return {
        "transcript_id": transcript_id,
        "db_source": db_source,
        "gene_id": hit.get("gene_id", ""),
        "gene_symbol": hit.get("gene_symbol", "") or hit.get("gene_id", ""),
        "organism": organism,
        "assembly_name": assembly_name,
        "assembly_accession": assembly_accession if assembly_accession is not None else pd.NA,
        "fasta_url": fasta_url if fasta_url is not None else pd.NA,
        "gtf_url": gtf_url if gtf_url is not None else pd.NA,
        "gtf_format": gtf_format if gtf_format is not None else pd.NA,
        "chrom": hit.get("chrom", ""),
        "start": hit.get("start", ""),
        "end": hit.get("end", ""),
        "strand": hit.get("strand", "."),
        "is_ambiguous": False,
    }


def build_unresolved_row(row: pd.Series, *, db_source: str, reason: str) -> dict:
    return {
        "transcript_id": str(row.get("transcript_id", "")),
        "raw_header": str(row.get("raw_header", "")),
        "source_file": str(row.get("source_file", "")),
        "db_source": db_source,
        "reason": reason,
    }


def resolve_classified_ids(
    classified_df: pd.DataFrame,
    *,
    db_source: str,
    organism: str,
    assembly_name: str,
    index: Dict[str, dict],
    candidates_fn: Callable[[str], list[str]],
    unresolved_reason: str,
    pre_resolve_fn: Optional[Callable[[pd.Series], Optional[dict]]] = None,
    assembly_accession=None,
    fasta_url=None,
    gtf_url=None,
    gtf_format=None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    source_df = classified_df[classified_df["db_source"].astype(str) == db_source].copy()
    resolved_rows = []
    unresolved_rows = []

    for _, row in source_df.iterrows():
        tid = str(row["transcript_id"])
        resolved = pre_resolve_fn(row) if pre_resolve_fn else None

        if resolved is None:
            for candidate in candidates_fn(tid):
                if candidate in index:
                    resolved = build_resolved_row(
                        tid,
                        db_source=db_source,
                        organism=organism,
                        assembly_name=assembly_name,
                        hit=index[candidate],
                        assembly_accession=assembly_accession,
                        fasta_url=fasta_url,
                        gtf_url=gtf_url,
                        gtf_format=gtf_format,
                    )
                    break

        if resolved is not None:
            resolved_rows.append(resolved)
        else:
            unresolved_rows.append(
                build_unresolved_row(row, db_source=db_source, reason=unresolved_reason)
            )

    return (
        pd.DataFrame(resolved_rows, columns=RESOLVED_COLS),
        pd.DataFrame(unresolved_rows, columns=UNRESOLVED_COLS),
    )
