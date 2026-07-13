"""
scripts/resolve_noncode_2016.py
Third-tier fallback for NONCODE IDs not found in NONCODEv5 BED/FASTA or
NONCODEv4 BED.

Resolution strategy
-------------------
Input : noncode_v4_unresolved.tsv — IDs exhausted by v5 and v4 lookups.
        NONCODEv5_Transcript2Gene — gene_id lookup (shared numbering).
        NONCODE2016.fa            — single combined FASTA; existence proof only.

For each unresolved ID:
  1. Build an ID set from NONCODE2016.fa headers (loaded once into memory).
  2. Try the full versioned ID (e.g. NONCELT024012.1).
  3. Try the base ID (version stripped, e.g. NONCELT024012).
  4. If found in either form → resolved with NA genomic coordinates.
     Gene-ID is pulled from NONCODEv5_Transcript2Gene (full then base probe).
  5. If not found → emitted as unresolved, reason = not_found_in_any_noncode.

Outputs
-------
noncode_2016_resolved.tsv   — RESOLVED_COLS schema; chrom/start/end/strand NA
noncode_2016_unresolved.tsv — transcript_id, raw_header, source_file, reason
"""

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger  # type: ignore[import-untyped]

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("resolve_noncode_2016", snakemake.log[0])  # type: ignore[name-defined]

unresolved_path: str = snakemake.input.noncode_v4_unresolved  # type: ignore[name-defined]
nc2016_fa: str = snakemake.params.nc2016_fa  # type: ignore[name-defined]
out_resolved: str = snakemake.output.resolved  # type: ignore[name-defined]
out_unresolved: str = snakemake.output.unresolved  # type: ignore[name-defined]

# ── Schema ────────────────────────────────────────────────────
RESOLVED_COLS = [
    "transcript_id",
    "db_source",
    "gene_id",
    "gene_symbol",
    "organism",
    "assembly_accession",
    "chrom",
    "start",
    "end",
    "strand",
    "is_ambiguous",
]
UNRESOLVED_COLS = ["transcript_id", "raw_header", "source_file", "reason"]

# ── Species metadata ─────────────────────────────────────────
# Keys are the 7-char NONCODE prefix; values: (organism, canonical assembly).
# Assembly is kept for downstream reference even though coords will be NA.
SPECIES_META: dict[str, tuple[str, str]] = {
    "NONDMET": ("Drosophila melanogaster", "dm6"),
    "NONDMEG": ("Drosophila melanogaster", "dm6"),
    "NONCELT": ("Caenorhabditis elegans", "ce10"),
    "NONCELG": ("Caenorhabditis elegans", "ce10"),
    "NONDRET": ("Danio rerio", "danRer10"),
    "NONDREG": ("Danio rerio", "danRer10"),
    "NONGGAT": ("Gallus gallus", "galGal4"),
    "NONGGAG": ("Gallus gallus", "galGal4"),
    "NONMDOT": ("Monodelphis domestica", "monDom5"),
    "NONMDOG": ("Monodelphis domestica", "monDom5"),
    "NONOANT": ("Ornithorhynchus anatinus", "ornAna1"),
    "NONOANG": ("Ornithorhynchus anatinus", "ornAna1"),
    "NONPPYT": ("Pongo abelii", "ponAbe2"),
    "NONPPYG": ("Pongo abelii", "ponAbe2"),
    "NONRATT": ("Rattus norvegicus", "rn6"),
    "NONRATG": ("Rattus norvegicus", "rn6"),
    "NONATHT": ("Arabidopsis thaliana", "tair10"),
    "NONATHG": ("Arabidopsis thaliana", "tair10"),
    "NONBTAT": ("Bos taurus", "bosTau6"),
    "NONBTAG": ("Bos taurus", "bosTau6"),
    "NONHSAT": ("Homo sapiens", "hg38"),
    "NONHSAG": ("Homo sapiens", "hg38"),
    "NONMMUT": ("Mus musculus", "mm10"),
    "NONMMUG": ("Mus musculus", "mm10"),
    "NONPTRT": ("Pan troglodytes", "panTro4"),
    "NONPTRG": ("Pan troglodytes", "panTro4"),
    "NONGGOT": ("Gorilla gorilla", "gorGor3"),
    "NONGROG": ("Gorilla gorilla", "gorGor3"),
    "NONMMLT": ("Macaca mulatta", "rheMac3"),
    "NONMMLG": ("Macaca mulatta", "rheMac3"),
    "NONSSCG": ("Sus scrofa", "susScr3"),
    "NONSSCT": ("Sus scrofa", "susScr3"),
}

_NON_PREFIX_RE = re.compile(r"^(NON[A-Z]{3}[TG])\d+\.\d+$")


def _noncode_prefix(tid: str) -> str | None:
    m = _NON_PREFIX_RE.match(tid)
    return m.group(1) if m else None


# ── Load unresolved IDs ───────────────────────────────────────
log.info(f"Loading unresolved NONCODE IDs from {unresolved_path}")
df_unres = pd.read_csv(unresolved_path, sep="\t", low_memory=False)
log.info(f"  Unresolved to attempt: {len(df_unres)}")

if df_unres.empty:
    log.info("No unresolved IDs — writing empty outputs.")
    pd.DataFrame(columns=RESOLVED_COLS).to_csv(out_resolved, sep="\t", index=False)
    pd.DataFrame(columns=UNRESOLVED_COLS).to_csv(out_unresolved, sep="\t", index=False)
    sys.exit(0)

# ── Load NONCODE2016 FASTA headers into memory ────────────────
log.info(f"Loading NONCODE2016 FASTA IDs from {nc2016_fa}")
nc2016_ids: set[str] = set()
nc2016_base_ids: set[str] = set()
fa_path = Path(nc2016_fa)
if not fa_path.exists():
    log.error(f"  NONCODE2016 FASTA not found: {nc2016_fa}")
    sys.exit(1)

with open(nc2016_fa) as fh:
    for line in fh:
        if line.startswith(">"):
            tid = line[1:].rstrip()
            nc2016_ids.add(tid)
            base = tid.rsplit(".", 1)[0] if "." in tid else tid
            nc2016_base_ids.add(base)

log.info(f"  NONCODE2016 IDs loaded: {len(nc2016_ids):,}")

# ── Resolve ───────────────────────────────────────────────────
resolved_rows: list[dict] = []
still_unresolved_rows: list[dict] = []

for _, row in df_unres.iterrows():
    tid: str = str(row["transcript_id"])
    raw_header: str = str(row.get("raw_header", ""))
    source_file: str = str(row.get("source_file", ""))

    prefix = _noncode_prefix(tid)
    if prefix is None:
        still_unresolved_rows.append(
            {
                "transcript_id": tid,
                "raw_header": raw_header,
                "source_file": source_file,
                "reason": "invalid_noncode_format",
            }
        )
        continue

    meta = SPECIES_META.get(prefix)
    if meta is None:
        still_unresolved_rows.append(
            {
                "transcript_id": tid,
                "raw_header": raw_header,
                "source_file": source_file,
                "reason": f"unknown_noncode_prefix:{prefix}",
            }
        )
        continue

    base_id = tid.rsplit(".", 1)[0] if "." in tid else tid

    # Probe NONCODE2016 FASTA: full ID first, then base ID.
    in_2016 = tid in nc2016_ids or base_id in nc2016_base_ids
    if not in_2016:
        still_unresolved_rows.append(
            {
                "transcript_id": tid,
                "raw_header": raw_header,
                "source_file": source_file,
                "reason": "not_found_in_any_noncode",
            }
        )
        continue

    # NONCODE2016 fallback confirms transcript existence only; coordinates are unavailable.
    still_unresolved_rows.append(
        {
            "transcript_id": tid,
            "raw_header": raw_header,
            "source_file": source_file,
            "reason": "matched_noncode2016_no_coordinates",
        }
    )

# ── Write outputs ─────────────────────────────────────────────
df_res = (
    pd.DataFrame(resolved_rows, columns=RESOLVED_COLS)
    if resolved_rows
    else pd.DataFrame(columns=RESOLVED_COLS)
)
df_still = (
    pd.DataFrame(still_unresolved_rows, columns=UNRESOLVED_COLS)
    if still_unresolved_rows
    else pd.DataFrame(columns=UNRESOLVED_COLS)
)

df_res.to_csv(out_resolved, sep="\t", index=False)
df_still.to_csv(out_unresolved, sep="\t", index=False)

log.info("=" * 60)
log.info(f"NONCODE2016 resolved   : {len(df_res)}")
log.info(f"NONCODE2016 unresolved : {len(df_still)}")
if not df_res.empty:
    by_org = df_res.groupby("organism").size().sort_values(ascending=False)
    for org, count in by_org.items():
        log.info(f"  {org:<35}: {count}")
if not df_still.empty:
    by_reason = df_still.groupby("reason").size().sort_values(ascending=False)
    for reason, count in by_reason.items():
        log.info(f"  [unresolved] {reason}: {count}")
log.info(f"Written {out_resolved}")
log.info(f"Written {out_unresolved}")
log.info("resolve_noncode_2016 complete.")
