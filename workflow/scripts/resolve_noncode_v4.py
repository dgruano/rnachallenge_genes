"""
scripts/resolve_noncode_v4.py
Resolve NONCODE IDs that were not found in NONCODEv5, using the
NONCODEv4 whole-species BED archive as a fallback.

Resolution strategy
-------------------
Input : noncode_unresolved.tsv — IDs not found in NONCODEv5 BED/FASTA.
        NONCODEv5_Transcript2Gene — still used for gene_id lookup; covers
        many IDs that are numbered the same across releases.

For each unresolved ID:
  1. Try the full versioned ID as-is in the v4 BED (e.g. NONCELT024266.1).
     NONCODEv4 carries no version suffixes, so this probe will always miss,
     but keeping it preserves the version-first convention from the user's
     request.
  2. Try the base ID (version stripped, e.g. NONCELT024266).
  3. Gene-ID: for transcript IDs (NON*T*) probe T2G with full ID first,
     then base ID.  For gene IDs (NON*G*) the ID itself is the gene_id.
  4. If still not found → emit as unresolved with reason
     not_found_in_noncode_v4.

V4 BED assembly names differ from v5 (dm3 vs dm6, danRer7 vs danRer10, etc.)
and are reflected in assembly_accession of resolved output rows.

Outputs
-------
noncode_v4_resolved.tsv   — RESOLVED_COLS schema
noncode_v4_unresolved.tsv — transcript_id, raw_header, source_file, reason
"""
import re
import sys
import zipfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger  # type: ignore[import-untyped]

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("resolve_noncode_v4", snakemake.log[0])  # type: ignore[name-defined]

unresolved_path: str = snakemake.input.noncode_unresolved   # type: ignore[name-defined]
transcript2gene_path: str = snakemake.input.transcript2gene  # type: ignore[name-defined]
v4_bed_zip: str = snakemake.params.v4_bed_zip               # type: ignore[name-defined]
out_resolved: str = snakemake.output.resolved               # type: ignore[name-defined]
out_unresolved: str = snakemake.output.unresolved           # type: ignore[name-defined]

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

# ── NONCODEv4 assembly map ────────────────────────────────────
# Maps 7-char NON prefix → (organism, v4 assembly, BED file name suffix)
# Only assemblies present in NONCODEv4_wholeSpecies_lncAndGene_bed.zip.
V4_META: dict[str, tuple[str, str]] = {
    "NONCELT": ("Caenorhabditis elegans",   "ce10"),
    "NONCELG": ("Caenorhabditis elegans",   "ce10"),
    "NONDMET": ("Drosophila melanogaster",  "dm3"),
    "NONDMEG": ("Drosophila melanogaster",  "dm3"),
    "NONDRET": ("Danio rerio",              "danRer7"),
    "NONDREG": ("Danio rerio",              "danRer7"),
    "NONGGAT": ("Gallus gallus",            "galGal3"),
    "NONGGAG": ("Gallus gallus",            "galGal3"),
    "NONBTAT": ("Bos taurus",               "bosTau6"),
    "NONBTAG": ("Bos taurus",               "bosTau6"),
    "NONHSAT": ("Homo sapiens",             "hg19"),
    "NONHSAG": ("Homo sapiens",             "hg19"),
    "NONMMUT": ("Mus musculus",             "mm9"),
    "NONMMUG": ("Mus musculus",             "mm9"),
}

_NON_PREFIX_RE = re.compile(r"^(NON[A-Z]{3}[TG])\d+\.\d+$")


def _noncode_prefix(tid: str) -> str | None:
    m = _NON_PREFIX_RE.match(tid)
    return m.group(1) if m else None


# ── Load NONCODEv5_Transcript2Gene ────────────────────────────
log.info(f"Loading NONCODEv5_Transcript2Gene from {transcript2gene_path}")
t2g: dict[str, str] = {}
with open(transcript2gene_path) as fh:
    for line in fh:
        parts = line.rstrip("\n").split()
        if len(parts) >= 2:
            t2g[parts[0]] = parts[1]
log.info(f"  Transcript2Gene entries: {len(t2g):,}")

# ── Load unresolved IDs ───────────────────────────────────────
log.info(f"Loading unresolved NONCODE IDs from {unresolved_path}")
df_unres = pd.read_csv(unresolved_path, sep="\t", low_memory=False)
log.info(f"  Unresolved to attempt: {len(df_unres)}")

if df_unres.empty:
    log.info("No unresolved IDs — writing empty outputs.")
    pd.DataFrame(columns=RESOLVED_COLS).to_csv(out_resolved, sep="\t", index=False)
    pd.DataFrame(columns=UNRESOLVED_COLS).to_csv(out_unresolved, sep="\t", index=False)
    sys.exit(0)

# ── Load NONCODEv4 BED files from zip ────────────────────────
# BED IDs carry NO version suffix.  We probe full ID first (always a miss),
# then base ID (version stripped).
BedCoords = tuple[str, int, int, str]
_bed_cache: dict[str, dict[str, BedCoords]] = {}


def _load_bed(asm: str) -> dict[str, BedCoords]:
    if asm in _bed_cache:
        return _bed_cache[asm]

    inner = f"NONCODEv4_wholeSpecies_lncAndGene_bed/NONCODEv4_{asm}.lncAndGene.bed"
    coords: dict[str, BedCoords] = {}
    try:
        with zipfile.ZipFile(v4_bed_zip) as zf:
            with zf.open(inner) as fh:
                for raw in fh:
                    parts = raw.decode().rstrip("\n").split("\t")
                    if len(parts) < 6:
                        continue
                    tid = parts[3]
                    coords[tid] = (
                        parts[0],           # chrom
                        int(parts[1]) + 1,  # start 0-based → 1-based
                        int(parts[2]),      # end
                        parts[5],           # strand
                    )
    except KeyError:
        log.warning(f"  v4 BED entry missing in zip: {inner}")
    except FileNotFoundError:
        log.error(f"  v4 BED zip not found: {v4_bed_zip}")

    _bed_cache[asm] = coords
    log.info(f"  Loaded v4 BED {asm}: {len(coords):,} entries")
    return coords


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
            {"transcript_id": tid, "raw_header": raw_header,
             "source_file": source_file, "reason": "invalid_noncode_format"}
        )
        continue

    v4_meta = V4_META.get(prefix)
    if v4_meta is None:
        # Species not covered by NONCODEv4 (rat, opossum, platypus, etc.)
        still_unresolved_rows.append(
            {"transcript_id": tid, "raw_header": raw_header,
             "source_file": source_file, "reason": "not_in_noncode_v4_species"}
        )
        continue

    organism, v4_asm = v4_meta
    bed_index = _load_bed(v4_asm)
    base_id = tid.rsplit(".", 1)[0] if "." in tid else tid

    # Gene-ID lookup: T2G (full versioned → base ID → fallback to base_id).
    if prefix[-1] == "T":
        gene_id = t2g.get(tid) or t2g.get(base_id) or ""
    else:
        gene_id = tid

    # Probe 1: full versioned ID  (always misses for v4, but checked first per convention)
    bed_entry = bed_index.get(tid)

    # Probe 2: base ID (version stripped)
    if bed_entry is None:
        bed_entry = bed_index.get(base_id)

    if bed_entry is None:
        still_unresolved_rows.append(
            {"transcript_id": tid, "raw_header": raw_header,
             "source_file": source_file, "reason": "not_found_in_noncode_v4"}
        )
        continue

    chrom, start, end, strand = bed_entry
    resolved_rows.append(
        {
            "transcript_id":      tid,
            "db_source":          "noncode_v4",
            "gene_id":            gene_id if gene_id else base_id,
            "gene_symbol":        gene_id if gene_id else base_id,
            "organism":           organism,
            "assembly_accession": v4_asm,
            "chrom":              chrom,
            "start":              start,
            "end":                end,
            "strand":             strand,
            "is_ambiguous":       False,
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
log.info(f"NONCODEv4 resolved   : {len(df_res)}")
log.info(f"NONCODEv4 unresolved : {len(df_still)}")
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
log.info("resolve_noncode_v4 complete.")
