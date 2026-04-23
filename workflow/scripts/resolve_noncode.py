"""
scripts/resolve_noncode.py
Resolve NONCODE v5 Transcript / Gene IDs
=========================================
Handles IDs classified as db_source == "noncode" by parse_ids.py.

Resolution strategy
-------------------
1.  For transcript IDs (NON*T*):  map to gene_id via NONCODEv5_Transcript2Gene.
    For gene IDs (NON*G*):        the ID itself is the gene_id.
2.  Genomic coordinates (chrom, start, end, strand) are looked up in the
    species-specific BED file (NONCODEv5_{assembly}.lncAndGene.bed.gz).
3.  FASTA files (NONCODEv5_{species}.fa.gz) are used to cross-check that
    uncoordinated IDs genuinely exist in the NONCODE release.

Outputs
-------
noncode_resolved.tsv   — RESOLVED_COLS schema
noncode_unresolved.tsv — transcript_id, raw_header, source_file, reason
"""
import gzip
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger  # type: ignore[import-untyped]

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("resolve_noncode", snakemake.log[0])  # type: ignore[name-defined]

classified_path: str = snakemake.input.classified        # type: ignore[name-defined]
transcript2gene_path: str = snakemake.input.transcript2gene  # type: ignore[name-defined]
bed_dir: str = snakemake.params.bed_dir                  # type: ignore[name-defined]
fa_dir: str = snakemake.params.fa_dir                    # type: ignore[name-defined]
out_resolved: str = snakemake.output.resolved            # type: ignore[name-defined]
out_unresolved: str = snakemake.output.unresolved        # type: ignore[name-defined]

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
# Keys are the 7-char NONCODE prefix (NON + 3-letter species code + T/G).
# Values: (organism, assembly_accession, BED assembly suffix, FASTA species name)
SPECIES_META: dict[str, tuple[str, str, str, str]] = {
    "NONDMET": ("Drosophila melanogaster",   "dm6",      "dm6",      "fruitfly"),
    "NONDMEG": ("Drosophila melanogaster",   "dm6",      "dm6",      "fruitfly"),
    "NONCELT": ("Caenorhabditis elegans",    "ce10",     "ce10",     "celegans"),
    "NONCELG": ("Caenorhabditis elegans",    "ce10",     "ce10",     "celegans"),
    "NONDRET": ("Danio rerio",               "danRer10", "danRer10", "zebrafish"),
    "NONDREG": ("Danio rerio",               "danRer10", "danRer10", "zebrafish"),
    "NONGGAT": ("Gallus gallus",             "galGal4",  "galGal4",  "chicken"),
    "NONGGAG": ("Gallus gallus",             "galGal4",  "galGal4",  "chicken"),
    "NONMDOT": ("Monodelphis domestica",     "monDom5",  "monDom5",  "opossum"),
    "NONMDOG": ("Monodelphis domestica",     "monDom5",  "monDom5",  "opossum"),
    "NONOANT": ("Ornithorhynchus anatinus",  "ornAna1",  "ornAna1",  "platypus"),
    "NONOANG": ("Ornithorhynchus anatinus",  "ornAna1",  "ornAna1",  "platypus"),
    "NONPPYT": ("Pongo abelii",              "ponAbe2",  "ponAbe2",  "orangutan"),
    "NONPPYG": ("Pongo abelii",              "ponAbe2",  "ponAbe2",  "orangutan"),
    "NONRATT": ("Rattus norvegicus",         "rn6",      "rn6",      "rat"),
    "NONRATG": ("Rattus norvegicus",         "rn6",      "rn6",      "rat"),
    "NONATHT": ("Arabidopsis thaliana",      "tair10",   "tair10",   "arabidopsis"),
    "NONATHG": ("Arabidopsis thaliana",      "tair10",   "tair10",   "arabidopsis"),
    "NONBTAT": ("Bos taurus",               "bosTau6",  "bosTau6",  "cow"),
    "NONBTAG": ("Bos taurus",               "bosTau6",  "bosTau6",  "cow"),
    # Human / Mouse / Chimp / Gorilla / Rhesus / Pig (if they appear)
    "NONHSAT": ("Homo sapiens",             "hg38",     "hg38",     "human"),
    "NONHSAG": ("Homo sapiens",             "hg38",     "hg38",     "human"),
    "NONMMUT": ("Mus musculus",             "mm10",     "mm10",     "mouse"),
    "NONMMUG": ("Mus musculus",             "mm10",     "mm10",     "mouse"),
    "NONPTRT": ("Pan troglodytes",          "panTro4",  "panTro4",  "chimp"),
    "NONPTRG": ("Pan troglodytes",          "panTro4",  "panTro4",  "chimp"),
    "NONGGOT": ("Gorilla gorilla",          "gorGor3",  "gorGor3",  "gorilla"),
    "NONGROG": ("Gorilla gorilla",          "gorGor3",  "gorGor3",  "gorilla"),
    "NONMMLT": ("Macaca mulatta",           "rheMac3",  "rheMac3",  "rhesus"),
    "NONMMLG": ("Macaca mulatta",           "rheMac3",  "rheMac3",  "rhesus"),
    "NONSSCG": ("Sus scrofa",               "susScr3",  "susScr3",  "pig"),
    "NONSSCT": ("Sus scrofa",               "susScr3",  "susScr3",  "pig"),
}

_NON_PREFIX_RE = re.compile(r"^(NON[A-Z]{3}[TG])\d+\.\d+$")


def _noncode_prefix(transcript_id: str) -> str | None:
    """Return the 7-char NON prefix (e.g. 'NONDMET', 'NONATHG') or None."""
    m = _NON_PREFIX_RE.match(transcript_id)
    return m.group(1) if m else None


# ── Load classified IDs filtered to noncode ───────────────────
log.info(f"Loading classified IDs from {classified_path}")
df_cls = pd.read_csv(classified_path, sep="\t", low_memory=False)
df_nc = df_cls[df_cls["db_source"] == "noncode"].copy()
log.info(f"  NONCODE IDs to resolve: {len(df_nc)}")

if df_nc.empty:
    log.warning("No NONCODE IDs found — writing empty outputs.")
    pd.DataFrame(columns=RESOLVED_COLS).to_csv(out_resolved, sep="\t", index=False)
    pd.DataFrame(columns=UNRESOLVED_COLS).to_csv(out_unresolved, sep="\t", index=False)
    sys.exit(0)

# ── Load Transcript2Gene ──────────────────────────────────────
log.info(f"Loading NONCODEv5_Transcript2Gene from {transcript2gene_path}")
t2g: dict[str, str] = {}
with open(transcript2gene_path) as fh:
    for line in fh:
        parts = line.rstrip("\n").split()
        if len(parts) >= 2:
            t2g[parts[0]] = parts[1]
log.info(f"  Transcript2Gene entries: {len(t2g):,}")

# ── Build BED index per assembly ──────────────────────────────
# Maps: bed_assembly → dict[transcript_id → (chrom, start_1based, end, strand)]
BedCoords = tuple[str, int, int, str]
bed_cache: dict[str, dict[str, BedCoords]] = {}
bed_base_cache: dict[str, dict[str, BedCoords]] = {}


def _strip_version(tid: str) -> str:
    """Strip trailing .N version suffix, e.g. NONDMET034404.2 → NONDMET034404"""
    return tid.rsplit(".", 1)[0] if "." in tid else tid


def _load_bed(bed_asm: str) -> dict[str, BedCoords]:
    """Load BED12 for a given assembly; cache and return coordinate dict."""
    if bed_asm in bed_cache:
        return bed_cache[bed_asm]

    bed_path = Path(bed_dir) / f"NONCODEv5_{bed_asm}.lncAndGene.bed.gz"
    if not bed_path.exists():
        log.warning(f"  BED file not found: {bed_path}")
        bed_cache[bed_asm] = {}
        bed_base_cache[bed_asm] = {}
        return {}

    coords: dict[str, BedCoords] = {}
    opener = gzip.open if bed_path.suffix == ".gz" else open
    with opener(bed_path, "rt") as fh:  # type: ignore[call-overload]
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom = parts[0]
            start_0 = int(parts[1])
            end = int(parts[2])
            tid = parts[3]
            strand = parts[5]
            coords[tid] = (chrom, start_0 + 1, end, strand)  # convert to 1-based start

    bed_cache[bed_asm] = coords
    bed_base_cache[bed_asm] = {_strip_version(k): v for k, v in coords.items()}
    log.info(f"  Loaded BED {bed_asm}: {len(coords):,} entries")
    return coords


# ── Build FASTA valid-ID sets per species name (optional check) ─
fa_cache: dict[str, set[str]] = {}


def _load_fa_ids(fa_name: str) -> set[str]:
    """Return the set of transcript IDs present in a NONCODE FASTA file."""
    if fa_name in fa_cache:
        return fa_cache[fa_name]

    fa_path = Path(fa_dir) / f"NONCODEv5_{fa_name}.fa.gz"
    if not fa_path.exists():
        log.warning(f"  FASTA file not found: {fa_path}")
        fa_cache[fa_name] = set()
        return set()

    ids: set[str] = set()
    opener = gzip.open if fa_path.suffix == ".gz" else open
    with opener(fa_path, "rt") as fh:  # type: ignore[call-overload]
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].rstrip())

    fa_cache[fa_name] = ids
    log.info(f"  Loaded FASTA {fa_name}: {len(ids):,} IDs")
    return ids


# ── Resolve each NONCODE ID ───────────────────────────────────
resolved_rows: list[dict] = []
unresolved_rows: list[dict] = []

for _, row in df_nc.iterrows():
    tid: str = str(row["transcript_id"])
    raw_header: str = str(row.get("raw_header", ""))
    source_file: str = str(row.get("source_file", ""))

    prefix = _noncode_prefix(tid)
    if prefix is None:
        unresolved_rows.append(
            {"transcript_id": tid, "raw_header": raw_header,
             "source_file": source_file, "reason": "invalid_noncode_format"}
        )
        continue

    meta = SPECIES_META.get(prefix)
    if meta is None:
        unresolved_rows.append(
            {"transcript_id": tid, "raw_header": raw_header,
             "source_file": source_file, "reason": f"unknown_noncode_prefix:{prefix}"}
        )
        continue

    organism, assembly, bed_asm, fa_name = meta

    # Gene-ID: for transcripts try T2G (full ID first, then version-stripped base).
    # For gene IDs the ID itself is the gene_id.
    base_id = tid.rsplit(".", 1)[0] if "." in tid else tid
    if prefix[-1] == "T":
        gene_id = t2g.get(tid) or t2g.get(base_id) or ""
    else:
        gene_id = tid

    # Fetch coordinates from NONCODEv5 BED (full versioned ID, then version-stripped fallback).
    bed_index = _load_bed(bed_asm)
    bed_entry = bed_index.get(tid)
    if bed_entry is None:
        bed_entry = bed_base_cache.get(bed_asm, {}).get(base_id)

    if bed_entry is not None:
        chrom, start, end, strand = bed_entry
    else:
        # Secondary: ID present in NONCODEv5 FASTA but no BED entry → NA coords.
        fa_ids = _load_fa_ids(fa_name)
        if tid in fa_ids:
            chrom, start, end, strand = "", pd.NA, pd.NA, ""
            log.debug(f"  {tid}: in v5 FASTA but no BED coords — coordinates NA")
        else:
            unresolved_rows.append(
                {"transcript_id": tid, "raw_header": raw_header,
                 "source_file": source_file, "reason": "not_found_in_noncode_v5"}
            )
            continue

    resolved_rows.append(
        {
            "transcript_id":      tid,
            "db_source":          "noncode",
            "gene_id":            gene_id if gene_id else base_id,
            "gene_symbol":        gene_id if gene_id else base_id,
            "organism":           organism,
            "assembly_accession": assembly,
            "chrom":              chrom,
            "start":              start,
            "end":                end,
            "strand":             strand,
            "is_ambiguous":       False,
        }
    )

# ── Write outputs ─────────────────────────────────────────────
df_resolved = pd.DataFrame(resolved_rows, columns=RESOLVED_COLS) if resolved_rows else pd.DataFrame(columns=RESOLVED_COLS)
df_unresolved = pd.DataFrame(unresolved_rows, columns=UNRESOLVED_COLS) if unresolved_rows else pd.DataFrame(columns=UNRESOLVED_COLS)

df_resolved.to_csv(out_resolved, sep="\t", index=False)
df_unresolved.to_csv(out_unresolved, sep="\t", index=False)

log.info("=" * 60)
log.info(f"NONCODE resolved     : {len(df_resolved)}")
log.info(f"NONCODE unresolved   : {len(df_unresolved)}")
if not df_resolved.empty:
    by_org = df_resolved.groupby("organism").size().sort_values(ascending=False)
    for org, count in by_org.items():
        log.info(f"  {org:<35}: {count}")
if not df_unresolved.empty:
    by_reason = df_unresolved.groupby("reason").size().sort_values(ascending=False)
    for reason, count in by_reason.items():
        log.info(f"  [unresolved] {reason}: {count}")
log.info(f"Written {out_resolved}")
log.info(f"Written {out_unresolved}")
log.info("resolve_noncode complete.")
