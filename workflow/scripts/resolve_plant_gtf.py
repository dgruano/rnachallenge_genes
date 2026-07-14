"""
scripts/resolve_plant_gtf.py
GTF-based resolver for plant transcript IDs
============================================
Consumes all classified plant IDs directly from classified_ids.tsv
rather than waiting for a BioMart step to emit them as unresolved.

For each species that has a GTF configured in plant_gtf_sources.yaml,
the script:
  1. Parses the gzipped GTF into a transcript → coordinates index
  2. Tries to match each unresolved transcript ID (and normalised
     variants) against the index
  3. Reports which GTF file and release resolved each transcript

Output columns match the biomart_plant_batch resolved schema so
merge_resolved can concatenate them directly.
"""

import gzip
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger
from resolution_guard import check_match_rates

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("resolve_plant_gtf", snakemake.log[0])

input_tsv = snakemake.input.classified  # classified_ids.tsv
gtf_files = snakemake.input.gtf_files  # list[str], one per species
out_resolved = snakemake.output.resolved
out_unresolved = snakemake.output.unresolved

# Config: {species: {url, assembly_name, release}}
GTF_SOURCES: Dict[str, dict] = snakemake.config.get("plant_gtf_sources", {})

# Build species → gtf_path mapping from the ordered gtf_files list.
# The rule expand() produces files in the same order as GTF_SOURCES.keys().
SPECIES_ORDER = list(GTF_SOURCES.keys())
SPECIES_TO_GTF: Dict[str, str] = {
    sp: gtf_files[i] for i, sp in enumerate(SPECIES_ORDER) if i < len(gtf_files)
}

RESOLVED_COLS = [
    "transcript_id",
    "db_source",
    "gene_id",
    "gene_symbol",
    "organism",
    "assembly_accession",
    "assembly_name",
    "chrom",
    "start",
    "end",
    "strand",
    "ensembl_plants_release",
    "is_ambiguous",
]
UNRESOLVED_COLS = ["transcript_id", "raw_header", "source_file", "reason"]


# ── ID prefix → species (mirrors biomart_plant_batch) ────────
PREFIX_TO_SPECIES = {
    "AT": "arabidopsis_thaliana",
    "Os": "oryza_sativa",
    "LOC_Os": "oryza_sativa",
    "OS": "oryza_sativa",
    "Zm": "zea_mays",
    "GRMZM": "zea_mays",
    "AC": "zea_mays",
    "Solyc": "solanum_lycopersicum",
    "Glyma": "glycine_max",
    "PGSC": "solanum_tuberosum",
    "VIT_": "vitis_vinifera",
    "GSVIVT": "vitis_vinifera",
    "GTVIVG": "vitis_vinifera",
    "orange": "citrus_sinensis",
}


def infer_species(transcript_id: str) -> Optional[str]:
    for prefix, species in PREFIX_TO_SPECIES.items():
        if transcript_id.startswith(prefix):
            return species
    return None


def generate_candidates(tid: str, species: str) -> List[str]:
    """Normalised ID variants — same logic as biomart_plant_batch."""
    candidates = [tid]

    if species == "oryza_sativa":
        if "_" in tid:
            candidates.append(tid.split("_")[0])
        if tid.startswith("OS"):
            candidates.append("Os" + tid[2:])
        m = re.search(r"t0?(\d+)", tid)
        if m and tid.endswith("t0" + m.group(1)):
            candidates.append(tid.replace("t0", "g0"))

    if species == "zea_mays":
        if "_" in tid:
            candidates.append(tid.split("_")[0])
        if "." in tid:
            parts = tid.split(".")
            while len(parts) > 1:
                parts = parts[:-1]
                candidates.append(".".join(parts))

    if species in ("solanum_lycopersicum", "solanum_tuberosum"):
        if tid.count(".") >= 2:
            candidates.append(tid.rsplit(".", 1)[0])
            candidates.append(tid.split(".")[0])
        elif "." in tid:
            candidates.append(tid.split(".")[0])
        if species == "solanum_tuberosum" and "DMT" in tid:
            candidates.append(tid.replace("DMT", "DMG"))

    if species == "glycine_max" and "." in tid:
        parts = tid.split(".")
        while len(parts) > 2:
            parts = parts[:-1]
            candidates.append(".".join(parts))

    return list(dict.fromkeys(candidates))


# ── GTF parser ────────────────────────────────────────────────

_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def _parse_attrs(attr_str: str) -> Dict[str, str]:
    return dict(_ATTR_RE.findall(attr_str))


def normalize_strand(value: str) -> str:
    if value in ("+", "1", "+1"):
        return "+"
    if value in ("-", "-1"):
        return "-"
    return "."


def build_gtf_index(gtf_path: str) -> Dict[str, dict]:
    """
    Parse a gzipped GTF and return a dict keyed by transcript_id.

    Only "transcript" feature lines are indexed; each entry holds:
        gene_id, gene_symbol, chrom, start, end, strand
    """
    index: Dict[str, dict] = {}
    opener = gzip.open if gtf_path.endswith(".gz") else open
    n_lines = 0

    with opener(gtf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            feature = parts[2]
            if feature != "transcript":
                continue

            attrs = _parse_attrs(parts[8])
            tid = attrs.get("transcript_id", "")
            if not tid:
                continue

            index[tid] = {
                "gene_id": attrs.get("gene_id", ""),
                "gene_symbol": attrs.get("gene_name", attrs.get("Name", "")),
                "chrom": parts[0],
                "start": parts[3],  # 1-based GTF coords
                "end": parts[4],
                "strand": normalize_strand(parts[6]),
            }
            n_lines += 1

    log.info(f"  Indexed {n_lines} transcripts from {gtf_path}")
    return index


# ── Main ─────────────────────────────────────────────────────
log.info("Loading classified plant IDs from classified_ids.tsv")

classified_df = pd.read_csv(input_tsv, sep="\t")
df = classified_df[classified_df["db_source"].astype(str) == "plant"].copy()
log.info(f"Plant IDs to resolve via GTF: {len(df)}")

if df.empty or not GTF_SOURCES:
    log.info("Nothing to resolve — writing empty outputs")
    pd.DataFrame(columns=RESOLVED_COLS).to_csv(out_resolved, sep="\t", index=False)
    pd.DataFrame(columns=UNRESOLVED_COLS).to_csv(out_unresolved, sep="\t", index=False)
    sys.exit(0)

# Lazy-load GTF indices: only parse files for species that appear in the input
df["inferred_species"] = df["transcript_id"].apply(infer_species)
needed_species = set(df["inferred_species"].dropna().unique()) & set(SPECIES_TO_GTF)

log.info(f"Building GTF indices for: {sorted(needed_species)}")
gtf_indices: Dict[str, Dict[str, dict]] = {}
for sp in needed_species:
    log.info(f"  Parsing GTF for {sp}: {SPECIES_TO_GTF[sp]}")
    gtf_indices[sp] = build_gtf_index(SPECIES_TO_GTF[sp])

resolved_rows: List[dict] = []
unresolved_rows: List[dict] = []

for _, row in df.iterrows():
    tid = str(row["transcript_id"])
    species = row.get("inferred_species")

    # Pass through IDs with unknown species or no GTF configured
    if not species or species not in gtf_indices:
        unresolved_rows.append(
            {
                "transcript_id": tid,
                "raw_header": row.get("raw_header", ""),
                "source_file": row.get("source_file", ""),
                "reason": (
                    "no_gtf_configured" if not species else f"no_gtf_for_{species}"
                ),
            }
        )
        continue

    index = gtf_indices[species]
    src_meta = GTF_SOURCES[species]

    matched = False
    for candidate in generate_candidates(tid, species):
        if candidate in index:
            hit = index[candidate]
            resolved_rows.append(
                {
                    "transcript_id": tid,
                    "db_source": "plant_gtf",
                    "gene_id": hit["gene_id"],
                    "gene_symbol": hit["gene_symbol"] or hit["gene_id"],
                    "organism": species,
                    "assembly_accession": "EnsemblPlants",
                    "assembly_name": src_meta["assembly_name"],
                    "chrom": hit["chrom"],
                    "start": hit["start"],
                    "end": hit["end"],
                    "strand": hit["strand"],
                    "ensembl_plants_release": src_meta["release"],
                    "is_ambiguous": False,
                }
            )
            matched = True
            break

    if not matched:
        unresolved_rows.append(
            {
                "transcript_id": tid,
                "raw_header": row.get("raw_header", ""),
                "source_file": row.get("source_file", ""),
                "reason": f"gtf_no_match_{species}",
            }
        )

# ── Write outputs ─────────────────────────────────────────────
res_df = pd.DataFrame(resolved_rows, columns=RESOLVED_COLS)
unres_df = pd.DataFrame(unresolved_rows, columns=UNRESOLVED_COLS)

res_df.to_csv(out_resolved, sep="\t", index=False)
unres_df.to_csv(out_unresolved, sep="\t", index=False)

log.info("GTF resolution complete:")
log.info(f"  Resolved:         {len(res_df)}")
log.info(f"  Still unresolved: {len(unres_df)}")

if not res_df.empty:
    for sp, grp in res_df.groupby("organism"):
        src = GTF_SOURCES.get(sp, {})
        log.info(
            f"  {sp}: {len(grp)} transcripts "
            f"[{src.get('assembly_name', '?')} / release {src.get('release', '?')}]"
        )

# ── Presence verification (namespace→assembly correctness guard) ──
strict = bool(snakemake.config.get("plant_resolution_strict", True))
min_rate = float(snakemake.config.get("plant_resolution_min_match_rate", 0.02))
attempted_counts = (
    df[df["inferred_species"].isin(gtf_indices)]
    .groupby("inferred_species")
    .size()
    .to_dict()
)
matched_counts = res_df.groupby("organism").size().to_dict() if not res_df.empty else {}
failures = check_match_rates(
    matched_counts, attempted_counts, min_rate=min_rate, log=log
)
if failures and strict:
    log.error(
        "Strict presence-check failed (plant_gtf); exiting non-zero. "
        "Set plant_resolution_strict: false in config to override."
    )
    sys.exit(1)
