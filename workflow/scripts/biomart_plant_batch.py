"""
scripts/biomart_plant_batch.py
Batch query Ensembl Plants BioMart for plant transcript IDs
===========================================================
Unlike REST API, BioMart can handle batch queries efficiently.
Queries multiple IDs at once grouped by species.

If a species query returns no results on the current Ensembl Plants release,
the script falls back through a ranked list of older releases and records
which release ultimately provided the data.
"""

import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("biomart_plant_batch", snakemake.log[0])
input_tsv    = snakemake.input.unresolved
out_resolved   = snakemake.output.resolved
out_unresolved = snakemake.output.unresolved
cfg = snakemake.config

MAX_RETRIES = int(cfg.get("max_retries", 3))
RETRY_WAIT  = int(cfg.get("retry_wait_seconds", 5))
BATCH_SIZE  = 500

# Ensembl Plants release fallback chain (newest → oldest).
# Format: (release_label, biomart_url)
# "current" uses the stable hostname; numbered releases use versioned hostnames.
ENSEMBL_PLANTS_RELEASES: List[Tuple[str, str]] = [
    ("current", "https://plants.ensembl.org/biomart/martservice"),
    ("release-60", "https://release-60-plants.ensembl.org/biomart/martservice"),
    ("release-59", "https://release-59-plants.ensembl.org/biomart/martservice"),
    ("release-58", "https://release-58-plants.ensembl.org/biomart/martservice"),
    ("release-57", "https://release-57-plants.ensembl.org/biomart/martservice"),
]

# Species to BioMart dataset mapping
SPECIES_TO_DATASET = {
    "arabidopsis_thaliana": "athaliana_eg_gene",
    "oryza_sativa":         "osativa_eg_gene",
    "zea_mays":             "zmays_eg_gene",
    "solanum_lycopersicum": "slycopersicum_eg_gene",
    "glycine_max":          "gmax_eg_gene",
    "solanum_tuberosum":    "stuberosum_eg_gene",
    "citrus_sinensis":      "csinensis_eg_gene",
}

# ID prefix to species mapping
PREFIX_TO_SPECIES = {
    "AT":     "arabidopsis_thaliana",
    "Os":     "oryza_sativa",
    "LOC_Os": "oryza_sativa",
    "OS":     "oryza_sativa",
    "Zm":     "zea_mays",
    "GRMZM":  "zea_mays",
    "AC":     "zea_mays",
    "Solyc":  "solanum_lycopersicum",
    "Glyma":  "glycine_max",
    "PGSC":   "solanum_tuberosum",
    "orange": "citrus_sinensis",
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


# ── Helpers ───────────────────────────────────────────────────

def normalize_strand(value) -> str:
    if value in (1, "1", "+", "+1"):
        return "+"
    if value in (-1, "-1", "-", "-1"):
        return "-"
    return "."


def get_species_from_id(transcript_id: str) -> Optional[str]:
    """Infer species from transcript ID prefix."""
    for prefix, species in PREFIX_TO_SPECIES.items():
        if transcript_id.startswith(prefix):
            return species
    return None


def generate_id_candidates(tid: str, species: str) -> List[str]:
    """Generate normalised ID variants for BioMart queries."""
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

    if species == "solanum_lycopersicum":
        if tid.count(".") >= 2:
            candidates.append(tid.rsplit(".", 1)[0])
            candidates.append(tid.split(".")[0])
        elif "." in tid:
            candidates.append(tid.split(".")[0])

    if species == "solanum_tuberosum":
        if "DMT" in tid:
            candidates.append(tid.replace("DMT", "DMG"))
        if "." in tid:
            candidates.append(tid.split(".")[0])

    if species == "glycine_max":
        if "." in tid:
            parts = tid.split(".")
            while len(parts) > 2:
                parts = parts[:-1]
                candidates.append(".".join(parts))

    return list(dict.fromkeys(candidates))  # deduplicate, preserve order


def build_biomart_query(dataset: str, transcript_ids: List[str]) -> str:
    """Build BioMart XML query for batch transcript lookup."""
    id_filters = ",".join(transcript_ids)
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        "<!DOCTYPE Query>\n"
        '<Query virtualSchemaName="plants_mart" formatter="TSV" '
        'header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">\n'
        f'    <Dataset name="{dataset}" interface="default">\n'
        f'        <Filter name="ensembl_transcript_id" value="{id_filters}"/>\n'
        '        <Attribute name="ensembl_transcript_id" />\n'
        '        <Attribute name="ensembl_gene_id" />\n'
        '        <Attribute name="external_gene_name" />\n'
        '        <Attribute name="assembly_name" />\n'
        '        <Attribute name="chromosome_name" />\n'
        '        <Attribute name="transcript_start" />\n'
        '        <Attribute name="transcript_end" />\n'
        '        <Attribute name="strand" />\n'
        "    </Dataset>\n"
        "</Query>"
    )


def _check_endpoint(url: str) -> bool:
    """Lightweight liveness check — avoids spending retry budget on dead endpoints."""
    try:
        resp = requests.get(url, params={"type": "registry"}, timeout=15)
        return resp.status_code == 200
    except (requests.exceptions.Timeout, requests.exceptions.ConnectionError):
        return False


def _post_query(url: str, query: str) -> Optional[pd.DataFrame]:
    """
    POST a BioMart XML query.

    Returns:
      None            — transient failure (network/HTTP error, BioMart ERROR body);
                        caller should retry.
      empty DataFrame — server is up but returned no results;
                        caller should NOT retry this URL, but may try next release.
      populated DataFrame — success.
    """
    try:
        resp = requests.post(url, data={"query": query}, timeout=120)
    except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as exc:
        log.warning(f"Request error: {exc}")
        return None

    if resp.status_code != 200:
        log.warning(f"HTTP {resp.status_code} from {url}")
        return None

    # BioMart returns HTTP 200 with "ERROR ..." body on malformed queries.
    text = resp.text.strip()
    if text.upper().startswith("ERROR"):
        log.warning(f"BioMart query error from {url}: {text[:200]}")
        return None  # retriable — may be a transient mart error

    lines = text.split("\n")
    if len(lines) < 2:
        return pd.DataFrame()  # server up, no results — don't retry

    header = lines[0].split("\t")
    data   = [line.split("\t") for line in lines[1:] if line.strip()]
    return pd.DataFrame(data, columns=header)


def query_biomart_with_fallback(
    dataset: str,
    transcript_ids: List[str],
) -> Tuple[pd.DataFrame, str]:
    """
    Query BioMart for *transcript_ids* in *dataset*, trying each release in
    ENSEMBL_PLANTS_RELEASES in order.  Returns (DataFrame, release_label).
    DataFrame is empty if all releases fail or return no matches.

    Release fallback strategy:
    - Skip release entirely if endpoint is unreachable (liveness check).
    - Retry only on None (transient failure); an empty-but-valid response
      means the server is up — no point retrying, but older releases may
      still have retired IDs so the fallback chain continues.
    - Return immediately on the first non-empty result.
    """
    query = build_biomart_query(dataset, transcript_ids)

    for release_label, url in ENSEMBL_PLANTS_RELEASES:
        log.info(
            f"Trying Ensembl Plants {release_label} for {dataset} "
            f"({len(transcript_ids)} IDs)"
        )

        if not _check_endpoint(url):
            log.warning(f"  {release_label} endpoint unreachable — skipping")
            continue

        df: Optional[pd.DataFrame] = None
        for attempt in range(1, MAX_RETRIES + 1):
            df = _post_query(url, query)
            if df is not None:
                break  # server responded (may be empty) — stop retrying this URL
            if attempt < MAX_RETRIES:
                wait = RETRY_WAIT * attempt
                log.warning(
                    f"  Attempt {attempt}/{MAX_RETRIES} failed for "
                    f"{release_label} — retrying in {wait}s"
                )
                time.sleep(wait)

        if df is None:
            # Endpoint failed all retries — try next release
            log.warning(
                f"  {release_label} failed after {MAX_RETRIES} attempts — trying next release"
            )
            continue

        if not df.empty:
            log.info(f"  Success on {release_label}: {len(df)} rows returned")
            return df, release_label

        # df is empty: server is live, IDs not in this release — try older releases
        log.info(f"  {release_label} returned no matches — trying older release")

    log.error(
        f"All Ensembl Plants releases exhausted for {dataset} — "
        f"marking {len(transcript_ids)} IDs as unresolved"
    )
    return pd.DataFrame(), "none"


# ── Main ─────────────────────────────────────────────────────
log.info("Starting BioMart batch plant ID resolution")

df = pd.read_csv(input_tsv, sep="\t")
log.info(f"Loaded {len(df)} unresolved IDs")

df["inferred_species"] = df["transcript_id"].apply(get_species_from_id)
plant_df    = df[df["inferred_species"].notna()].copy()
non_plant_df = df[df["inferred_species"].isna()]
log.info(f"Found {len(plant_df)} plant IDs to query via BioMart")

resolved_rows:   List[dict] = []
unresolved_rows: List[dict] = []

for species, group in plant_df.groupby("inferred_species"):
    if species not in SPECIES_TO_DATASET:
        log.warning(f"No BioMart dataset configured for species: {species}")
        for _, row in group.iterrows():
            unresolved_rows.append({
                "transcript_id": row["transcript_id"],
                "raw_header":    row.get("raw_header", ""),
                "source_file":   row.get("source_file", ""),
                "reason":        f"no_biomart_dataset_{species}",
            })
        continue

    dataset = SPECIES_TO_DATASET[species]
    log.info(f"Processing {len(group)} IDs for {species} (dataset: {dataset})")

    # Build candidate ID list
    all_candidates: Dict[str, dict] = {}
    for _, row in group.iterrows():
        tid = str(row["transcript_id"])
        all_candidates[tid] = {
            "candidates": generate_id_candidates(tid, species),
            "raw_header": row.get("raw_header", ""),
            "source_file": row.get("source_file", ""),
        }

    unique_candidates = list(
        dict.fromkeys(c for info in all_candidates.values() for c in info["candidates"])
    )
    log.info(f"  {len(unique_candidates)} unique candidate IDs generated")

    # Query in batches, collecting all results
    all_results: List[pd.DataFrame] = []
    used_release = "none"
    for i in range(0, len(unique_candidates), BATCH_SIZE):
        batch = unique_candidates[i : i + BATCH_SIZE]
        batch_df, release_label = query_biomart_with_fallback(dataset, batch)
        if not batch_df.empty:
            batch_df["_release"] = release_label
            all_results.append(batch_df)
            used_release = release_label  # last successful release (may vary per batch)
        time.sleep(1)  # rate-limit courtesy pause

    if not all_results:
        log.warning(f"  No BioMart results for {species}")
        for tid, info in all_candidates.items():
            unresolved_rows.append({
                "transcript_id": tid,
                "raw_header":    info["raw_header"],
                "source_file":   info["source_file"],
                "reason":        f"biomart_no_match_{species}_all_releases_exhausted",
            })
        continue

    biomart_df = pd.concat(all_results, ignore_index=True)

    # Normalise column names (BioMart returns human-readable headers)
    col_map = {
        "Transcript stable ID":       "queried_id",
        "Gene stable ID":             "gene_id",
        "Gene name":                  "gene_symbol",
        "Assembly":                   "assembly_name",
        "Chromosome/scaffold name":   "chrom",
        "Transcript start (bp)":      "start",
        "Transcript end (bp)":        "end",
        "Strand":                     "strand",
    }
    biomart_df = biomart_df.rename(columns=col_map)

    # Build lookup: queried_id → result row
    biomart_lookup: Dict[str, dict] = {}
    for _, row in biomart_df.iterrows():
        qid = str(row.get("queried_id", ""))
        if qid and qid not in biomart_lookup:
            biomart_lookup[qid] = {
                "gene_id":       row.get("gene_id", ""),
                "gene_symbol":   row.get("gene_symbol", ""),
                "assembly_name": row.get("assembly_name", ""),
                "chrom":         row.get("chrom", ""),
                "start":         row.get("start", ""),
                "end":           row.get("end", ""),
                "strand":        normalize_strand(row.get("strand", "")),
                "release":       row.get("_release", used_release),
            }

    # Match original IDs → resolved / unresolved
    for tid, info in all_candidates.items():
        matched = False
        for candidate in info["candidates"]:
            if candidate in biomart_lookup:
                m = biomart_lookup[candidate]
                resolved_rows.append({
                    "transcript_id":        tid,
                    "db_source":            "plant",
                    "gene_id":              m["gene_id"],
                    "gene_symbol":          m["gene_symbol"] or m["gene_id"],
                    "organism":             species,
                    "assembly_accession":   "EnsemblPlants",
                    "assembly_name":        m["assembly_name"],
                    "chrom":                m["chrom"],
                    "start":                m["start"],
                    "end":                  m["end"],
                    "strand":               m["strand"],
                    "ensembl_plants_release": m["release"],
                    "is_ambiguous":         False,
                })
                matched = True
                break

        if not matched:
            unresolved_rows.append({
                "transcript_id": tid,
                "raw_header":    info["raw_header"],
                "source_file":   info["source_file"],
                "reason":        f"biomart_no_match_{species}",
            })

# Non-plant IDs pass through as unresolved
for _, row in non_plant_df.iterrows():
    unresolved_rows.append({
        "transcript_id": row["transcript_id"],
        "raw_header":    row.get("raw_header", ""),
        "source_file":   row.get("source_file", ""),
        "reason":        "not_plant_id",
    })

# ── Write outputs ─────────────────────────────────────────────
res_df   = pd.DataFrame(resolved_rows,   columns=RESOLVED_COLS)
unres_df = pd.DataFrame(unresolved_rows, columns=UNRESOLVED_COLS)

res_df.to_csv(out_resolved,   sep="\t", index=False)
unres_df.to_csv(out_unresolved, sep="\t", index=False)

log.info("BioMart batch resolution complete:")
log.info(f"  Resolved:         {len(res_df)}")
log.info(f"  Still unresolved: {len(unres_df)}")

if not res_df.empty:
    release_summary = res_df.groupby("ensembl_plants_release").size()
    log.info("  Resolved by release:")
    for rel, count in release_summary.items():
        log.info(f"    {rel}: {count}")
    assembly_summary = res_df.groupby(["organism", "assembly_name"]).size()
    log.info("  Assembly used per species:")
    for (org, asm), count in assembly_summary.items():
        log.info(f"    {org}: {asm} ({count} transcripts)")
