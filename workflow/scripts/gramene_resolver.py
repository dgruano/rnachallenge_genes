"""
scripts/gramene_resolver.py
Resolve legacy plant IDs using Gramene API
==========================================
Uses Gramene search endpoint to map legacy IDs to modern Ensembl IDs
"""

import sys
import time
from pathlib import Path
from typing import Dict, Optional, Tuple

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("gramene_resolver", snakemake.log[0])
unresolved_tsv = Path(snakemake.input[0])
resolved_output = Path(snakemake.output.resolved)
unresolved_output = Path(snakemake.output.unresolved)
cfg = snakemake.config

GRAMENE_API = "https://data.gramene.org/v69/search"
REQUEST_DELAY = 0.1  # seconds between requests
MAX_RETRIES = 3
GRAMENE_FIELDS = "id,synonyms,taxon_id,system_name,biotype,chrom,seq_region_name,start,end,strand,region"


def normalize_strand(value: object) -> str:
    if value in (1, "+", "+1", "1"):
        return "+"
    if value in (-1, "-", "-1"):
        return "-"
    return ""


def parse_coordinates(doc: Dict) -> Dict[str, object]:
    """Extract coordinate fields from a Gramene Solr document."""
    chrom = doc.get("chrom") or doc.get("seq_region_name") or ""
    start = doc.get("start")
    end = doc.get("end")
    strand = normalize_strand(doc.get("strand"))

    # Some records expose coordinates as region="chr:start-end".
    region = doc.get("region")
    if isinstance(region, str) and region:
        region_text = region.split(":")
        if len(region_text) >= 2:
            if not chrom:
                chrom = region_text[0]
            if start is None or end is None:
                span = region_text[1].split("-")
                if len(span) == 2:
                    try:
                        start = int(span[0])
                        end = int(span[1])
                    except ValueError:
                        pass
            if not strand and len(region_text) >= 3:
                strand = normalize_strand(region_text[2])

    return {
        "chrom": chrom,
        "start": start,
        "end": end,
        "strand": strand,
    }


def build_query_candidates(transcript_id: str) -> list[str]:
    """
    Build ordered Gramene query candidates from transcript ID.

    Examples:
        GRMZM5G842623_T01 → GRMZM5G842623
        Os03t0779300_01 → Os03t0779300
        Solyc01g005000.2.1 → Solyc01g005000.2, Solyc01g005000
    """
    candidates: list[str] = [transcript_id]

    gene_id = transcript_id
    if "_T" in gene_id:
        gene_id = gene_id.split("_T")[0]
        candidates.append(gene_id)
    elif "_" in gene_id:
        head, tail = gene_id.rsplit("_", 1)
        if tail.isdigit():
            gene_id = head
            candidates.append(gene_id)

    if "." in gene_id:
        parts = gene_id.split(".")
        if len(parts) > 2:
            candidates.append(".".join(parts[:-1]))
        candidates.append(parts[0])

    # de-duplicate while preserving order
    seen = set()
    uniq = []
    for candidate in candidates:
        if candidate not in seen:
            seen.add(candidate)
            uniq.append(candidate)
    return uniq


def infer_species_hint(row: Dict[str, object]) -> str:
    species_value = row.get("inferred_species")
    if pd.notna(species_value):
        return str(species_value)
    reason = str(row.get("reason", ""))
    if reason.startswith("biomart_no_match_"):
        return reason.replace("biomart_no_match_", "")
    return ""


def query_gramene(session: requests.Session, gene_id: str) -> Optional[Dict]:
    """
    Query Gramene search API for a gene ID.

    Returns:
        Dict with keys: modern_id, synonyms, taxon_id, system_name, region, start, end, strand
        None if not found or error
    """
    for attempt in range(MAX_RETRIES):
        try:
            # Restrict payload to required fields and first hit to reduce transfer/parse cost.
            params = {"q": gene_id, "fl": GRAMENE_FIELDS, "rows": 1}
            response = session.get(GRAMENE_API, params=params, timeout=30)

            if response.status_code == 200:
                data = response.json()
                num_found = data["response"]["numFound"]

                if num_found == 0:
                    return None

                # Take first result (highest score)
                doc = data["response"]["docs"][0]

                coords = parse_coordinates(doc)

                return {
                    "modern_id": doc["id"],
                    "synonyms": doc.get("synonyms", []),
                    "taxon_id": doc.get("taxon_id"),
                    "system_name": doc.get("system_name"),
                    "biotype": doc.get("biotype"),
                    "num_found": num_found,
                    "assembly_name": doc.get("system_name"),
                    **coords,
                }
            else:
                log.warning(
                    f"Gramene API returned HTTP {response.status_code} for {gene_id}"
                )

        except requests.RequestException as exc:
            log.warning(
                f"Gramene request failed (attempt {attempt+1}/{MAX_RETRIES}): {exc}"
            )
            if attempt < MAX_RETRIES - 1:
                time.sleep(1 * (attempt + 1))  # Exponential backoff
        except Exception as exc:
            log.error(f"Unexpected error querying {gene_id}: {exc}")
            return None

    return None


def resolve_via_gramene(
    unresolved_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Resolve transcript IDs via Gramene API.

    Returns:
        (resolved_df, still_unresolved_df)
    """
    resolved_records = []
    unresolved_records = []
    query_cache: Dict[str, Optional[Dict]] = {}
    api_calls = 0

    total = len(unresolved_df)
    log.info(f"Querying Gramene for {total} IDs...")

    session = requests.Session()

    for idx, row in enumerate(unresolved_df.to_dict("records"), start=1):
        transcript_id = row["transcript_id"]
        species_hint = infer_species_hint(row)

        candidates = build_query_candidates(transcript_id)
        result = None
        matched_candidate = ""
        for candidate in candidates:
            if candidate in query_cache:
                result = query_cache[candidate]
            else:
                result = query_gramene(session, candidate)
                query_cache[candidate] = result
                api_calls += 1
                time.sleep(REQUEST_DELAY)
            if result:
                matched_candidate = candidate
                break

        if result:
            # Map modern gene ID back to transcript coordinate
            modern_gene_id = result["modern_id"]

            synonyms = result["synonyms"]
            if isinstance(synonyms, list):
                synonyms_text = ";".join(synonyms)
            else:
                synonyms_text = str(synonyms or "")

            resolved_records.append(
                {
                    "transcript_id": transcript_id,
                    "gene_id": modern_gene_id,
                    "db_source": "gramene",
                    "organism": result["system_name"],
                    "assembly_name": result.get("assembly_name")
                    or result["system_name"],
                    "assembly_accession": pd.NA,
                    "chrom": result.get("chrom", ""),
                    "start": result.get("start", pd.NA),
                    "end": result.get("end", pd.NA),
                    "strand": result.get("strand", ""),
                    "is_ambiguous": False,
                    "species": result["system_name"],
                    "original_species_hint": species_hint,
                    "query_id": matched_candidate,
                    "gramene_synonyms": synonyms_text,
                    "biotype": result["biotype"],
                    "num_matches": result["num_found"],
                }
            )

            if idx % 50 == 0:
                log.info(f"Progress: {idx}/{total} ({len(resolved_records)} resolved)")
        else:
            unresolved_records.append(
                {
                    "transcript_id": transcript_id,
                    "normalized_gene_id": candidates[-1],
                    "inferred_species": species_hint,
                    "reason": "not_found_in_gramene",
                }
            )

    if total > 0:
        resolved_pct = 100 * len(resolved_records) / total
        log.info(
            f"Gramene resolution complete: {len(resolved_records)}/{total} "
            f"resolved ({resolved_pct:.1f}%), API calls={api_calls}, cache size={len(query_cache)}"
        )
    else:
        log.info("Gramene resolution complete: no IDs to resolve")

    resolved_df = pd.DataFrame(resolved_records)
    still_unresolved_df = pd.DataFrame(unresolved_records)

    return resolved_df, still_unresolved_df


# ── Main ───────────────────────────────────────────────────────

log.info(f"Reading unresolved IDs from {unresolved_tsv}")
unresolved_df = pd.read_csv(unresolved_tsv, sep="\t")
log.info(f"Loaded {len(unresolved_df)} unresolved IDs")

# Show breakdown by species
if "inferred_species" in unresolved_df.columns:
    species_counts = unresolved_df["inferred_species"].value_counts()
    log.info(f"Species breakdown:\n{species_counts}")

# Resolve via Gramene
resolved_df, still_unresolved_df = resolve_via_gramene(unresolved_df)

# Save results
resolved_output.parent.mkdir(parents=True, exist_ok=True)
unresolved_output.parent.mkdir(parents=True, exist_ok=True)

if len(resolved_df) > 0:
    resolved_df.to_csv(resolved_output, sep="\t", index=False)
    log.info(f"Saved {len(resolved_df)} resolved IDs to {resolved_output}")
else:
    # Create empty file with header
    pd.DataFrame(columns=["transcript_id", "gene_id", "db_source"]).to_csv(
        resolved_output, sep="\t", index=False
    )
    log.warning("No IDs resolved via Gramene")

if len(still_unresolved_df) > 0:
    still_unresolved_df.to_csv(unresolved_output, sep="\t", index=False)
    log.info(
        f"Saved {len(still_unresolved_df)} still-unresolved IDs to {unresolved_output}"
    )
else:
    # All resolved!
    pd.DataFrame(columns=["transcript_id", "reason"]).to_csv(
        unresolved_output, sep="\t", index=False
    )
    log.info("All IDs resolved!")

# Summary statistics
if len(resolved_df) > 0:
    log.info("\n=== Gramene Resolution Summary ===")
    log.info(f"Total queried: {len(unresolved_df)}")
    log.info(
        f"Resolved: {len(resolved_df)} ({100*len(resolved_df)/len(unresolved_df):.1f}%)"
    )
    log.info(
        f"Still unresolved: {len(still_unresolved_df)} ({100*len(still_unresolved_df)/len(unresolved_df):.1f}%)"
    )

    if "species" in resolved_df.columns:
        species_resolved = resolved_df["species"].value_counts()
        log.info(f"\nResolved by species:\n{species_resolved}")
