"""
scripts/resolve_external_ids.py
Resolve external transcript IDs (plant)
================================================================
Input: classified TSV with columns transcript_id, raw_header, source_file, db_source.
Outputs: external_resolved.tsv, external_ambiguous.tsv, external_unresolved.tsv.

Resolution strategy by db_source
---------------------------------

Plant  (db_source == "plant")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Species is inferred from transcript ID prefix (PLANT_PREFIX_TO_SPECIES dict).
Candidates are produced by plant_candidates(), which extends generic_candidates()
with species-specific normalizations:
  • GRMZM — gene ID is the prefix before the first underscore
  • OS/LOC_Os (rice) — case-insensitive prefix variants and LOC_Os → Os conversions
  • Solyc (tomato) — strips one or both trailing dot-delimited version components
  • Glyma (soybean) — converts "Glyma.10G011600.1" → "GLYMA_10G011600" (dot→underscore, uppercase)

Resolution steps (in order):
1. Local metadata table  — if species has an entry in config.external_metadata_tables,
   look up each candidate key in the preloaded dict; returns coordinates without any
   REST call.  Used for plants because plants.rest.ensembl.org is frequently unavailable.
2. Ensembl lookup/id    — skipped when SKIP_REST_FOR_PLANTS=true; tries each candidate
   via /lookup/id/{cand}.
3. Ensembl xrefs/symbol — skipped when SKIP_REST_FOR_PLANTS=true; tries each candidate
   via /xrefs/symbol/{species}/{cand}, then follows each xref through lookup/id.
Reason codes: plant_id_not_resolved

Global controls
---------------
SKIP_REST_LOOKUPS      — skip all REST calls (dry-run / offline mode).
SKIP_REST_FOR_PLANTS   — skip REST for plant species only; metadata tables still used.
MAX_RETRIES / RETRY_WAIT — retry budget per request (exponential back-off).
MAX_503_TOTAL / MAX_503_CONSECUTIVE — abort thresholds to avoid infinite 503 loops.
"""

import re
import sys
import time
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("resolve_external_ids", snakemake.log[0])
input_tsv = snakemake.input.classified
out_resolved = snakemake.output.resolved
out_ambig = snakemake.output.ambiguous
out_unresolved = snakemake.output.unresolved
cfg = snakemake.config

REST_BASE = cfg.get("ensembl_rest_base", "https://rest.ensembl.org")
SKIP_REST_LOOKUPS = cfg.get("skip_rest_lookups", False)
SKIP_REST_FOR_PLANTS = cfg.get("skip_rest_for_plants", False)
MAX_RETRIES = int(cfg.get("max_retries", 3))
RETRY_WAIT = int(cfg.get("retry_wait_seconds", 5))
MAX_503_TOTAL = int(cfg.get("max_503_total", 50))
MAX_503_CONSECUTIVE = int(cfg.get("max_503_consecutive", 15))

# Create a persistent session for REST API connection pooling
REST_SESSION = requests.Session()

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
AMBIG_COLS = [
    "transcript_id",
    "db_source",
    "chosen_gene_id",
    "alternative_gene_id",
    "alternative_gene_symbol",
    "organism",
    "assembly_accession",
    "chrom",
    "start",
    "end",
    "strand",
]
UNRESOLVED_COLS = ["transcript_id", "raw_header", "source_file", "reason"]

PLANT_PREFIX_TO_SPECIES = {
    "Solyc": ["solanum_lycopersicum"],
    "OS": ["oryza_sativa"],
    "Glyma": ["glycine_max"],
    "AT": ["arabidopsis_thaliana"],
    "Zm": ["zea_mays"],
    "GRMZM": ["zea_mays"],
    "LOC_Os": ["oryza_sativa"],
    "Bradi": ["brachypodium_distachyon"],
    "TraesCS": ["triticum_aestivum"],
    "PGSC": ["solanum_tuberosum"],
    "Potri": ["populus_trichocarpa"],
    "Sobic": ["sorghum_bicolor"],
    "VIT_": ["vitis_vinifera"],
    "Bra": ["brassica_rapa"],
    "BnaA": ["brassica_napus"],
    "BnaC": ["brassica_napus"],
    "Bo": ["brassica_oleracea"],
    "AET": ["aegilops_tauschii"],
    "Amtr": ["amborella_trichopoda"],
    "evm.model": ["amborella_trichopoda"],
    "Cre": ["chlamydomonas_reinhardtii"],
    "Pp": ["physcomitrella_patens"],
    "Medtr": ["medicago_truncatula"],
    "GSMUA": ["musa_acuminata"],
    "OB": ["oryza_brachyantha"],
    "Si": ["setaria_italica"],
    "Thecc1EG": ["theobroma_cacao"],
    "orange1.1": ["citrus_sinensis"],
    "cassava": ["manihot_esculenta"],
    "AC": ["zea_mays"],
}

REST_503_TOTAL = 0
REST_503_CONSEC = 0


# ── Helpers ───────────────────────────────────────────────────
def get_rest_base_for_species(species: str) -> str:
    # Use common REST_BASE for all species (plants subdomain frequently unavailable)
    return REST_BASE


def is_plant_species(species: str) -> bool:
    """Check if a species is a plant (should use metadata, not REST)."""
    plant_species = {
        "arabidopsis_thaliana",
        "oryza_sativa",
        "zea_mays",
        "solanum_lycopersicum",
        "glycine_max",
        "brachypodium_distachyon",
        "triticum_aestivum",
        "solanum_tuberosum",
        "populus_trichocarpa",
        "sorghum_bicolor",
        "vitis_vinifera",
        "brassica_rapa",
        "brassica_napus",
        "brassica_oleracea",
        "aegilops_tauschii",
        "amborella_trichopoda",
        "chlamydomonas_reinhardtii",
        "physcomitrella_patens",
        "medicago_truncatula",
        "musa_acuminata",
        "oryza_brachyantha",
        "setaria_italica",
        "theobroma_cacao",
        "citrus_sinensis",
        "manihot_esculenta",
    }
    return species in plant_species


def rest_get(
    path: str,
    params: Optional[dict] = None,
    base_url: Optional[str] = None,
    species: Optional[str] = None,
) -> Optional[dict]:
    global REST_503_TOTAL
    global REST_503_CONSEC

    if SKIP_REST_LOOKUPS:
        return None

    # Skip REST for plant species if configured (Ensembl Plants down)
    if SKIP_REST_FOR_PLANTS and species and is_plant_species(species):
        return None

    url = f"{(base_url or REST_BASE)}{path}"
    headers = {"Content-Type": "application/json"}
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = REST_SESSION.get(url, params=params, headers=headers, timeout=30)
            if resp.status_code == 200:
                REST_503_CONSEC = 0
                return resp.json()
            if resp.status_code in (400, 404):
                REST_503_CONSEC = 0
                return None
            if resp.status_code == 503:
                REST_503_TOTAL += 1
                REST_503_CONSEC += 1
                if "Domain not existing" in resp.text:
                    log.error(
                        f"  REST {path} attempt {attempt} failed: 503 Domain not existing from {resp.url}"
                    )
                if (
                    REST_503_TOTAL >= MAX_503_TOTAL
                    or REST_503_CONSEC >= MAX_503_CONSECUTIVE
                ):
                    raise RuntimeError(
                        "Too many 503 responses from Ensembl REST; aborting to avoid endless retries."
                    )
            summary = " ".join(resp.text.splitlines())[:200]
            log.warning(
                f"  REST {path} attempt {attempt} failed: {resp.status_code} {summary}"
            )
        except Exception as exc:
            log.warning(f"  REST {path} attempt {attempt} failed: {exc}")
        if attempt < MAX_RETRIES:
            time.sleep(RETRY_WAIT * attempt)
    return None


def load_metadata_table(path: str) -> dict:
    """Load a metadata table keyed by transcript_id for fast local resolution.

    Expected columns:
      transcript_id, gene_id, gene_symbol, assembly_accession, chrom, start, end, strand

    Returns a dict indexed by both full transcript IDs and gene ID variants to improve
    candidate matching when version suffixes differ.
    """
    try:
        table = pd.read_csv(path, sep="\t")
    except Exception as exc:
        log.warning(f"Failed to read metadata table {path}: {exc}")
        return {}

    required = {
        "transcript_id",
        "gene_id",
        "gene_symbol",
        "assembly_accession",
        "chrom",
        "start",
        "end",
        "strand",
    }
    missing = required - set(table.columns)
    if missing:
        log.warning(
            f"Metadata table {path} missing columns: {', '.join(sorted(missing))}"
        )
        return {}

    index = {}
    for _, row in table.iterrows():
        metadata = {
            "gene_id": row["gene_id"],
            "gene_symbol": row["gene_symbol"],
            "assembly_accession": row["assembly_accession"],
            "chrom": row["chrom"],
            "start": row["start"],
            "end": row["end"],
            "strand": row["strand"],
        }

        # Index by full transcript ID
        transcript_id = str(row["transcript_id"])
        index[transcript_id] = metadata

        # Also index by gene ID (base without transcript version)
        # This allows matching Solyc06g068790.2.1 → Solyc06g068790.3.1 via gene base
        gene_id = str(row["gene_id"])
        if gene_id and gene_id not in index:
            index[gene_id] = metadata

        # Index by transcript base (strip trailing version)
        # e.g., Solyc06g068790.3.1 → Solyc06g068790.3 and Solyc06g068790
        if "." in transcript_id:
            parts = transcript_id.split(".")
            # Strip last component (Solyc06g068790.3.1 → Solyc06g068790.3)
            base_variant = ".".join(parts[:-1])
            if base_variant not in index:
                index[base_variant] = metadata
            # Strip all version components (Solyc06g068790.3.1 → Solyc06g068790)
            base = parts[0]
            if base not in index:
                index[base] = metadata

    return index


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


def plant_candidates(tid: str) -> list[str]:
    candidates = generic_candidates(tid)

    # GRMZM (Maize): Gene ID is prefix before underscore
    if tid.startswith("GRMZM"):
        gene_id = tid.split("_")[0]
        candidates.append(gene_id)

    # Os (Rice): Gene ID is prefix before underscore
    if tid.startswith("Os") and "_" in tid:
        gene_id = tid.split("_")[0]
        candidates.append(gene_id)

    # Solyc (Tomato): Two dot suffixes - last is transcript, first is gene version
    # Strip last suffix first, then both if needed
    if tid.startswith("Solyc"):
        if tid.count(".") >= 2:
            # Strip last suffix only (e.g., Solyc01g005000.2.1 → Solyc01g005000.2)
            candidates.append(tid.rsplit(".", 1)[0])
            # Strip both suffixes (e.g., Solyc01g005000.2.1 → Solyc01g005000)
            candidates.append(tid.split(".")[0])
        elif tid.count(".") == 1:
            # Single dot - strip it
            candidates.append(tid.split(".")[0])

    # Rice OS prefix variants
    if tid.startswith("OS"):
        candidates.append("Os" + tid[2:])
        match = re.match(r"OS(\d+)T(\d+)", tid, re.IGNORECASE)
        if match:
            candidates.append(f"LOC_Os{match.group(1)}g{match.group(2)}")
            candidates.append(f"Os{match.group(1)}t{match.group(2)}-01")
            candidates.append(f"Os{match.group(1)}t{match.group(2)}")

    # Glyma (Soybean): Format conversion needed
    # Input format: Glyma.10G011600.1 (dots, mixed case)
    # Metadata format: GLYMA_10G011600 (underscores, uppercase)
    if tid.startswith(("Glyma.", "GLYMA.")):
        # Convert dots to underscores and uppercase
        normalized = tid.replace(".", "_").upper()
        candidates.append(normalized)
        # Strip version suffixes: GLYMA_10G011600_1 → GLYMA_10G011600
        if normalized.count("_") >= 2:
            parts = normalized.split("_")
            # Keep first two parts (GLYMA_10G011600)
            gene_id = "_".join(parts[:2])
            candidates.append(gene_id)

    return [c for c in dict.fromkeys(candidates) if c]


def lookup_ensembl_id(
    tid: str, base_url: Optional[str] = None, species: Optional[str] = None
) -> Optional[dict]:
    return rest_get(f"/lookup/id/{tid}", base_url=base_url, species=species)


def lookup_by_symbol(
    species: str, symbol: str, base_url: Optional[str] = None
) -> list[dict]:
    data = rest_get(
        f"/xrefs/symbol/{species}/{symbol}", base_url=base_url, species=species
    )
    return data if isinstance(data, list) else []


def lookup_xrefs_id(
    identifier: str, base_url: Optional[str] = None, species: Optional[str] = None
) -> list[dict]:
    data = rest_get(f"/xrefs/id/{identifier}", base_url=base_url, species=species)
    return data if isinstance(data, list) else []


def build_resolved_row(
    transcript_id: str, db_source: str, lookup: dict, is_ambiguous: bool
) -> dict:
    gene_id = lookup.get("Parent") or lookup.get("id") or ""
    gene_symbol = lookup.get("display_name") or ""
    return {
        "transcript_id": transcript_id,
        "db_source": db_source,
        "gene_id": gene_id,
        "gene_symbol": gene_symbol,
        "organism": lookup.get("species", ""),
        "assembly_accession": lookup.get("assembly_name", ""),
        "chrom": lookup.get("seq_region_name", ""),
        "start": lookup.get("start", ""),
        "end": lookup.get("end", ""),
        "strand": normalize_strand(lookup.get("strand")),
        "is_ambiguous": bool(is_ambiguous),
    }


def build_ambig_rows(transcript_id: str, db_source: str, chosen: dict, alts: list):
    rows = []
    for alt in alts:
        rows.append(
            {
                "transcript_id": transcript_id,
                "db_source": db_source,
                "chosen_gene_id": chosen.get("gene_id", ""),
                "alternative_gene_id": alt.get("gene_id", ""),
                "alternative_gene_symbol": alt.get("gene_symbol", ""),
                "organism": alt.get("organism", ""),
                "assembly_accession": alt.get("assembly_accession", ""),
                "chrom": alt.get("chrom", ""),
                "start": alt.get("start", ""),
                "end": alt.get("end", ""),
                "strand": alt.get("strand", "."),
            }
        )
    return rows


# ── Main ─────────────────────────────────────────────────────
log.info("Resolving external plant transcript IDs")

df = pd.read_csv(input_tsv, sep="\t")
if df.empty:
    pd.DataFrame(columns=RESOLVED_COLS).to_csv(out_resolved, sep="\t", index=False)
    pd.DataFrame(columns=AMBIG_COLS).to_csv(out_ambig, sep="\t", index=False)
    pd.DataFrame(columns=UNRESOLVED_COLS).to_csv(out_unresolved, sep="\t", index=False)
    log.info("No classified IDs found; wrote empty outputs.")
    sys.exit(0)

unresolved_input = None
try:
    if "unresolved" in snakemake.input:
        unresolved_input = snakemake.input["unresolved"]
except Exception:
    unresolved_input = None

metadata_tables = {}
for species, table_path in (cfg.get("external_metadata_tables", {}) or {}).items():
    metadata_tables[species] = load_metadata_table(table_path)

external_df = df[df["db_source"].isin(["plant"])].copy()
if unresolved_input:
    unresolved_df = pd.read_csv(unresolved_input, sep="\t")
    unresolved_ids = set(unresolved_df["transcript_id"].astype(str))
    external_df = external_df[
        external_df["transcript_id"].astype(str).isin(unresolved_ids)
    ].copy()
    log.info(f"Filtering to unresolved IDs: {len(external_df)}")
log.info(f"External IDs to resolve: {len(external_df)}")

resolved_rows = []
ambig_rows = []
unresolved_rows = []

total_ids = len(external_df)
processed_count = 0

for _, row in external_df.iterrows():
    processed_count += 1
    if processed_count % 100 == 0 or processed_count == total_ids:
        resolved_count = len(resolved_rows)
        ambig_count = len(ambig_rows)
        unresolved_count = len(unresolved_rows)
        success_rate = (
            100 * resolved_count / processed_count if processed_count > 0 else 0
        )
        log.info(
            f"Progress: {processed_count}/{total_ids} IDs ({100*processed_count/total_ids:.1f}%) | "
            f"Resolved: {resolved_count} ({success_rate:.1f}%), Ambiguous: {ambig_count}, "
            f"Unresolved: {unresolved_count}"
        )
        sys.stdout.flush()
        sys.stderr.flush()

    transcript_id = str(row["transcript_id"])
    raw_header = str(row["raw_header"])
    source_file = str(row["source_file"])
    db_source = str(row["db_source"])

    if db_source == "plant":
        resolved = None
        resolved_flag = False
        candidates = plant_candidates(transcript_id)
        species = None
        for key in PLANT_PREFIX_TO_SPECIES:
            if transcript_id.startswith(key):
                species = PLANT_PREFIX_TO_SPECIES[key][0]
                break

        if species and species in metadata_tables:
            table = metadata_tables[species]
            for cand in candidates:
                if cand in table:
                    meta = table[cand]
                    resolved = {
                        "transcript_id": transcript_id,
                        "db_source": db_source,
                        "gene_id": meta.get("gene_id", ""),
                        "gene_symbol": meta.get("gene_symbol", ""),
                        "organism": species,
                        "assembly_accession": meta.get("assembly_accession", ""),
                        "chrom": meta.get("chrom", ""),
                        "start": meta.get("start", ""),
                        "end": meta.get("end", ""),
                        "strand": normalize_strand(meta.get("strand")),
                        "is_ambiguous": False,
                    }
                    break

        # Skip REST lookups for plants (plants.rest.ensembl.org is down, metadata-only)
        if not SKIP_REST_FOR_PLANTS:
            base_url = get_rest_base_for_species(species) if species else REST_BASE
            for cand in candidates:
                lookup = lookup_ensembl_id(cand, base_url=base_url, species=species)
                if lookup:
                    resolved = build_resolved_row(
                        transcript_id, db_source, lookup, False
                    )
                    break
        if resolved:
            resolved_rows.append(resolved)
            resolved_flag = True
            continue

        # Try species-specific symbol lookup if prefix suggests a species
        if not SKIP_REST_FOR_PLANTS:
            prefix = None
            for key in PLANT_PREFIX_TO_SPECIES:
                if transcript_id.startswith(key):
                    prefix = key
                    break
            if prefix:
                for species in PLANT_PREFIX_TO_SPECIES[prefix]:
                    base_url = get_rest_base_for_species(species)
                    xrefs = []
                    for cand in candidates:
                        xrefs = lookup_by_symbol(species, cand, base_url=base_url)
                        if xrefs:
                            break
                    if not xrefs:
                        continue
                    lookups = []
                    for xref in xrefs:
                        lookup = lookup_ensembl_id(
                            xref.get("id", ""), base_url=base_url, species=species
                        )
                        if lookup:
                            lookups.append(
                                build_resolved_row(
                                    transcript_id, db_source, lookup, False
                                )
                            )
                    if lookups:
                        chosen = lookups[0]
                        resolved_rows.append(chosen)
                        if len(lookups) > 1:
                            ambig_rows.extend(
                                build_ambig_rows(
                                    transcript_id, db_source, chosen, lookups[1:]
                                )
                            )
                        resolved_flag = True
                        break

        if not resolved_flag:
            unresolved_rows.append(
                {
                    "transcript_id": transcript_id,
                    "raw_header": raw_header,
                    "source_file": source_file,
                    "reason": "plant_id_not_resolved",
                }
            )
        continue

    unresolved_rows.append(
        {
            "transcript_id": transcript_id,
            "raw_header": raw_header,
            "source_file": source_file,
            "reason": "unsupported_db_source",
        }
    )

# ── Write outputs ─────────────────────────────────────────────
res_df = pd.DataFrame(resolved_rows, columns=RESOLVED_COLS)
amb_df = pd.DataFrame(ambig_rows, columns=AMBIG_COLS)
unres_df = pd.DataFrame(unresolved_rows, columns=UNRESOLVED_COLS)

res_df.to_csv(out_resolved, sep="\t", index=False)
amb_df.to_csv(out_ambig, sep="\t", index=False)
unres_df.to_csv(out_unresolved, sep="\t", index=False)

log.info(
    f"Resolution complete: {len(res_df)} resolved, {len(amb_df)} ambiguous, {len(unres_df)} unresolved"
)

log.info("=" * 60)
log.info(f"Resolved external IDs : {len(res_df)}")
log.info(f"Ambiguous external IDs: {len(amb_df)}")
log.info(f"Unresolved external IDs: {len(unres_df)}")
log.info(f"Written external_resolved.tsv   → {out_resolved}")
log.info(f"Written external_ambiguous.tsv  → {out_ambig}")
log.info(f"Written external_unresolved.tsv → {out_unresolved}")
log.info("resolve_external_ids complete.")
