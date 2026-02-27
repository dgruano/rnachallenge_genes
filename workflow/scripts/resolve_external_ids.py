"""
scripts/resolve_external_ids.py
Resolve plant and WormBase transcript IDs
=========================================
Uses Ensembl REST lookups for plant IDs and parses WormBase-encoded
headers when coordinates are embedded. Remaining IDs are written to
external_unresolved.tsv for transparency.
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
MAX_RETRIES = int(cfg.get("max_retries", 3))
RETRY_WAIT = int(cfg.get("retry_wait_seconds", 5))

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

WORMBASE_COORD_RE = re.compile(
    r"^(?P<transcript>[^_]+)_wormbase:known_chromosome:(?P<assembly>WBcel\d+):"
    r"(?P<chrom>[^:]+):(?P<start>\d+):(?P<end>\d+):.*gene:(?P<gene_id>WBGene\d+)",
    re.IGNORECASE,
)

WORMBASE_ASM_TO_ORG = {
    "WBcel235": "caenorhabditis_elegans",
}

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

FLYBASE_COORD_RE = re.compile(
    r"chromosome:(?P<assembly>[^:]+):(?P<chrom>[^:]+):(?P<start>\d+):(?P<end>\d+):(?P<strand>[-_]?1)"
    r"\s+gene:(?P<gene_id>FBgn\d+).*?gene_symbol:(?P<gene_symbol>[^\s]+)",
    re.IGNORECASE,
)

SGD_ACC_RE = re.compile(r"S\d{9}")


# ── Helpers ───────────────────────────────────────────────────
def rest_get(path: str, params: Optional[dict] = None) -> Optional[dict]:
    url = f"{REST_BASE}{path}"
    headers = {"Content-Type": "application/json"}
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.get(url, params=params, headers=headers, timeout=30)
            if resp.status_code == 200:
                return resp.json()
            if resp.status_code in (400, 404):
                return None
            log.warning(
                f"  REST {path} attempt {attempt} failed: {resp.status_code} {resp.text}"
            )
        except Exception as exc:
            log.warning(f"  REST {path} attempt {attempt} failed: {exc}")
        if attempt < MAX_RETRIES:
            time.sleep(RETRY_WAIT * attempt)
    return None


def normalize_strand(value) -> str:
    if value in (1, "1", "+", "+1"):
        return "+"
    if value in (-1, "-1", "-", "-1"):
        return "-"
    return "."


def plant_candidates(tid: str) -> list[str]:
    candidates = [tid]
    if tid.endswith("_cdna"):
        candidates.append(tid[: -len("_cdna")])
    candidates.append(re.sub(r"_\d+_cdna$", "", tid))
    candidates.append(re.sub(r"_\d+$", "", tid))
    if re.search(r"\.\d+(?:\.\d+)?$", tid):
        parts = tid.split(".")
        while len(parts) > 1:
            parts = parts[:-1]
            candidates.append(".".join(parts))
    if tid.startswith("OS"):
        candidates.append("Os" + tid[2:])
        match = re.match(r"OS(\d+)T(\d+)", tid, re.IGNORECASE)
        if match:
            candidates.append(f"LOC_Os{match.group(1)}g{match.group(2)}")
            candidates.append(f"Os{match.group(1)}t{match.group(2)}-01")
            candidates.append(f"Os{match.group(1)}t{match.group(2)}")
    return [c for c in dict.fromkeys(candidates) if c]


def wormbase_candidates(tid: str) -> list[str]:
    candidates = [tid]
    candidates.append(re.sub(r"\.[a-z]?\d+$", "", tid))
    candidates.append(re.sub(r"\.[a-z]\.\d+$", "", tid))
    return [c for c in dict.fromkeys(candidates) if c]


def lookup_ensembl_id(tid: str) -> Optional[dict]:
    return rest_get(f"/lookup/id/{tid}")


def lookup_by_symbol(species: str, symbol: str) -> list[dict]:
    data = rest_get(f"/xrefs/symbol/{species}/{symbol}")
    return data if isinstance(data, list) else []


def lookup_xrefs_id(identifier: str) -> list[dict]:
    data = rest_get(f"/xrefs/id/{identifier}")
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
log.info("Resolving external (plant/wormbase) transcript IDs")

df = pd.read_csv(input_tsv, sep="\t")
if df.empty:
    pd.DataFrame(columns=RESOLVED_COLS).to_csv(out_resolved, sep="\t", index=False)
    pd.DataFrame(columns=AMBIG_COLS).to_csv(out_ambig, sep="\t", index=False)
    pd.DataFrame(columns=UNRESOLVED_COLS).to_csv(out_unresolved, sep="\t", index=False)
    log.info("No classified IDs found; wrote empty outputs.")
    sys.exit(0)

external_df = df[df["db_source"].isin(["plant", "wormbase", "flybase", "sgd"])].copy()
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
        success_rate = 100 * resolved_count / processed_count if processed_count > 0 else 0
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

    if db_source == "wormbase":
        match = WORMBASE_COORD_RE.match(transcript_id)
        if match:
            assembly = match.group("assembly")
            organism = WORMBASE_ASM_TO_ORG.get(assembly, "caenorhabditis_elegans")
            resolved_rows.append(
                {
                    "transcript_id": transcript_id,
                    "db_source": db_source,
                    "gene_id": match.group("gene_id"),
                    "gene_symbol": match.group("transcript"),
                    "organism": organism,
                    "assembly_accession": assembly,
                    "chrom": match.group("chrom"),
                    "start": int(match.group("start")),
                    "end": int(match.group("end")),
                    "strand": ".",
                    "is_ambiguous": False,
                }
            )
            continue

        # Try Ensembl REST xrefs for gene-style WormBase IDs
        species = "caenorhabditis_elegans"
        xrefs = []
        for cand in wormbase_candidates(transcript_id):
            xrefs = lookup_by_symbol(species, cand)
            if xrefs:
                break
        if xrefs:
            lookups = []
            for xref in xrefs:
                lookup = lookup_ensembl_id(xref.get("id", ""))
                if lookup:
                    lookups.append(build_resolved_row(transcript_id, db_source, lookup, False))
            if lookups:
                chosen = lookups[0]
                resolved_rows.append(chosen)
                if len(lookups) > 1:
                    ambig_rows.extend(build_ambig_rows(transcript_id, db_source, chosen, lookups[1:]))
                continue

        unresolved_rows.append(
            {
                "transcript_id": transcript_id,
                "raw_header": raw_header,
                "source_file": source_file,
                "reason": "wormbase_id_not_resolved",
            }
        )
        continue

    if db_source == "plant":
        resolved = None
        resolved_flag = False
        for cand in plant_candidates(transcript_id):
            lookup = lookup_ensembl_id(cand)
            if lookup:
                resolved = build_resolved_row(transcript_id, db_source, lookup, False)
                break
        if resolved:
            resolved_rows.append(resolved)
            resolved_flag = True
            continue

        # Try species-specific symbol lookup if prefix suggests a species
        prefix = None
        for key in PLANT_PREFIX_TO_SPECIES:
            if transcript_id.startswith(key):
                prefix = key
                break
        if prefix:
            for species in PLANT_PREFIX_TO_SPECIES[prefix]:
                xrefs = lookup_by_symbol(species, transcript_id)
                if not xrefs:
                    continue
                lookups = []
                for xref in xrefs:
                    lookup = lookup_ensembl_id(xref.get("id", ""))
                    if lookup:
                        lookups.append(build_resolved_row(transcript_id, db_source, lookup, False))
                if lookups:
                    chosen = lookups[0]
                    resolved_rows.append(chosen)
                    if len(lookups) > 1:
                        ambig_rows.extend(build_ambig_rows(transcript_id, db_source, chosen, lookups[1:]))
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

    if db_source == "flybase":
        match = FLYBASE_COORD_RE.search(raw_header)
        if match:
            strand_val = match.group("strand")
            strand = "+" if strand_val in ("1", "+1") else "-"
            resolved_rows.append(
                {
                    "transcript_id": transcript_id,
                    "db_source": db_source,
                    "gene_id": match.group("gene_id"),
                    "gene_symbol": match.group("gene_symbol"),
                    "organism": "drosophila_melanogaster",
                    "assembly_accession": match.group("assembly"),
                    "chrom": match.group("chrom"),
                    "start": int(match.group("start")),
                    "end": int(match.group("end")),
                    "strand": strand,
                    "is_ambiguous": False,
                }
            )
            continue

        lookup = lookup_ensembl_id(transcript_id)
        if lookup:
            resolved_rows.append(build_resolved_row(transcript_id, db_source, lookup, False))
            continue

        xrefs = lookup_by_symbol("drosophila_melanogaster", transcript_id)
        if not xrefs:
            xrefs = lookup_xrefs_id(transcript_id)
        if xrefs:
            lookups = []
            for xref in xrefs:
                lookup = lookup_ensembl_id(xref.get("id", ""))
                if lookup:
                    lookups.append(build_resolved_row(transcript_id, db_source, lookup, False))
            if lookups:
                chosen = lookups[0]
                resolved_rows.append(chosen)
                if len(lookups) > 1:
                    ambig_rows.extend(build_ambig_rows(transcript_id, db_source, chosen, lookups[1:]))
                continue

        unresolved_rows.append(
            {
                "transcript_id": transcript_id,
                "raw_header": raw_header,
                "source_file": source_file,
                "reason": "flybase_id_not_resolved",
            }
        )
        continue

    if db_source == "sgd":
        symbol = transcript_id.replace("_A", "").replace("_B", "")
        xrefs = lookup_by_symbol("saccharomyces_cerevisiae", symbol)
        if not xrefs:
            sgd_match = SGD_ACC_RE.search(transcript_id) or SGD_ACC_RE.search(raw_header)
            if sgd_match:
                xrefs = lookup_xrefs_id(sgd_match.group(0))
        if not xrefs and transcript_id.startswith("Q"):
            xrefs = lookup_by_symbol("saccharomyces_cerevisiae", transcript_id)

        if xrefs:
            lookups = []
            for xref in xrefs:
                lookup = lookup_ensembl_id(xref.get("id", ""))
                if lookup:
                    lookups.append(build_resolved_row(transcript_id, db_source, lookup, False))
            if lookups:
                chosen = lookups[0]
                resolved_rows.append(chosen)
                if len(lookups) > 1:
                    ambig_rows.extend(build_ambig_rows(transcript_id, db_source, chosen, lookups[1:]))
                continue

        unresolved_rows.append(
            {
                "transcript_id": transcript_id,
                "raw_header": raw_header,
                "source_file": source_file,
                "reason": "sgd_id_not_resolved",
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

log.info(f"Resolution complete: {len(res_df)} resolved, {len(amb_df)} ambiguous, {len(unres_df)} unresolved")

log.info("=" * 60)
log.info(f"Resolved external IDs : {len(res_df)}")
log.info(f"Ambiguous external IDs: {len(amb_df)}")
log.info(f"Unresolved external IDs: {len(unres_df)}")
log.info(f"Written external_resolved.tsv   → {out_resolved}")
log.info(f"Written external_ambiguous.tsv  → {out_ambig}")
log.info(f"Written external_unresolved.tsv → {out_unresolved}")
log.info("resolve_external_ids complete.")
