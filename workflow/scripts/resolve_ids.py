"""
scripts/resolve_ids.py
Stage 2 — Resolve Transcript IDs to Gene Coordinates
=====================================================
For each classified transcript ID, queries the appropriate API:
  - NCBI  : Entrez efetch / esummary (Biopython Entrez)
  - Ensembl: REST API /lookup/id (batch)
  - UCSC  : UCSC REST API /list/tracks + /getData/track

For ambiguous IDs (multiple gene mappings), picks the primary /
canonical entry and logs all alternatives to ambiguous.tsv.

Output schema (resolved_ids.tsv)
---------------------------------
transcript_id | db_source | gene_id | gene_symbol | organism |
assembly_accession | chrom | start | end | strand | is_ambiguous
"""

import json
import sys
import time
from pathlib import Path
from typing import Optional

import pandas as pd
import requests
from Bio import Entrez

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("resolve_ids", snakemake.log[0])
input_tsv = snakemake.input.classified
out_resolved = snakemake.output.resolved
out_ambig = snakemake.output.ambiguous
cfg = snakemake.config

Entrez.api_key = cfg["ncbi_api_key"]
Entrez.email = cfg["ncbi_email"]

MAX_RETRIES = int(cfg.get("max_retries", 3))
RETRY_WAIT = int(cfg.get("retry_wait_seconds", 5))
NCBI_BATCH = int(cfg.get("ncbi_batch_size", 50))
ENSEMBL_BATCH = int(cfg.get("ensembl_batch_size", 50))
ENSEMBL_REST_URL = "https://rest.ensembl.org"
UCSC_REST_URL = "https://api.genome.ucsc.edu"

# ── Unified output columns ────────────────────────────────────
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


# ── Retry helper ─────────────────────────────────────────────
def with_retry(fn, *args, label="request", **kwargs):
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            return fn(*args, **kwargs)
        except Exception as exc:
            log.warning(f"  [{label}] attempt {attempt}/{MAX_RETRIES} failed: {exc}")
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_WAIT * attempt)
    log.error(f"  [{label}] all {MAX_RETRIES} attempts failed — skipping")
    return None


# ════════════════════════════════════════════════════════════
# NCBI resolver
# ════════════════════════════════════════════════════════════
def resolve_ncbi_batch(transcript_ids: list[str]) -> tuple[list[dict], list[dict]]:
    """
    Resolve a batch of NCBI RefSeq transcript IDs.
    Returns (resolved_rows, ambiguous_rows).
    """
    resolved_rows = []
    ambig_rows = []

    log.info(f"  NCBI: resolving {len(transcript_ids)} IDs in batches of {NCBI_BATCH}")

    for i in range(0, len(transcript_ids), NCBI_BATCH):
        batch = transcript_ids[i : i + NCBI_BATCH]
        log.debug(
            f"  NCBI batch {i//NCBI_BATCH + 1}: {batch[:3]}{'...' if len(batch)>3 else ''}"
        )

        # Step 1: accession → GI (nucleotide db)
        handle = with_retry(
            Entrez.esearch,
            db="nucleotide",
            term=" OR ".join(f"{acc}[ACCN]" for acc in batch),
            retmax=len(batch) * 5,
            label="ncbi_esearch",
        )
        if handle is None:
            for tid in batch:
                log.warning(f"  NCBI esearch failed for batch containing {tid}")
            continue
        search_res = Entrez.read(handle)
        handle.close()

        gis = search_res.get("IdList", [])
        if not gis:
            log.warning(
                f"  NCBI esearch returned no GIs for batch starting at index {i}"
            )
            continue

        # Step 2: GI → gene link
        link_handle = with_retry(
            Entrez.elink,
            dbfrom="nucleotide",
            db="gene",
            id=",".join(gis),
            label="ncbi_elink",
        )
        if link_handle is None:
            continue
        link_results = Entrez.read(link_handle)
        link_handle.close()

        # Build GI → gene_ids mapping
        gi_to_genes: dict[str, list[str]] = {}
        for link_set in link_results:
            src_gi = link_set["IdList"][0] if link_set["IdList"] else None
            if src_gi is None:
                continue
            gene_ids = [
                link["Id"]
                for db_links in link_set.get("LinkSetDb", [])
                if db_links["LinkName"] == "nuccore_gene"
                for link in db_links["Link"]
            ]
            gi_to_genes[src_gi] = gene_ids

        # Step 3: fetch gene summaries
        all_gene_ids = list({g for genes in gi_to_genes.values() for g in genes})
        if not all_gene_ids:
            log.warning("  No gene links found for this batch")
            continue

        gene_summary_handle = with_retry(
            Entrez.esummary,
            db="gene",
            id=",".join(all_gene_ids),
            label="ncbi_gene_esummary",
        )
        if gene_summary_handle is None:
            continue
        gene_summaries = Entrez.read(gene_summary_handle)
        gene_summary_handle.close()

        gene_info: dict[str, dict] = {}
        for doc in gene_summaries.get("DocumentSummarySet", {}).get(
            "DocumentSummary", []
        ):
            gid = doc.attributes.get("uid", "")
            genomic_info = doc.get("GenomicInfo", [{}])
            loc = genomic_info[0] if genomic_info else {}
            gene_info[gid] = {
                "gene_id": gid,
                "gene_symbol": doc.get("Name", ""),
                "organism": doc.get("Organism", {}).get("ScientificName", ""),
                "assembly_accession": loc.get("ChrAccVer", ""),
                "chrom": loc.get("ChrLoc", ""),
                "start": int(loc.get("ChrStart", 0)),
                "end": int(loc.get("ChrStop", 0)),
                "strand": "+" if loc.get("ExonCount", 0) >= 0 else "-",
            }

        # Step 4: fetch nucleotide summaries to map accession → GI
        nucl_handle = with_retry(
            Entrez.esummary,
            db="nucleotide",
            id=",".join(gis),
            label="ncbi_nucl_esummary",
        )
        if nucl_handle is None:
            continue
        nucl_summaries = Entrez.read(nucl_handle)
        nucl_handle.close()

        accn_to_gi: dict[str, str] = {}
        for doc in nucl_summaries:
            acc = doc.get("Caption", "")
            gi = str(doc.get("Gi", ""))
            if acc and gi:
                accn_to_gi[acc] = gi
                # Also map versioned accession
                acc_v = doc.get("AccessionVersion", acc)
                accn_to_gi[acc_v] = gi

        # Step 5: assemble per-transcript rows
        for tid in batch:
            # Try both versioned and unversioned accession
            gi = accn_to_gi.get(tid) or accn_to_gi.get(tid.split(".")[0])
            if gi is None:
                log.warning(f"  NCBI: could not map {tid} to a GI")
                continue

            gene_ids_for_tid = gi_to_genes.get(gi, [])
            if not gene_ids_for_tid:
                log.warning(f"  NCBI: no gene link for transcript {tid} (GI={gi})")
                continue

            is_ambiguous = len(gene_ids_for_tid) > 1
            if is_ambiguous:
                log.info(
                    f"  Ambiguous: {tid} maps to {len(gene_ids_for_tid)} genes — picking primary"
                )

            primary_gid = gene_ids_for_tid[0]
            primary = gene_info.get(primary_gid, {})

            resolved_rows.append(
                {
                    "transcript_id": tid,
                    "db_source": "ncbi",
                    "gene_id": primary.get("gene_id", ""),
                    "gene_symbol": primary.get("gene_symbol", ""),
                    "organism": primary.get("organism", ""),
                    "assembly_accession": primary.get("assembly_accession", ""),
                    "chrom": primary.get("chrom", ""),
                    "start": primary.get("start", 0),
                    "end": primary.get("end", 0),
                    "strand": primary.get("strand", "+"),
                    "is_ambiguous": is_ambiguous,
                }
            )

            for alt_gid in gene_ids_for_tid[1:]:
                alt = gene_info.get(alt_gid, {})
                ambig_rows.append(
                    {
                        "transcript_id": tid,
                        "db_source": "ncbi",
                        "chosen_gene_id": primary_gid,
                        "alternative_gene_id": alt.get("gene_id", ""),
                        "alternative_gene_symbol": alt.get("gene_symbol", ""),
                        "organism": alt.get("organism", ""),
                        "assembly_accession": alt.get("assembly_accession", ""),
                        "chrom": alt.get("chrom", ""),
                        "start": alt.get("start", 0),
                        "end": alt.get("end", 0),
                        "strand": alt.get("strand", "+"),
                    }
                )

        time.sleep(0.12)  # ~8 req/sec with API key

    return resolved_rows, ambig_rows


# ════════════════════════════════════════════════════════════
# Ensembl resolver
# ════════════════════════════════════════════════════════════
def _ensembl_post(endpoint: str, payload: dict, label: str) -> Optional[dict]:
    """POST to Ensembl REST and return JSON."""
    url = f"{ENSEMBL_REST_URL}{endpoint}"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    def _do():
        resp = requests.post(url, headers=headers, data=json.dumps(payload), timeout=60)
        resp.raise_for_status()
        return resp.json()

    return with_retry(_do, label=label)


def _ensembl_get(endpoint: str, label: str) -> Optional[dict]:
    """GET from Ensembl REST and return JSON."""
    url = f"{ENSEMBL_REST_URL}{endpoint}"
    headers = {"Accept": "application/json"}

    def _do():
        resp = requests.get(url, headers=headers, timeout=30)
        resp.raise_for_status()
        return resp.json()

    return with_retry(_do, label=label)


def resolve_ensembl_batch(transcript_ids: list[str]) -> tuple[list[dict], list[dict]]:
    resolved_rows = []
    ambig_rows = []

    log.info(
        f"  Ensembl: resolving {len(transcript_ids)} IDs in batches of {ENSEMBL_BATCH}"
    )

    for i in range(0, len(transcript_ids), ENSEMBL_BATCH):
        batch = transcript_ids[i : i + ENSEMBL_BATCH]
        log.debug(
            f"  Ensembl batch {i//ENSEMBL_BATCH + 1}: {batch[:3]}{'...' if len(batch)>3 else ''}"
        )

        # Batch lookup
        result = _ensembl_post(
            "/lookup/id",
            {"ids": batch, "expand": 1},
            label="ensembl_lookup",
        )
        if result is None:
            log.warning(f"  Ensembl batch lookup failed for batch at index {i}")
            continue

        for tid in batch:
            data = result.get(tid)
            if data is None or "error" in data:
                log.warning(f"  Ensembl: {tid} not found — {data}")
                continue

            # data is transcript-level; get parent gene
            parent_id = data.get("Parent") or data.get("gene_id")
            if not parent_id:
                log.warning(f"  Ensembl: no parent gene for {tid}")
                continue

            gene_data = _ensembl_get(
                f"/lookup/id/{parent_id}?expand=0", label=f"ensembl_gene_{parent_id}"
            )
            if gene_data is None:
                log.warning(f"  Ensembl: could not fetch gene {parent_id}")
                continue

            # Strand: Ensembl uses 1 / -1
            strand = "+" if gene_data.get("strand", 1) == 1 else "-"
            # Assembly: infer from species + coord_system_version
            assembly = gene_data.get("assembly_name", "")

            resolved_rows.append(
                {
                    "transcript_id": tid,
                    "db_source": "ensembl",
                    "gene_id": parent_id,
                    "gene_symbol": gene_data.get("display_name", ""),
                    "organism": gene_data.get("species", "").replace("_", " "),
                    "assembly_accession": assembly,
                    "chrom": gene_data.get("seq_region_name", ""),
                    "start": int(gene_data.get("start", 0)),
                    "end": int(gene_data.get("end", 0)),
                    "strand": strand,
                    "is_ambiguous": False,  # Ensembl transcript→gene is 1:1
                }
            )

        time.sleep(0.1)  # be polite to Ensembl REST

    return resolved_rows, ambig_rows


# ════════════════════════════════════════════════════════════
# UCSC resolver
# ════════════════════════════════════════════════════════════
def resolve_ucsc_batch(transcript_ids: list[str]) -> tuple[list[dict], list[dict]]:
    """
    Resolve UCSC transcript IDs via the UCSC REST API.
    UCSC IDs encode the assembly (e.g. uc001aaa.3 → hg19).
    We query the knownGene / refGene track for coordinates.
    """
    resolved_rows = []
    ambig_rows = []

    log.info(f"  UCSC: resolving {len(transcript_ids)} IDs individually")

    # UCSC IDs don't encode species reliably; we try common assemblies
    UCSC_ASSEMBLIES = ["hg38", "hg19", "mm39", "mm10", "rn7", "dm6", "danRer11"]

    for tid in transcript_ids:
        found = False
        for assembly in UCSC_ASSEMBLIES:
            url = (
                f"{UCSC_REST_URL}/getData/track"
                f"?genome={assembly}&track=knownGene&name={tid}"
            )

            def _do(u=url):
                resp = requests.get(u, timeout=30)
                resp.raise_for_status()
                return resp.json()

            data = with_retry(_do, label=f"ucsc_{tid}_{assembly}")
            if data is None:
                continue

            hits = data.get("knownGene", [])
            if not hits:
                continue

            hit = hits[0]
            is_ambiguous = len(hits) > 1
            if is_ambiguous:
                log.info(
                    f"  UCSC: {tid} → {len(hits)} hits in {assembly} — picking first"
                )

            strand = hit.get("strand", "+")
            resolved_rows.append(
                {
                    "transcript_id": tid,
                    "db_source": "ucsc",
                    "gene_id": hit.get("name2", hit.get("name", "")),
                    "gene_symbol": hit.get("name2", ""),
                    "organism": "",  # UCSC REST doesn't return species name directly
                    "assembly_accession": assembly,
                    "chrom": hit.get("chrom", ""),
                    "start": int(hit.get("txStart", 0)),
                    "end": int(hit.get("txEnd", 0)),
                    "strand": strand,
                    "is_ambiguous": is_ambiguous,
                }
            )

            for alt in hits[1:]:
                ambig_rows.append(
                    {
                        "transcript_id": tid,
                        "db_source": "ucsc",
                        "chosen_gene_id": hit.get("name2", ""),
                        "alternative_gene_id": alt.get("name2", ""),
                        "alternative_gene_symbol": alt.get("name2", ""),
                        "organism": "",
                        "assembly_accession": assembly,
                        "chrom": alt.get("chrom", ""),
                        "start": int(alt.get("txStart", 0)),
                        "end": int(alt.get("txEnd", 0)),
                        "strand": alt.get("strand", "+"),
                    }
                )

            found = True
            break  # stop trying assemblies once found

        if not found:
            log.warning(f"  UCSC: {tid} not found in any tested assembly")

        time.sleep(0.05)

    return resolved_rows, ambig_rows


# ── Main ─────────────────────────────────────────────────────
log.info("Stage 2: Resolving transcript IDs to gene coordinates")

df_cls = pd.read_csv(input_tsv, sep="\t")
log.info(f"Loaded {len(df_cls)} classified transcript IDs")

all_resolved: list[dict] = []
all_ambig: list[dict] = []

for db_source, group in df_cls.groupby("db_source"):
    ids = group["transcript_id"].tolist()
    log.info(f"Resolving {len(ids)} {db_source.upper()} transcripts...")

    if db_source == "ncbi":
        res, amb = resolve_ncbi_batch(ids)
    elif db_source == "ensembl":
        res, amb = resolve_ensembl_batch(ids)
    elif db_source == "ucsc":
        res, amb = resolve_ucsc_batch(ids)
    else:
        log.warning(f"No resolver for db_source={db_source!r} — skipping")
        continue

    log.info(f"  → {len(res)} resolved, {len(amb)} ambiguous alternatives recorded")
    all_resolved.extend(res)
    all_ambig.extend(amb)

df_resolved = pd.DataFrame(all_resolved, columns=RESOLVED_COLS)
df_ambig = pd.DataFrame(all_ambig, columns=AMBIG_COLS)

df_resolved.to_csv(out_resolved, sep="\t", index=False)
df_ambig.to_csv(out_ambig, sep="\t", index=False)

# ── Summary ──────────────────────────────────────────────────
total_in = len(df_cls)
total_out = len(df_resolved)
missed = total_in - total_out

log.info("=" * 60)
log.info(f"Input transcripts            : {total_in}")
log.info(f"Successfully resolved        : {total_out}")
log.info(f"Failed to resolve            : {missed}")
log.info(f"Ambiguous (alternatives)     : {len(df_ambig)}")
log.info(f"Written resolved    → {out_resolved}")
log.info(f"Written ambiguous   → {out_ambig}")
log.info("Stage 2 complete.")
