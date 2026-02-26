"""
scripts/resolve_ids.py
Resolve NCBI + UCSC Transcript IDs  (Ensembl handled via BioMart wrapper)
=========================================================================
Queries NCBI Entrez (esearch → elink → esummary) and the UCSC REST API
to map transcript IDs → gene + genomic coordinates.

Ensembl IDs are intentionally SKIPPED here; they are resolved by the
separate BioMart sub-DAG:
  detect_ensembl_species → biomart_lookup (wrapper) → join_ensembl_results

Output schema (ncbi_ucsc_resolved.tsv) — same unified schema as Ensembl:
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
UCSC_REST = "https://api.genome.ucsc.edu"

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
# NCBI resolver (unchanged from v1)
# ════════════════════════════════════════════════════════════
def resolve_ncbi_batch(transcript_ids: list[str]) -> tuple[list[dict], list[dict]]:
    resolved_rows = []
    ambig_rows = []

    log.info(f"  NCBI: resolving {len(transcript_ids)} IDs in batches of {NCBI_BATCH}")

    for i in range(0, len(transcript_ids), NCBI_BATCH):
        batch = transcript_ids[i : i + NCBI_BATCH]
        log.debug(
            f"  NCBI batch {i//NCBI_BATCH + 1}: {batch[:3]}{'...' if len(batch)>3 else ''}"
        )

        handle = with_retry(
            Entrez.esearch,
            db="nucleotide",
            term=" OR ".join(f"{acc}[ACCN]" for acc in batch),
            retmax=len(batch) * 5,
            label="ncbi_esearch",
        )
        if handle is None:
            continue
        search_res = Entrez.read(handle)
        handle.close()

        gis = search_res.get("IdList", [])
        if not gis:
            log.warning(f"  NCBI esearch returned no GIs for batch at index {i}")
            continue

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

        all_gene_ids = list({g for genes in gi_to_genes.values() for g in genes})
        if not all_gene_ids:
            log.warning("  No gene links found for this NCBI batch")
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
            acc_v = doc.get("AccessionVersion", acc)
            if acc and gi:
                accn_to_gi[acc] = gi
                accn_to_gi[acc_v] = gi

        for tid in batch:
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
                    f"  Ambiguous: {tid} → {len(gene_ids_for_tid)} genes, picking primary"
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

        time.sleep(0.12)

    return resolved_rows, ambig_rows


# ════════════════════════════════════════════════════════════
# UCSC resolver (unchanged from v1)
# ════════════════════════════════════════════════════════════
def resolve_ucsc_batch(transcript_ids: list[str]) -> tuple[list[dict], list[dict]]:
    resolved_rows = []
    ambig_rows = []

    log.info(f"  UCSC: resolving {len(transcript_ids)} IDs")

    UCSC_ASSEMBLIES = ["hg38", "hg19", "mm39", "mm10", "rn7", "dm6", "danRer11"]

    for tid in transcript_ids:
        found = False
        for assembly in UCSC_ASSEMBLIES:
            url = (
                f"{UCSC_REST}/getData/track"
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

            resolved_rows.append(
                {
                    "transcript_id": tid,
                    "db_source": "ucsc",
                    "gene_id": hit.get("name2", hit.get("name", "")),
                    "gene_symbol": hit.get("name2", ""),
                    "organism": "",
                    "assembly_accession": assembly,
                    "chrom": hit.get("chrom", ""),
                    "start": int(hit.get("txStart", 0)),
                    "end": int(hit.get("txEnd", 0)),
                    "strand": hit.get("strand", "+"),
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
            break

        if not found:
            log.warning(f"  UCSC: {tid} not found in any tested assembly")

        time.sleep(0.05)

    return resolved_rows, ambig_rows


# ── Main ─────────────────────────────────────────────────────
log.info("resolve_ids: resolving NCBI and UCSC transcript IDs")
log.info("NOTE: Ensembl IDs are resolved separately via the BioMart wrapper sub-DAG")

df_cls = pd.read_csv(input_tsv, sep="\t")

# Explicitly exclude Ensembl — handled by BioMart pipeline
df_cls = df_cls[df_cls["db_source"] != "ensembl"]
log.info(f"IDs to resolve (NCBI + UCSC): {len(df_cls)}")

all_resolved: list[dict] = []
all_ambig: list[dict] = []

for db_source, group in df_cls.groupby("db_source"):
    ids = group["transcript_id"].tolist()
    log.info(f"Resolving {len(ids)} {db_source.upper()} transcripts...")

    if db_source == "ncbi":
        res, amb = resolve_ncbi_batch(ids)
    elif db_source == "ucsc":
        res, amb = resolve_ucsc_batch(ids)
    else:
        log.warning(f"No resolver for db_source={db_source!r} — skipping")
        continue

    log.info(f"  → {len(res)} resolved, {len(amb)} ambiguous alternatives")
    all_resolved.extend(res)
    all_ambig.extend(amb)

df_resolved = pd.DataFrame(all_resolved, columns=RESOLVED_COLS)
df_ambig = pd.DataFrame(all_ambig, columns=AMBIG_COLS)

df_resolved.to_csv(out_resolved, sep="\t", index=False)
df_ambig.to_csv(out_ambig, sep="\t", index=False)

log.info("=" * 60)
log.info(f"NCBI + UCSC input   : {len(df_cls)}")
log.info(f"Resolved            : {len(df_resolved)}")
log.info(f"Failed              : {len(df_cls) - len(df_resolved)}")
log.info(f"Ambiguous alts      : {len(df_ambig)}")
log.info(f"Written → {out_resolved}")
log.info("resolve_ids complete.")
