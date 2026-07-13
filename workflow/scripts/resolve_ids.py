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
out_unresolved = snakemake.output.unresolved
cfg = snakemake.config

Entrez.api_key = cfg["ncbi_api_key"]
Entrez.email = cfg["ncbi_email"]

MAX_RETRIES = int(cfg.get("max_retries", 3))
RETRY_WAIT = int(cfg.get("retry_wait_seconds", 5))
NCBI_BATCH = int(cfg.get("ncbi_batch_size", 50))
NCBI_EFETCH_BATCH = int(cfg.get("ncbi_efetch_batch_size", 400))
UCSC_REST = "https://api.genome.ucsc.edu"

# Create a persistent session for connection pooling
UCSC_SESSION = requests.Session()

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
UNRESOLVED_COLS = ["transcript_id", "db_source", "reason"]


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


def _parse_efetch_feature_table(ft_text: str) -> dict[str, dict]:
    """Parse efetch feature-table text into transcript -> {gene_id, gene_symbol}."""
    mapping: dict[str, dict] = {}
    current_tid: Optional[str] = None
    current_gene_id = ""
    current_gene_symbol = ""

    def flush_current():
        nonlocal current_tid, current_gene_id, current_gene_symbol
        if current_tid and current_gene_id:
            mapping[current_tid] = {
                "gene_id": current_gene_id,
                "gene_symbol": current_gene_symbol,
            }

    for raw in ft_text.splitlines():
        line = raw.rstrip("\n")
        if line.startswith(">Feature ref|"):
            flush_current()
            parts = line.split("|")
            current_tid = parts[1] if len(parts) > 1 else None
            current_gene_id = ""
            current_gene_symbol = ""
            continue

        if not current_tid:
            continue

        stripped = line.strip()
        if (
            stripped.startswith("db_xref")
            and "GeneID:" in stripped
            and not current_gene_id
        ):
            current_gene_id = stripped.split("GeneID:", 1)[1].strip()
        elif stripped.startswith("gene") and not current_gene_symbol:
            tokens = stripped.split()
            if len(tokens) > 1:
                current_gene_symbol = tokens[1]

    flush_current()
    return mapping


def resolve_via_efetch_ft_batch(transcript_ids: list[str]) -> dict[str, dict]:
    """Fast fallback: batch efetch feature-table parsing for suppressed/deleted IDs."""
    if not transcript_ids:
        return {}

    resolved: dict[str, dict] = {}
    for i in range(0, len(transcript_ids), NCBI_EFETCH_BATCH):
        chunk = transcript_ids[i : i + NCBI_EFETCH_BATCH]
        handle = with_retry(
            Entrez.efetch,
            db="nucleotide",
            id=",".join(chunk),
            rettype="ft",
            retmode="text",
            label=f"ncbi_efetch_ft_{i//NCBI_EFETCH_BATCH + 1}",
        )
        if handle is None:
            continue
        ft_text = handle.read()
        handle.close()

        parsed = _parse_efetch_feature_table(ft_text)
        resolved.update(parsed)

    return resolved


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
            term=" OR ".join(
                f"{acc}" for acc in batch
            ),  # Removed [ACCN] to include suppressed records
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
            # Convert Entrez StringElement to plain string
            src_gi_str = str(src_gi)
            gene_ids = [
                str(link["Id"])
                for db_links in link_set.get("LinkSetDb", [])
                if db_links["LinkName"] == "nuccore_gene"
                for link in db_links["Link"]
            ]
            gi_to_genes[src_gi_str] = gene_ids

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
                "strand": str(loc.get("ChrStrand", "")).strip() or "+",
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

        # Build TWO mappings:
        # 1. accn_to_gi: accession string → GI (for when NCBI returns matching accessions)
        # 2. gi_to_accessions: GI → list of accessions (reverse mapping for debugging)
        accn_to_gi: dict[str, str] = {}
        gi_to_accessions: dict[str, list[str]] = {}
        for doc in nucl_summaries:
            acc = doc.get("Caption", "")
            gi = str(doc.get("Gi", ""))
            acc_v = doc.get("AccessionVersion", acc)
            if acc and gi:
                accn_to_gi[acc] = gi
                accn_to_gi[acc_v] = gi
                if gi not in gi_to_accessions:
                    gi_to_accessions[gi] = []
                gi_to_accessions[gi].extend([acc, acc_v])

        log.debug(
            f"  Built accn_to_gi with {len(accn_to_gi)} keys from {len(nucl_summaries)} docs"
        )

        # Fast fallback map for suppressed/deleted records and no-gene-link records.
        # One batched efetch call per chunk is much faster than per-ID efetch.
        efetch_map = resolve_via_efetch_ft_batch(batch)
        efetch_symbols = {
            info["gene_id"]: info.get("gene_symbol", "")
            for info in efetch_map.values()
            if info.get("gene_id")
        }

        for tid in batch:
            # Try full ID first, then base ID (without version suffix)
            tid_base = tid.split(".")[0] if "." in tid else tid
            gi = accn_to_gi.get(tid) or accn_to_gi.get(tid_base)

            # FALLBACK: If ID not in accn_to_gi, try direct esearch
            if gi is None:
                log.debug(f"  NCBI: {tid} not in accn_to_gi, trying direct esearch")
                fallback_handle = with_retry(
                    Entrez.esearch,
                    db="nucleotide",
                    term=f"{tid}[ACCN]",
                    retmax=3,
                    label=f"ncbi_esearch_direct_{tid}",
                )
                if fallback_handle is not None:
                    try:
                        fallback_res = Entrez.read(fallback_handle)
                        fallback_handle.close()
                        fallback_gis = fallback_res.get("IdList", [])
                        if fallback_gis:
                            gi = str(fallback_gis[0])
                            log.debug(
                                f"  NCBI: {tid} found via direct search → GI={gi}"
                            )

                            # Need to get gene links for this newly found GI
                            link_handle = with_retry(
                                Entrez.elink,
                                dbfrom="nucleotide",
                                db="gene",
                                id=gi,
                                label=f"ncbi_elink_fallback_{gi}",
                            )
                            if link_handle is not None:
                                try:
                                    link_res = Entrez.read(link_handle)
                                    link_handle.close()
                                    for link_set in link_res:
                                        gene_ids = [
                                            str(link["Id"])
                                            for db_links in link_set.get(
                                                "LinkSetDb", []
                                            )
                                            if db_links["LinkName"] == "nuccore_gene"
                                            for link in db_links["Link"]
                                        ]
                                        if gene_ids:
                                            gi_to_genes[gi] = gene_ids
                                            log.debug(
                                                f"  NCBI: linked GI={gi} → {len(gene_ids)} genes"
                                            )
                                except Exception as exc:
                                    log.warning(
                                        f"  NCBI: elink fallback failed for GI={gi}: {exc}"
                                    )
                    except Exception as exc:
                        log.warning(f"  NCBI: direct esearch failed for {tid}: {exc}")

            # Try base ID without version if still not found
            if gi is None and tid != tid_base:
                log.debug(f"  NCBI: trying base ID {tid_base}")
                fallback_handle = with_retry(
                    Entrez.esearch,
                    db="nucleotide",
                    term=f"{tid_base}[ACCN]",
                    retmax=3,
                    label=f"ncbi_esearch_base_{tid_base}",
                )
                if fallback_handle is not None:
                    try:
                        fallback_res = Entrez.read(fallback_handle)
                        fallback_handle.close()
                        fallback_gis = fallback_res.get("IdList", [])
                        if fallback_gis:
                            gi = str(fallback_gis[0])
                            log.debug(f"  NCBI: {tid_base} found → GI={gi}")

                            # Get gene links for newly found GI
                            link_handle = with_retry(
                                Entrez.elink,
                                dbfrom="nucleotide",
                                db="gene",
                                id=gi,
                                label=f"ncbi_elink_base_{gi}",
                            )
                            if link_handle is not None:
                                try:
                                    link_res = Entrez.read(link_handle)
                                    link_handle.close()
                                    for link_set in link_res:
                                        gene_ids = [
                                            str(link["Id"])
                                            for db_links in link_set.get(
                                                "LinkSetDb", []
                                            )
                                            if db_links["LinkName"] == "nuccore_gene"
                                            for link in db_links["Link"]
                                        ]
                                        if gene_ids:
                                            gi_to_genes[gi] = gene_ids
                                            log.info(
                                                f"  NCBI: {tid} resolved via base ID {tid_base} → {len(gene_ids)} genes"
                                            )
                                except Exception as exc:
                                    log.warning(
                                        f"  NCBI: elink base fallback failed for GI={gi}: {exc}"
                                    )
                    except Exception as exc:
                        log.warning(
                            f"  NCBI: base ID esearch failed for {tid_base}: {exc}"
                        )

            if gi is None:
                gene_ids_for_tid = []
                efetch_info = efetch_map.get(tid)
                if efetch_info and efetch_info.get("gene_id"):
                    synthetic_gid = efetch_info["gene_id"]
                    gene_ids_for_tid = [synthetic_gid]
                    gi = "efetch"
                    log.info(
                        f"  NCBI: {tid} resolved via batched efetch fallback (suppressed record) → gene {synthetic_gid}"
                    )
                else:
                    log.warning(
                        f"  NCBI: could not resolve {tid} even via efetch fallback"
                    )
                    continue
            else:
                # Try gene link for the found GI
                gene_ids_for_tid = gi_to_genes.get(gi, [])
            if not gene_ids_for_tid:
                efetch_info = efetch_map.get(tid)
                if efetch_info and efetch_info.get("gene_id"):
                    synthetic_gid = efetch_info["gene_id"]
                    gene_ids_for_tid = [synthetic_gid]
                    log.info(
                        f"  NCBI: {tid} resolved via batched efetch fallback → gene {synthetic_gid}"
                    )
                if not gene_ids_for_tid:
                    log.warning(
                        f"  NCBI: no gene link for transcript {tid} (GI={gi}) even after efetch fallback"
                    )
                    continue

            is_ambiguous = len(gene_ids_for_tid) > 1
            if is_ambiguous:
                log.info(
                    f"  Ambiguous: {tid} → {len(gene_ids_for_tid)} genes, picking primary"
                )

            primary_gid = gene_ids_for_tid[0]

            # If gene info not yet fetched (e.g., from fallback), fetch it now
            if primary_gid not in gene_info:
                log.debug(
                    f"  Fetching gene info for newly discovered gene {primary_gid}"
                )
                missing_gene_ids = [
                    gid for gid in gene_ids_for_tid if gid not in gene_info
                ]
                if missing_gene_ids:
                    gene_fetch_handle = with_retry(
                        Entrez.esummary,
                        db="gene",
                        id=",".join(missing_gene_ids),
                        label=f"ncbi_gene_esummary_fallback",
                    )
                    if gene_fetch_handle is not None:
                        try:
                            gene_fetch_summaries = Entrez.read(gene_fetch_handle)
                            gene_fetch_handle.close()
                            for doc in gene_fetch_summaries.get(
                                "DocumentSummarySet", {}
                            ).get("DocumentSummary", []):
                                gid = doc.attributes.get("uid", "")
                                genomic_info = doc.get("GenomicInfo", [{}])
                                loc = genomic_info[0] if genomic_info else {}
                                gene_info[gid] = {
                                    "gene_id": gid,
                                    "gene_symbol": doc.get("Name", ""),
                                    "organism": doc.get("Organism", {}).get(
                                        "ScientificName", ""
                                    ),
                                    "assembly_accession": loc.get("ChrAccVer", ""),
                                    "chrom": loc.get("ChrLoc", ""),
                                    "start": int(loc.get("ChrStart", 0)),
                                    "end": int(loc.get("ChrStop", 0)),
                                    "strand": str(loc.get("ChrStrand", "")).strip()
                                    or "+",
                                }
                        except Exception as exc:
                            log.warning(
                                f"  Failed to fetch gene info for {missing_gene_ids}: {exc}"
                            )

            primary = gene_info.get(primary_gid, {})
            if not primary:
                primary = {
                    "gene_id": primary_gid,
                    "gene_symbol": efetch_symbols.get(primary_gid, ""),
                    "organism": "",
                    "assembly_accession": "",
                    "chrom": "",
                    "start": 0,
                    "end": 0,
                    "strand": "+",
                }

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
                resp = UCSC_SESSION.get(u, timeout=30)
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

# Species-first routing: this resolver handles only RefSeq/NCBI and UCSC.
df_cls = df_cls[df_cls["db_source"].isin(["ncbi", "ucsc"])]
log.info(f"IDs to resolve (NCBI + UCSC): {len(df_cls)}")

all_resolved: list[dict] = []
all_ambig: list[dict] = []
all_unresolved: list[dict] = []

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

    resolved_ids = {r["transcript_id"] for r in res}
    not_found = [tid for tid in ids if tid not in resolved_ids]
    reason = "not_found_in_ncbi" if db_source == "ncbi" else "not_found_in_ucsc"
    all_unresolved.extend(
        [
            {"transcript_id": tid, "db_source": db_source, "reason": reason}
            for tid in not_found
        ]
    )

    log.info(
        f"  → {len(res)} resolved, {len(amb)} ambiguous alternatives, {len(not_found)} not found"
    )
    all_resolved.extend(res)
    all_ambig.extend(amb)

df_resolved = pd.DataFrame(all_resolved, columns=RESOLVED_COLS)
df_ambig = pd.DataFrame(all_ambig, columns=AMBIG_COLS)
df_unresolved = pd.DataFrame(all_unresolved, columns=UNRESOLVED_COLS)

df_resolved.to_csv(out_resolved, sep="\t", index=False)
df_ambig.to_csv(out_ambig, sep="\t", index=False)
df_unresolved.to_csv(out_unresolved, sep="\t", index=False)

log.info("=" * 60)
log.info(f"NCBI + UCSC input   : {len(df_cls)}")
log.info(f"Resolved            : {len(df_resolved)}")
log.info(f"Matched but not found: {len(df_unresolved)}")
log.info(f"Ambiguous alts      : {len(df_ambig)}")
log.info(f"Written → {out_resolved}")
log.info("resolve_ids complete.")
