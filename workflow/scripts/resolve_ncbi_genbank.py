"""
scripts/resolve_ncbi_genbank.py
Second-pass NCBI resolver using EPost → EFetch → GenBank parsing
=================================================================
Consumes ncbi_ucsc_unresolved.tsv (NCBI rows only), fetches full
GenBank records via NCBIGenBankFetcher, then looks up genomic
coordinates from the NCBI Gene DB for any gene ID discovered.

Output schema matches ncbi_ucsc_resolved.tsv so it can feed
directly into merge_resolved.
"""

import re
import sys
import time
from pathlib import Path

import pandas as pd
from Bio import Entrez

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger
from ncbi_genbank_fetcher import NCBIGenBankFetcher

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("resolve_ncbi_genbank", snakemake.log[0])  # type: ignore[name-defined]
input_unresolved: str = snakemake.input.unresolved  # type: ignore[name-defined]
out_resolved: str = snakemake.output.resolved  # type: ignore[name-defined]
out_unresolved: str = snakemake.output.unresolved  # type: ignore[name-defined]
cfg = snakemake.config  # type: ignore[name-defined]

Entrez.email = cfg["ncbi_email"]
Entrez.api_key = cfg["ncbi_api_key"]

MAX_RETRIES = int(cfg.get("max_retries", 3))
RETRY_WAIT = float(cfg.get("retry_wait_seconds", 5))
BATCH_SIZE = int(cfg.get("ncbi_efetch_batch_size", 500))
GENE_BATCH = 200  # esummary gene IDs per request

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
UNRESOLVED_COLS = ["transcript_id", "db_source", "reason"]

# Real NCBI accession shape: 1-2 letters, optional underscore, 5+ digits,
# optional version. Junk like `xm_003`/`nc_201` (3 digits) fails this and, if
# eposted, makes NCBI reject the whole batch — poisoning the valid IDs with it.
ACCESSION_RE = re.compile(r"^[A-Z]{1,2}_?\d{5,}(\.\d+)?$")


# ── Gene-info lookup ──────────────────────────────────────────────────────────


def _fetch_gene_info(gene_ids: list[str]) -> dict[str, dict]:
    """Batch-esummary a list of NCBI Gene IDs → genomic coordinate dicts."""
    info: dict[str, dict] = {}
    for i in range(0, len(gene_ids), GENE_BATCH):
        chunk = gene_ids[i : i + GENE_BATCH]
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                handle = Entrez.esummary(db="gene", id=",".join(chunk))
                summaries = Entrez.read(handle)
                handle.close()
                break
            except Exception as exc:
                log.warning(f"gene esummary attempt {attempt}/{MAX_RETRIES}: {exc}")
                if attempt < MAX_RETRIES:
                    time.sleep(RETRY_WAIT * attempt)
        else:
            log.error(f"gene esummary failed for chunk at index {i}")
            continue

        for doc in summaries.get("DocumentSummarySet", {}).get("DocumentSummary", []):
            gid = doc.attributes.get("uid", "")
            loc = (doc.get("GenomicInfo") or [{}])[0]
            info[gid] = {
                "gene_id": gid,
                "gene_symbol": doc.get("Name", ""),
                "organism": doc.get("Organism", {}).get("ScientificName", ""),
                "assembly_accession": loc.get("ChrAccVer", ""),
                "chrom": loc.get("ChrLoc", ""),
                "start": int(loc.get("ChrStart", 0)),
                "end": int(loc.get("ChrStop", 0)),
                "strand": str(loc.get("ChrStrand", "")).strip() or "+",
            }
        time.sleep(0.12)
    return info


# ── Main ──────────────────────────────────────────────────────────────────────

log.info("resolve_ncbi_genbank: second-pass EPost → EFetch → GenBank resolver")

df_in = pd.read_csv(input_unresolved, sep="\t")
df_ncbi = df_in[df_in["db_source"] == "ncbi"].copy()

if df_ncbi.empty:
    log.warning("No NCBI IDs in unresolved file — writing empty outputs")
    pd.DataFrame(columns=RESOLVED_COLS).to_csv(out_resolved, sep="\t", index=False)
    pd.DataFrame(columns=UNRESOLVED_COLS).to_csv(out_unresolved, sep="\t", index=False)
    log.info("resolve_ncbi_genbank complete (nothing to do).")
    sys.exit(0)

all_ncbi_ids: list[str] = df_ncbi["transcript_id"].tolist()
accessions: list[str] = []
junk_rows: list[dict] = []
for tid in all_ncbi_ids:
    if ACCESSION_RE.match(str(tid).strip().upper()):
        accessions.append(tid)
    else:
        junk_rows.append(
            {"transcript_id": tid, "db_source": "ncbi", "reason": "invalid_accession"}
        )
if junk_rows:
    log.warning(
        f"Quarantined {len(junk_rows)} non-accession IDs before epost: "
        f"{', '.join(str(r['transcript_id']) for r in junk_rows)}"
    )

log.info(f"Fetching GenBank records for {len(accessions)} unresolved NCBI IDs")

fetcher = NCBIGenBankFetcher(
    email=cfg["ncbi_email"],
    api_key=cfg["ncbi_api_key"],
    batch_size=BATCH_SIZE,
    log=log,
    max_retries=MAX_RETRIES,
    retry_wait=RETRY_WAIT,
)
df_gb = fetcher.fetch(accessions)

# ── Collect gene IDs discovered from GenBank records ─────────────────────────
gene_ids_found = (
    df_gb.loc[df_gb["gene_id"].notna() & (df_gb["gene_id"] != "WITHDRAWN"), "gene_id"]
    .unique()
    .tolist()
)
log.info(
    f"Fetching genomic coordinates for {len(gene_ids_found)} gene IDs via esummary"
)
gene_info = _fetch_gene_info(gene_ids_found)

# ── Build resolved / unresolved rows ─────────────────────────────────────────
resolved_rows: list[dict] = []
unresolved_rows: list[dict] = list(junk_rows)

for _, row in df_gb.iterrows():
    acc = row["accession"]
    gid = row["gene_id"]
    gsym = row["gene_symbol"] or ""
    organism = row["organism"] or ""

    if not gid or gid == "WITHDRAWN":
        reason = (
            "withdrawn_or_suppressed" if gid == "WITHDRAWN" else "no_gene_id_in_genbank"
        )
        unresolved_rows.append(
            {"transcript_id": acc, "db_source": "ncbi", "reason": reason}
        )
        log.warning(f"  {acc}: {reason}")
        continue

    coords = gene_info.get(str(gid), {})
    resolved_rows.append(
        {
            "transcript_id": acc,
            "db_source": "ncbi",
            "gene_id": gid,
            "gene_symbol": coords.get("gene_symbol", gsym),
            "organism": coords.get("organism", organism),
            "assembly_accession": coords.get("assembly_accession", ""),
            "chrom": coords.get("chrom", ""),
            "start": coords.get("start", 0),
            "end": coords.get("end", 0),
            "strand": coords.get("strand", "+"),
            "is_ambiguous": False,
        }
    )

# Any accession in input not covered at all
returned = {r["transcript_id"] for r in resolved_rows} | {
    r["transcript_id"] for r in unresolved_rows
}
for acc in accessions:
    if acc not in returned:
        unresolved_rows.append(
            {
                "transcript_id": acc,
                "db_source": "ncbi",
                "reason": "not_returned_by_efetch",
            }
        )

# ── Write outputs ─────────────────────────────────────────────────────────────
df_resolved = pd.DataFrame(resolved_rows, columns=RESOLVED_COLS)
df_unresolved = pd.DataFrame(unresolved_rows, columns=UNRESOLVED_COLS)

df_resolved.to_csv(out_resolved, sep="\t", index=False)
df_unresolved.to_csv(out_unresolved, sep="\t", index=False)

log.info("=" * 60)
log.info(f"Input NCBI IDs      : {len(accessions)}")
log.info(f"Resolved            : {len(df_resolved)}")
log.info(f"Still unresolved    : {len(df_unresolved)}")
log.info(f"Written → {out_resolved}")
log.info("resolve_ncbi_genbank complete.")
