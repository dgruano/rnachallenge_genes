"""
scripts/biomart_arabidopsis_download.py
Arabidopsis thaliana BioMart downloader (Ensembl Plants)
=========================================================

Mirrors biomart_download.py but queries the Ensembl Plants BioMart endpoint
(plants.ensembl.org) which is not supported by the official Snakemake wrapper
(that wrapper hardcodes useEnsembl(), which only reaches vertebrate Ensembl).

The Arabidopsis dataset name in Ensembl Plants is "athaliana_eg_gene".
"""

import gzip
import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

if "snakemake" in dir():
    log = get_logger("biomart_arabidopsis_download", snakemake.log[0])
    release = snakemake.params.release
    output_file = snakemake.output.table
    max_retries = int(snakemake.config.get("max_retries", 3))
    retry_wait = int(snakemake.config.get("retry_wait_seconds", 5))
else:
    release = sys.argv[1] if len(sys.argv) > 1 else "60"
    output_file = sys.argv[2] if len(sys.argv) > 2 else "athaliana_biomart.tsv.gz"
    log_file = sys.argv[3] if len(sys.argv) > 3 else None
    max_retries = 3
    retry_wait = 5
    log = get_logger("biomart_arabidopsis_download", log_file)

# Ensembl Plants uses a separate BioMart host and mart name.
# Dataset: athaliana_eg_gene  (TAIR10 reference)
BIOMART_URL = "https://plants.ensembl.org/biomart/martservice"
DATASET = "athaliana_eg_gene"

BIOMART_ATTRIBUTES = [
    "ensembl_transcript_id",
    "ensembl_gene_id",
    "external_gene_name",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
]


def build_query(dataset: str, attributes: list) -> str:
    attr_xml = "\n        ".join(f'<Attribute name="{a}" />' for a in attributes)
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        "<!DOCTYPE Query>\n"
        '<Query virtualSchemaName="plants_mart" formatter="TSV" '
        'header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">\n'
        f'    <Dataset name="{dataset}" interface="default">\n'
        f"        {attr_xml}\n"
        "    </Dataset>\n"
        "</Query>"
    )


def query_biomart(max_retries: int, retry_wait: int) -> pd.DataFrame:
    query = build_query(DATASET, BIOMART_ATTRIBUTES)

    for attempt in range(1, max_retries + 1):
        try:
            log.info(
                f"Querying Ensembl Plants BioMart ({DATASET}), attempt {attempt}/{max_retries}"
            )
            resp = requests.post(BIOMART_URL, data={"query": query}, timeout=300)

            if resp.status_code != 200:
                log.warning(f"HTTP {resp.status_code}: {resp.text[:300]}")
                _maybe_retry(attempt, max_retries, retry_wait)
                continue

            lines = resp.text.strip().split("\n")
            if len(lines) < 2:
                log.warning(f"Empty response ({len(lines)} lines)")
                _maybe_retry(attempt, max_retries, retry_wait)
                continue

            header = lines[0].split("\t")
            data = [l.split("\t") for l in lines[1:] if l.strip()]
            df = pd.DataFrame(data, columns=header)
            log.info(f"Fetched {len(df)} rows")
            return df

        except requests.exceptions.Timeout:
            log.error(f"Timeout on attempt {attempt}")
            _maybe_retry(attempt, max_retries, retry_wait)
        except requests.exceptions.ConnectionError as exc:
            log.error(f"Connection error on attempt {attempt}: {exc}")
            _maybe_retry(attempt, max_retries, retry_wait)

    raise RuntimeError(f"BioMart query failed after {max_retries} attempts")


def _maybe_retry(attempt: int, max_retries: int, retry_wait: int) -> None:
    if attempt < max_retries:
        wait = retry_wait * (2 ** (attempt - 1))
        log.info(f"Retrying in {wait}s...")
        time.sleep(wait)
    else:
        raise RuntimeError("All retries exhausted")


def write_gzipped_tsv(df: pd.DataFrame, path: str) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tsv")
    df.to_csv(tmp, sep="\t", index=False)
    with open(tmp, "rb") as fi, gzip.open(out, "wb") as fo:
        fo.writelines(fi)
    tmp.unlink()
    log.info(f"Wrote {len(df)} rows → {out}")


def main():
    log.info("=== Arabidopsis thaliana BioMart Download (Ensembl Plants) ===")
    log.info(f"Release: {release}  Dataset: {DATASET}")
    log.info(f"Output: {output_file}")

    df = query_biomart(max_retries, retry_wait)

    missing = set(BIOMART_ATTRIBUTES) - set(df.columns)
    if missing:
        log.warning(f"Missing columns in response: {missing}")

    write_gzipped_tsv(df, output_file)
    log.info("=== Done ===")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log.error(f"Fatal: {exc}", exc_info=True)
        sys.exit(1)
