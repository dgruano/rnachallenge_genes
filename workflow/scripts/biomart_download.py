"""
scripts/biomart_download.py
Custom Python-based Ensembl BioMart downloader for vertebrate species
=====================================================================

Replaces the snakemake-wrappers R-based biomart wrapper to avoid variable name
bugs in the R code. This script directly queries the Ensembl BioMart REST API
via XML queries to download full transcript annotation tables.

Features:
  - Queries BioMart to fetch all transcripts for a species (no pre-filtering)
  - Handles network errors with exponential backoff retries
  - Outputs gzipped TSV for efficient storage and caching
  - Supports all vertebrate species in Ensembl
  - Logs detailed progress and statistics
"""

import gzip
import sys
import time
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake integration ──────────────────────────────────────
# When called from Snakemake, these globals are available:
#   snakemake.params.species       (e.g., "homo_sapiens")
#   snakemake.params.release       (e.g., "115")
#   snakemake.output.table         (path to output .tsv.gz file)
#   snakemake.log[0]               (log file path)
#   snakemake.config               (pipeline config)
#
# When called standalone (for testing):
#   python biomart_download.py <species> <release> <output_file> [log_file]

if "snakemake" in dir():
    # Running under Snakemake
    log = get_logger("biomart_download", snakemake.log[0])
    species = snakemake.params.species
    release = snakemake.params.release
    output_file = snakemake.output.table
    max_retries = int(snakemake.config.get("max_retries", 3))
    retry_wait = int(snakemake.config.get("retry_wait_seconds", 5))
    on_failure = snakemake.config.get("biomart_on_failure", "fail")
else:
    # Running standalone
    if len(sys.argv) < 4:
        print("Usage: biomart_download.py <species> <release> <output_file> [log_file]")
        print("Example: biomart_download.py homo_sapiens 115 out.tsv.gz")
        sys.exit(1)

    species = sys.argv[1]
    release = sys.argv[2]
    output_file = sys.argv[3]
    log_file = sys.argv[4] if len(sys.argv) > 4 else None
    max_retries = 3
    retry_wait = 5
    on_failure = "fail"

    log = get_logger("biomart_download", log_file)

# ── Ensembl BioMart configuration ──────────────────────────────
BIOMART_URL = "https://www.ensembl.org/biomart/martservice"

# Mapping from Ensembl species names to BioMart dataset names
# Reference: https://www.ensembl.org/Help/Powered?action=push;target=FTPtop
SPECIES_TO_DATASET = {
    "homo_sapiens": "hsapiens_gene_ensembl",
    "mus_musculus": "mmusculus_gene_ensembl",
    "rattus_norvegicus": "rnorvegicus_gene_ensembl",
    "danio_rerio": "drerio_gene_ensembl",
    "gallus_gallus": "ggallus_gene_ensembl",
    "sus_scrofa": "sscrofa_gene_ensembl",
    "bos_taurus": "btaurus_gene_ensembl",
    "canis_lupus_familiaris": "cfamiliaris_gene_ensembl",
    "equus_caballus": "ecaballus_gene_ensembl",
    "pan_troglodytes": "ptroglodytes_gene_ensembl",
    "macaca_mulatta": "mmulatta_gene_ensembl",
    "ovis_aries": "oaries_gene_ensembl",
    "capra_hircus": "chircus_gene_ensembl",
    "chlorocebus_aethiops": "caethiops_gene_ensembl",
    "papio_anubis": "panubis_gene_ensembl",
    "mandrillus_leucophaeus": "mleucophaus_gene_ensembl",
    "heterocephalus_glaber": "hglaber_gene_ensembl",
    "ochotona_princeps": "oprinceps_gene_ensembl",
    "oryctolagus_cuniculus": "ocuniculus_gene_ensembl",
    "cavia_porcellus": "cporcellus_gene_ensembl",
    "tursiops_truncatus": "ttruncatus_gene_ensembl",
    "monodelphis_domestica": "mdomestica_gene_ensembl",
    "sarcophilus_harrisii": "sharrisii_gene_ensembl",
    "ornithorhynchus_anatinus": "oanatinus_gene_ensembl",
    "taeniopygia_guttata": "tguttata_gene_ensembl",
    "parus_major": "pmajor_gene_ensembl",
    "anolis_carolinensis": "acarolinensis_gene_ensembl",
    "pogona_vitticeps": "pvitticeps_gene_ensembl",
    "python_bivittatus": "pbivittatus_gene_ensembl",
    "xenopus_tropicalis": "xtropicalis_gene_ensembl",
    "xenopus_laevis": "xlaevis_gene_ensembl",
    "latimeria_chalumnae": "lchalumnae_gene_ensembl",
    "petromyzon_marinus": "pmarinus_gene_ensembl",
    "callorhinchus_milii": "cmilii_gene_ensembl",
}

# Attributes to fetch from BioMart
# These match the attributes requested in the Snakefile rule
BIOMART_ATTRIBUTES = [
    "ensembl_transcript_id",
    "ensembl_transcript_id_version",
    "ensembl_gene_id",
    "external_gene_name",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
]

# ── Helper functions ──────────────────────────────────────────


def escape_xml(text: str) -> str:
    """Escape special characters for XML attributes."""
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )


def build_biomart_query(dataset: str, attributes: list) -> str:
    """
    Build BioMart XML query to fetch all transcripts for a species.

    Parameters
    ----------
    dataset : str
        BioMart dataset name (e.g., "hsapiens_gene_ensembl")
    attributes : list
        List of attribute names to fetch (e.g., "ensembl_transcript_id")

    Returns
    -------
    str
        XML query string
    """
    attr_xml = "\n        ".join(f'<Attribute name="{attr}" />' for attr in attributes)

    query = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
    <Dataset name="{dataset}" interface="default">
        {attr_xml}
    </Dataset>
</Query>"""

    return query


def query_biomart(
    dataset: str, attributes: list, max_retries: int = 3, retry_wait: int = 5
) -> pd.DataFrame:
    """
    Query Ensembl BioMart for all transcripts of a given species.

    Parameters
    ----------
    dataset : str
        BioMart dataset name
    attributes : list
        List of attributes to fetch
    max_retries : int
        Number of retry attempts on network failure
    retry_wait : int
        Seconds to wait between retries (with exponential backoff)

    Returns
    -------
    pd.DataFrame
        Query results with transcript data
    """
    query = build_biomart_query(dataset, attributes)

    for attempt in range(1, max_retries + 1):
        try:
            log.info(f"Querying BioMart {dataset} (attempt {attempt}/{max_retries})")

            response = requests.post(
                BIOMART_URL,
                data={"query": query},
                timeout=300,  # 5-minute timeout for large queries
            )

            if response.status_code == 200:
                # Parse TSV response
                lines = response.text.strip().split("\n")

                if len(lines) < 2:
                    log.warning(
                        f"BioMart returned empty result (only {len(lines)} lines)"
                    )
                    if attempt < max_retries:
                        wait_time = retry_wait * (2 ** (attempt - 1))
                        log.info(f"Retrying in {wait_time} seconds...")
                        time.sleep(wait_time)
                        continue
                    else:
                        raise RuntimeError(
                            "BioMart query returned no data after all retries"
                        )

                # Parse header and data
                header = lines[0].split("\t")
                data = [line.split("\t") for line in lines[1:] if line.strip()]

                if not data:
                    log.warning("BioMart returned only header, no data rows")
                    if attempt < max_retries:
                        wait_time = retry_wait * (2 ** (attempt - 1))
                        log.info(f"Retrying in {wait_time} seconds...")
                        time.sleep(wait_time)
                        continue
                    else:
                        raise RuntimeError(
                            "BioMart query returned no data rows after all retries"
                        )

                df = pd.DataFrame(data, columns=header)
                log.info(
                    f"BioMart query successful: fetched {len(df)} transcript records"
                )
                return df

            else:
                log.warning(f"BioMart returned HTTP {response.status_code}")
                if response.text:
                    log.debug(f"Response body: {response.text[:500]}")

                if attempt < max_retries:
                    wait_time = retry_wait * (2 ** (attempt - 1))
                    log.info(f"Retrying in {wait_time} seconds...")
                    time.sleep(wait_time)
                    continue
                else:
                    raise RuntimeError(f"BioMart returned HTTP {response.status_code}")

        except requests.exceptions.Timeout:
            log.error(f"BioMart query timed out (attempt {attempt}/{max_retries})")
            if attempt < max_retries:
                wait_time = retry_wait * (2 ** (attempt - 1))
                log.info(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
                continue
            else:
                raise RuntimeError("BioMart query timed out after all retries")

        except requests.exceptions.ConnectionError as exc:
            log.error(
                f"Connection error to BioMart (attempt {attempt}/{max_retries}): {exc}"
            )
            if attempt < max_retries:
                wait_time = retry_wait * (2 ** (attempt - 1))
                log.info(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
                continue
            else:
                raise RuntimeError(
                    f"Could not connect to BioMart after {max_retries} attempts"
                )

        except Exception as exc:
            log.error(
                f"Unexpected error during BioMart query (attempt {attempt}/{max_retries}): {exc}"
            )
            if attempt < max_retries:
                wait_time = retry_wait * (2 ** (attempt - 1))
                log.info(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
                continue
            else:
                raise RuntimeError(
                    f"BioMart query failed after {max_retries} attempts: {exc}"
                )


def write_gzipped_tsv(df: pd.DataFrame, output_path: str) -> None:
    """
    Write DataFrame to gzipped TSV file.

    Parameters
    ----------
    df : pd.DataFrame
        Data to write
    output_path : str
        Path to output .tsv.gz file
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write to temporary uncompressed file first
    temp_path = output_path.with_suffix(".tsv")
    df.to_csv(temp_path, sep="\t", index=False)

    # Compress to .gz
    with open(temp_path, "rb") as f_in:
        with gzip.open(output_path, "wb") as f_out:
            f_out.writelines(f_in)

    # Clean up temporary file
    temp_path.unlink()

    log.info(f"Wrote {len(df)} rows to {output_path}")


# ── Main ────────────────────────────────────────────────────────


def main():
    """Main entry point."""
    log.info(f"=== Ensembl BioMart Download ===")
    log.info(f"Species: {species}")
    log.info(f"Release: {release}")
    log.info(f"Output: {output_file}")

    # Validate species
    if species not in SPECIES_TO_DATASET:
        log.error(f"Unknown species: {species}")
        log.error(f"Supported species: {', '.join(sorted(SPECIES_TO_DATASET.keys()))}")
        raise ValueError(f"Unknown species: {species}")

    dataset = SPECIES_TO_DATASET[species]
    log.info(f"Dataset: {dataset}")

    # Query BioMart
    try:
        df = query_biomart(dataset, BIOMART_ATTRIBUTES, max_retries, retry_wait)
    except Exception as exc:
        if on_failure == "warn_continue":
            log.warning(f"BioMart query failed (development mode): {exc}")
            log.warning(f"Creating empty placeholder output for {species}")
            df = pd.DataFrame(columns=BIOMART_ATTRIBUTES)
        else:
            log.error(f"BioMart download failed: {exc}")
            raise

    # Validate results
    if df.empty:
        if on_failure == "warn_continue":
            log.warning(
                "BioMart query returned no data (development mode: continuing with empty output)"
            )
        else:
            log.error("BioMart query returned no data")
            raise RuntimeError("BioMart query returned no data")

    log.info(f"Fetched {len(df)} transcripts for {species}")

    # Check for expected columns
    expected_cols = set(BIOMART_ATTRIBUTES)
    actual_cols = set(df.columns)
    missing_cols = expected_cols - actual_cols

    if missing_cols:
        log.warning(f"Missing columns in BioMart response: {missing_cols}")
        log.warning(f"Available columns: {list(actual_cols)}")
        # Don't fail here; downstream rules may tolerate missing columns

    # Summary statistics
    log.info(f"Column count: {len(df.columns)}")
    log.info(f"Row count: {len(df)}")

    # Sample of data
    log.debug(f"Sample rows:\n{df.head(3).to_string()}")

    # Write output
    try:
        write_gzipped_tsv(df, output_file)
    except Exception as exc:
        log.error(f"Failed to write output file: {exc}")
        raise

    log.info(f"=== Download complete ===")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log.error(f"Fatal error: {exc}", exc_info=True)
        sys.exit(1)
