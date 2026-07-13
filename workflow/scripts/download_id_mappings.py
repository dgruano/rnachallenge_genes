"""
scripts/download_id_mappings.py
Download ID mapping files for legacy plant nomenclatures
=========================================================
Rice: MSU → IRGSP/RAP-DB
Maize: GRMZM → Zm (modern B73)
"""

import gzip
import sys
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("download_id_mappings", snakemake.log[0])
output_dir = Path(snakemake.output[0]).parent
cfg = snakemake.config

output_dir.mkdir(parents=True, exist_ok=True)

# ── Rice MSU → RAP-DB mapping ────────────────────────────────
log.info("Downloading Rice MSU → RAP-DB ID mapping...")

# RAP-DB provides MSU-IRGSP mapping
# https://rapdb.dna.affrc.go.jp/download/irgsp1.html
RICE_MAPPING_URL = "https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_annotation_2021-11-11.tsv.gz"

try:
    log.info(f"Fetching {RICE_MAPPING_URL}")
    response = requests.get(RICE_MAPPING_URL, timeout=60)

    if response.status_code == 200:
        # Decompress and parse
        import io

        content = gzip.decompress(response.content).decode("utf-8")

        # Parse TSV
        lines = content.strip().split("\n")
        header = lines[0].split("\t")
        data = [line.split("\t") for line in lines[1:]]

        df = pd.DataFrame(data, columns=header)

        # Save relevant columns: transcript ID, gene ID, MSU locus
        # Columns vary, but typically include: Locus_ID, RAP-DB Gene ID, MSU Gene ID
        rice_output = output_dir / "rice_msu_irgsp_mapping.tsv"

        # Try to extract MSU mapping columns
        if "MSU" in "\t".join(header):
            log.info(f"Found MSU columns in RAP-DB annotation")
            # Save full file for manual inspection
            df.to_csv(rice_output, sep="\t", index=False)
            log.info(f"Saved Rice mapping to {rice_output} ({len(df)} rows)")
        else:
            log.warning("No MSU column found in RAP-DB file - checking alternatives")
            # Save anyway for inspection
            df.to_csv(rice_output, sep="\t", index=False)
    else:
        log.error(f"Failed to download Rice mapping: HTTP {response.status_code}")

except Exception as exc:
    log.error(f"Rice mapping download failed: {exc}")
    # Try alternative source: Ensembl Plants FTP
    log.info("Trying alternative: Ensembl Plants BioMart export...")
    alt_url = "https://plants.ensembl.org/biomart/martservice"

    # Query for all rice transcripts with external references
    query = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="plants_mart" formatter="TSV" header="1" uniqueRows="1">
    <Dataset name="osativa_eg_gene" interface="default">
        <Attribute name="ensembl_transcript_id" />
        <Attribute name="ensembl_gene_id" />
        <Attribute name="external_gene_name" />
        <Attribute name="external_transcript_name" />
        <Attribute name="description" />
    </Dataset>
</Query>"""

    try:
        response = requests.post(alt_url, data={"query": query}, timeout=300)
        if response.status_code == 200:
            rice_output = output_dir / "rice_ensembl_export.tsv"
            rice_output.write_text(response.text)
            log.info(f"Downloaded Ensembl Plants rice export to {rice_output}")
    except Exception as e2:
        log.error(f"Alternative rice download also failed: {e2}")

# ── Maize GRMZM → Zm mapping ──────────────────────────────────
log.info("Downloading Maize GRMZM → Zm ID mapping...")

# MaizeGDB provides B73 v3 → v4/v5 mappings
# Try Ensembl Plants FTP for cross-references
MAIZE_XREF_URL = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/tsv/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.110.entrezgene.tsv.gz"

try:
    log.info(f"Fetching {MAIZE_XREF_URL}")
    response = requests.get(MAIZE_XREF_URL, timeout=60)

    if response.status_code == 200:
        content = gzip.decompress(response.content).decode("utf-8")
        maize_output = output_dir / "maize_entrez_xref.tsv"
        maize_output.write_text(content)
        log.info(f"Saved Maize EntrezGene xrefs to {maize_output}")

        # Parse to see if we have useful mappings
        df = pd.read_csv(io.StringIO(content), sep="\t")
        log.info(f"Maize xref columns: {', '.join(df.columns)}")
        log.info(f"Maize xref rows: {len(df)}")
    else:
        log.warning(f"Maize xref download failed: HTTP {response.status_code}")

except Exception as exc:
    log.error(f"Maize mapping download failed: {exc}")
    # Try BioMart export
    log.info("Trying alternative: Ensembl Plants BioMart for maize...")
    alt_url = "https://plants.ensembl.org/biomart/martservice"

    query = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="plants_mart" formatter="TSV" header="1" uniqueRows="1">
    <Dataset name="zmays_eg_gene" interface="default">
        <Attribute name="ensembl_transcript_id" />
        <Attribute name="ensembl_gene_id" />
        <Attribute name="external_gene_name" />
        <Attribute name="external_transcript_name" />
        <Attribute name="description" />
    </Dataset>
</Query>"""

    try:
        response = requests.post(alt_url, data={"query": query}, timeout=300)
        if response.status_code == 200:
            maize_output = output_dir / "maize_ensembl_export.tsv"
            maize_output.write_text(response.text)
            log.info(f"Downloaded Ensembl Plants maize export to {maize_output}")
    except Exception as e2:
        log.error(f"Alternative maize download also failed: {e2}")

log.info("ID mapping download complete. Check resources/id_mappings/ for files.")
log.info("Note: These may require manual parsing to extract GRMZM/MSU mappings.")
