"""
scripts/download_metadata_table.py
Download metadata tables from Ensembl Plants BioMart
====================================================
Fetches transcript-to-gene mappings with genomic coordinates
for plant species to enable fast local ID resolution.
"""

import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("download_metadata_table", snakemake.log[0])
species = snakemake.wildcards.species
output_file = snakemake.output[0]
cfg = snakemake.config

BIOMART_URL = "https://plants.ensembl.org/biomart/martservice"
MAX_RETRIES = int(cfg.get("max_retries", 3))
RETRY_WAIT = int(cfg.get("retry_wait_seconds", 5))

# Map species to BioMart dataset names
SPECIES_TO_DATASET = {
    "arabidopsis_thaliana": "athaliana_eg_gene",
    "oryza_sativa": "osativa_eg_gene",
    "zea_mays": "zmays_eg_gene",
    "solanum_lycopersicum": "slycopersicum_eg_gene",
    "glycine_max": "gmax_eg_gene",
    "triticum_aestivum": "taestivum_eg_gene",
    "brachypodium_distachyon": "bdistachyon_eg_gene",
    "sorghum_bicolor": "sbicolor_eg_gene",
    "vitis_vinifera": "vvinifera_eg_gene",
}

def build_biomart_query(dataset: str) -> str:
    """Build BioMart XML query for transcript metadata."""
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="plants_mart" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
    <Dataset name="{dataset}" interface="default">
        <Attribute name="ensembl_transcript_id" />
        <Attribute name="ensembl_gene_id" />
        <Attribute name="external_gene_name" />
        <Attribute name="chromosome_name" />
        <Attribute name="transcript_start" />
        <Attribute name="transcript_end" />
        <Attribute name="strand" />
    </Dataset>
</Query>"""

def download_biomart_table(dataset: str) -> pd.DataFrame:
    """Download metadata table from BioMart with retries."""
    query = build_biomart_query(dataset)
    
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            log.info(f"Querying BioMart for {dataset} (attempt {attempt}/{MAX_RETRIES})")
            response = requests.post(
                BIOMART_URL,
                data={"query": query},
                timeout=300,  # 5 minutes
            )
            
            if response.status_code == 200:
                # Parse TSV response
                lines = response.text.strip().split("\n")
                if len(lines) < 2:
                    log.error(f"Empty response from BioMart for {dataset}")
                    return pd.DataFrame()
                
                # Header is first line
                header = lines[0].split("\t")
                data = [line.split("\t") for line in lines[1:]]
                
                df = pd.DataFrame(data, columns=header)
                log.info(f"Downloaded {len(df)} records from BioMart")
                return df
            else:
                log.warning(f"BioMart request failed with status {response.status_code}")
                
        except Exception as exc:
            log.warning(f"BioMart request attempt {attempt} failed: {exc}")
        
        if attempt < MAX_RETRIES:
            time.sleep(RETRY_WAIT * attempt)
    
    return pd.DataFrame()

def normalize_strand(value: str) -> str:
    """Normalize strand to +/-/."""
    if value in ("1", "+1", "+"):
        return "+"
    elif value in ("-1", "-"):
        return "-"
    else:
        return "."

# ── Main ─────────────────────────────────────────────────────
log.info(f"Downloading metadata table for {species}")

if species not in SPECIES_TO_DATASET:
    log.error(f"Unknown species {species}. Supported: {', '.join(SPECIES_TO_DATASET.keys())}")
    sys.exit(1)

dataset = SPECIES_TO_DATASET[species]
df = download_biomart_table(dataset)

if df.empty:
    log.error(f"Failed to download metadata for {species}")
    sys.exit(1)

# Rename columns to match expected format
column_mapping = {
    "Gene stable ID": "gene_id",
    "Transcript stable ID": "transcript_id",
    "Gene name": "gene_symbol",
    "Chromosome/scaffold name": "chrom",
    "Transcript start (bp)": "start",
    "Transcript end (bp)": "end",
    "Strand": "strand",
}

df = df.rename(columns=column_mapping)

# Fill missing gene_symbol with gene_id
if "gene_symbol" in df.columns:
    df["gene_symbol"] = df["gene_symbol"].fillna(df["gene_id"])

# Ensure all required columns exist
required_cols = ["transcript_id", "gene_id", "gene_symbol", "chrom", "start", "end", "strand"]
missing = set(required_cols) - set(df.columns)
if missing:
    log.error(f"Missing columns after renaming: {', '.join(missing)}")
    log.error(f"Available columns: {', '.join(df.columns)}")
    sys.exit(1)

# Add assembly_accession (we'll use the species name for now, or could query separately)
df["assembly_accession"] = "EnsemblPlants"  # Placeholder - ideally we'd get the actual assembly name

# Normalize strand
df["strand"] = df["strand"].astype(str).apply(normalize_strand)

# Select and reorder columns
output_cols = ["transcript_id", "gene_id", "gene_symbol", "assembly_accession", "chrom", "start", "end", "strand"]
df = df[output_cols]

# Remove rows with missing critical info
df = df.dropna(subset=["transcript_id", "gene_id"])

# Save
output_path = Path(output_file)
output_path.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(output_file, sep="\t", index=False)

log.info(f"Wrote {len(df)} records to {output_file}")
log.info("download_metadata_table complete.")
