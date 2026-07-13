#!/usr/bin/env python3
"""Quick test of Gramene API for sample IDs"""

import time

import pandas as pd
import requests

# Load sample of unresolved IDs
df = pd.read_csv("results/biomart_unresolved.tsv", sep="\t")
print(f"Total unresolved: {len(df)}")

# Take first 20 for quick test
sample = df.head(20)


def normalize_gene_id(transcript_id):
    """Strip transcript suffixes"""
    gene_id = transcript_id
    if "_T" in gene_id:
        gene_id = gene_id.split("_T")[0]
    elif gene_id.count("_") > 0:
        parts = gene_id.split("_")
        if parts[-1].isdigit() or (len(parts[-1]) == 2 and parts[-1][0].isdigit()):
            gene_id = "_".join(parts[:-1])
    if "." in gene_id:
        parts = gene_id.split(".")
        gene_id = ".".join(parts[:-1]) if len(parts) > 1 else parts[0]
    return gene_id


resolved = 0
for idx, row in sample.iterrows():
    transcript_id = row["transcript_id"]
    gene_id = normalize_gene_id(transcript_id)

    url = f"https://data.gramene.org/v69/search?q={gene_id}"
    response = requests.get(url, timeout=30)

    if response.status_code == 200:
        data = response.json()
        num_found = data["response"]["numFound"]

        if num_found > 0:
            doc = data["response"]["docs"][0]
            modern_id = doc["id"]
            synonyms = doc.get("synonyms", [])
            resolved += 1
            print(
                f"✓ {transcript_id} → {gene_id} → {modern_id} (synonyms: {synonyms[:2]})"
            )
        else:
            print(f"✗ {transcript_id} → {gene_id} (not found)")
    else:
        print(f"✗ {transcript_id} (HTTP {response.status_code})")

    time.sleep(0.1)

print(f"\nResolved {resolved}/{len(sample)} ({100*resolved/len(sample):.1f}%)")
