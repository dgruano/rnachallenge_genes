import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).parent.parent / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from url_fill_utils import fill_urls_from_table

FILL_COLS = ["assembly_accession", "fasta_url", "gtf_url", "gtf_format"]


def test_fill_urls_from_assembly_name_first():
    df = pd.DataFrame(
        {
            "assembly_name": ["TAIR10"],
            "organism": ["arabidopsis_thaliana"],
            "assembly_accession": [pd.NA],
            "fasta_url": [pd.NA],
            "gtf_url": [pd.NA],
            "gtf_format": [pd.NA],
        }
    )
    url_table = pd.DataFrame(
        {
            "assembly_name": ["TAIR10"],
            "organism": ["arabidopsis_thaliana"],
            "assembly_accession": ["GCF_000001735.4"],
            "fasta_url": ["https://example.org/Arabidopsis_thaliana.TAIR10.fa.gz"],
            "gtf_url": ["https://example.org/Arabidopsis_thaliana.TAIR10.gtf.gz"],
            "gtf_format": ["gtf"],
        }
    )

    filled, count = fill_urls_from_table(df, url_table, fill_cols=FILL_COLS)

    assert count == 1
    assert (
        filled.loc[0, "fasta_url"]
        == "https://example.org/Arabidopsis_thaliana.TAIR10.fa.gz"
    )
    assert filled.loc[0, "assembly_accession"] == "GCF_000001735.4"


def test_fill_urls_falls_back_to_organism():
    df = pd.DataFrame(
        {
            "assembly_name": [pd.NA],
            "organism": ["oryza_sativa"],
            "assembly_accession": [pd.NA],
            "fasta_url": [pd.NA],
            "gtf_url": [pd.NA],
            "gtf_format": [pd.NA],
        }
    )
    url_table = pd.DataFrame(
        {
            "assembly_name": ["IRGSP-1.0"],
            "organism": ["oryza_sativa"],
            "assembly_accession": [pd.NA],
            "fasta_url": ["https://example.org/Oryza_sativa.IRGSP-1.0.fa.gz"],
            "gtf_url": ["https://example.org/Oryza_sativa.IRGSP-1.0.gtf.gz"],
            "gtf_format": ["gtf"],
        }
    )

    filled, count = fill_urls_from_table(
        df,
        url_table,
        fill_cols=FILL_COLS,
        fallback_on_organism=True,
    )

    assert count == 1
    assert (
        filled.loc[0, "fasta_url"] == "https://example.org/Oryza_sativa.IRGSP-1.0.fa.gz"
    )
