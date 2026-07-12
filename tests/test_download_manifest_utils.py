import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).parent.parent / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from cache_key_utils import build_cache_key_from_url
from download_manifest_utils import build_download_manifest


def test_build_cache_key_from_url_uses_basename_and_is_stable():
    url = (
        "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/"
        "arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
    )
    key1 = build_cache_key_from_url(url)
    key2 = build_cache_key_from_url(url)

    assert key1 == key2
    assert key1.startswith("url_Arabidopsis_thaliana.TAIR10.dna.toplevel_")


def test_build_download_manifest_includes_ncbi_and_url_backed_rows():
    df = pd.DataFrame(
        {
            "assembly_accession": [
                "GCF_000001405.40",
                pd.NA,
                "Phytozome",
            ],
            "fasta_url": [
                pd.NA,
                "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz",
                pd.NA,
            ],
        }
    )

    manifest = build_download_manifest(df)

    assert "GCF_000001405.40" in manifest["cache_key"].values
    plant_key = build_cache_key_from_url(
        "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz"
    )
    assert plant_key in manifest["cache_key"].values
    assert "Phytozome" not in manifest["cache_key"].values


def test_build_download_manifest_deduplicates_cache_keys():
    plant_url = (
        "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/"
        "zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz"
    )
    key = build_cache_key_from_url(plant_url)
    df = pd.DataFrame(
        {
            "assembly_accession": [pd.NA, pd.NA],
            "fasta_url": [plant_url, plant_url],
        }
    )

    manifest = build_download_manifest(df)
    rows = manifest[manifest["cache_key"] == key]

    assert len(rows) == 1
