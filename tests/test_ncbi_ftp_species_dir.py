"""Unit tests for ncbi_assembly_utils.ncbi_ftp_species_dir (FTP path math)."""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.parent / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from ncbi_assembly_utils import NCBI_FTP_ALL_BASE, ncbi_ftp_species_dir


def test_gcf_accession_path():
    assert (
        ncbi_ftp_species_dir("GCF_000001735.4")
        == f"{NCBI_FTP_ALL_BASE}/GCF/000/001/735/"
    )


def test_gca_accession_path():
    assert (
        ncbi_ftp_species_dir("GCA_036512215.1")
        == f"{NCBI_FTP_ALL_BASE}/GCA/036/512/215/"
    )


def test_version_suffix_ignored():
    # Different .N versions share the same species dir.
    assert ncbi_ftp_species_dir("GCF_000001735.9") == ncbi_ftp_species_dir(
        "GCF_000001735.4"
    )
