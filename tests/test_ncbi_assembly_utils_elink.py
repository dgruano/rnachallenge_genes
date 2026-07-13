"""
tests/test_ncbi_assembly_utils_elink.py

Unit tests for map_genomic_to_assembly_elink() added to ncbi_assembly_utils.
"""

import sys
from pathlib import Path
from unittest.mock import patch

import pytest

SCRIPT_DIR = Path(__file__).parent.parent / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from ncbi_assembly_utils import map_genomic_to_assembly_elink

# ---------------------------------------------------------------------------
# Helpers — minimal Biopython-style mock objects
# ---------------------------------------------------------------------------


def _make_search_result(uid: str) -> dict:
    return {"IdList": [uid], "Count": "1", "RetMax": "1", "RetStart": "0"}


def _make_empty_search_result() -> dict:
    return {"IdList": [], "Count": "0", "RetMax": "0", "RetStart": "0"}


def _make_elink_result(from_uid: str, to_uids: list[str]) -> list:
    """One LinkSet entry for a single nuccore→assembly link."""
    links = [{"Id": uid} for uid in to_uids]
    return [
        {
            "DbFrom": "nuccore",
            "IdList": [from_uid],
            "LinkSetDb": [
                {"DbTo": "assembly", "LinkName": "nuccore_assembly", "Link": links}
            ],
        }
    ]


def _make_empty_elink_result(from_uid: str) -> list:
    return [{"DbFrom": "nuccore", "IdList": [from_uid], "LinkSetDb": []}]


class _DocSummary(dict):
    """Minimal mock for Biopython DocumentSummary (dict with .attributes)."""

    def __init__(self, data: dict, uid: str):
        super().__init__(data)
        self.attributes = {"uid": uid}


def _make_esummary_result(
    uid: str, accession: str, organism: str = "Homo sapiens"
) -> dict:
    return {
        "DocumentSummarySet": {
            "DocumentSummary": [
                _DocSummary({"AssemblyAccession": accession, "Organism": organism}, uid)
            ]
        }
    }


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestMapGenomicToAssemblyElink:

    def test_single_accession_mapped(self):
        """NC_000001.11 → GCF_000001405.40 via esearch/elink/esummary chain."""
        with (
            patch("ncbi_assembly_utils.Entrez.esearch") as mock_search,
            patch("ncbi_assembly_utils.Entrez.elink") as mock_elink,
            patch("ncbi_assembly_utils.Entrez.esummary") as mock_summary,
            patch("ncbi_assembly_utils.Entrez.read") as mock_read,
            patch("ncbi_assembly_utils.time.sleep"),
        ):
            # esearch returns nuccore UID for NC_000001.11
            # elink returns assembly UID 999
            # esummary returns GCF_000001405.40 for assembly UID 999
            mock_read.side_effect = [
                _make_search_result("111111"),  # esearch
                _make_elink_result("111111", ["999"]),  # elink
                _make_esummary_result("999", "GCF_000001405.40"),  # esummary
            ]

            result = map_genomic_to_assembly_elink(["NC_000001.11"])

        assert result == {"NC_000001.11": "GCF_000001405.40"}

    def test_accession_not_found_in_esearch(self):
        """Returns None for accessions with no nuccore hit."""
        with (
            patch("ncbi_assembly_utils.Entrez.esearch"),
            patch("ncbi_assembly_utils.Entrez.elink"),
            patch("ncbi_assembly_utils.Entrez.esummary"),
            patch("ncbi_assembly_utils.Entrez.read") as mock_read,
            patch("ncbi_assembly_utils.time.sleep"),
        ):
            mock_read.side_effect = [_make_empty_search_result()]

            result = map_genomic_to_assembly_elink(["FAKE_000.1"])

        assert result == {"FAKE_000.1": None}

    def test_no_assembly_link(self):
        """Returns None when elink finds no assembly for the nuccore record."""
        with (
            patch("ncbi_assembly_utils.Entrez.esearch"),
            patch("ncbi_assembly_utils.Entrez.elink"),
            patch("ncbi_assembly_utils.Entrez.esummary"),
            patch("ncbi_assembly_utils.Entrez.read") as mock_read,
            patch("ncbi_assembly_utils.time.sleep"),
        ):
            mock_read.side_effect = [
                _make_search_result("111111"),
                _make_empty_elink_result("111111"),
            ]

            result = map_genomic_to_assembly_elink(["NT_033779.5"])

        assert result == {"NT_033779.5": None}

    def test_multiple_accessions_batched(self):
        """Two accessions map independently; only 1 elink batch + 1 esummary batch needed."""
        with (
            patch("ncbi_assembly_utils.Entrez.esearch"),
            patch("ncbi_assembly_utils.Entrez.elink"),
            patch("ncbi_assembly_utils.Entrez.esummary"),
            patch("ncbi_assembly_utils.Entrez.read") as mock_read,
            patch("ncbi_assembly_utils.time.sleep"),
        ):
            # Two esearches (one per accession), then batched elink + esummary
            combined_elink = _make_elink_result("111", ["901"]) + _make_elink_result(
                "222", ["902"]
            )
            mock_read.side_effect = [
                _make_search_result("111"),  # esearch NC_000001.11
                _make_search_result("222"),  # esearch NC_000002.12
                combined_elink,  # elink batch
                {
                    "DocumentSummarySet": {
                        "DocumentSummary": [
                            _DocSummary(
                                {
                                    "AssemblyAccession": "GCF_000001405.40",
                                    "Organism": "Homo sapiens",
                                },
                                "901",
                            ),
                            _DocSummary(
                                {
                                    "AssemblyAccession": "GCF_000001405.40",
                                    "Organism": "Homo sapiens",
                                },
                                "902",
                            ),
                        ]
                    }
                },
            ]

            result = map_genomic_to_assembly_elink(["NC_000001.11", "NC_000002.12"])

        assert result["NC_000001.11"] == "GCF_000001405.40"
        assert result["NC_000002.12"] == "GCF_000001405.40"

    def test_empty_input(self):
        """Empty input returns empty dict without making any API calls."""
        with patch("ncbi_assembly_utils.Entrez.esearch") as mock_search:
            result = map_genomic_to_assembly_elink([])
        mock_search.assert_not_called()
        assert result == {}

    def test_deduplicates_assembly_uids_before_esummary(self):
        """Multiple chromosomes pointing to same assembly use only one esummary call."""
        with (
            patch("ncbi_assembly_utils.Entrez.esearch"),
            patch("ncbi_assembly_utils.Entrez.elink"),
            patch("ncbi_assembly_utils.Entrez.esummary") as mock_esummary_call,
            patch("ncbi_assembly_utils.Entrez.read") as mock_read,
            patch("ncbi_assembly_utils.time.sleep"),
        ):
            same_asm_uid = "999"
            mock_read.side_effect = [
                _make_search_result("111"),
                _make_search_result("222"),
                _make_elink_result("111", [same_asm_uid])
                + _make_elink_result("222", [same_asm_uid]),
                _make_esummary_result(same_asm_uid, "GCF_000001405.40"),
            ]

            result = map_genomic_to_assembly_elink(["NC_000001.11", "NC_000002.12"])

        # esummary called once (deduplicated UID)
        assert mock_esummary_call.call_count == 1
        assert result["NC_000001.11"] == "GCF_000001405.40"
        assert result["NC_000002.12"] == "GCF_000001405.40"
