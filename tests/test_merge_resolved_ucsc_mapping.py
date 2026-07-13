"""
tests/test_merge_resolved_ucsc_mapping.py

Tests for EXTENDED_UCSC_TO_GCF and apply_ucsc_to_gcf_mapping in ncbi_assembly_utils.
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPT_DIR = Path(__file__).parent.parent / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from ncbi_assembly_utils import EXTENDED_UCSC_TO_GCF, apply_ucsc_to_gcf_mapping


class TestExtendedUCSCMapping:

    def test_ce10_maps_to_gcf(self):
        assert EXTENDED_UCSC_TO_GCF.get("CE10") == "GCF_000002985.6"

    def test_dm6_maps_to_gcf(self):
        assert EXTENDED_UCSC_TO_GCF.get("DM6") == "GCF_000001215.4"

    def test_danrer7_maps_to_gcf(self):
        """danRer7 (zebrafish GRCz10) — only in noncode_v4, must be in extended map."""
        assert EXTENDED_UCSC_TO_GCF.get("DANRER7") == "GCF_000002035.5"

    def test_dm3_maps_to_gcf(self):
        """dm3 (BDGP5) — only in noncode_v4."""
        assert EXTENDED_UCSC_TO_GCF.get("DM3") == "GCF_000001215.3"

    def test_galgal3_maps_to_gcf(self):
        """galGal3 — only in noncode_v4."""
        assert EXTENDED_UCSC_TO_GCF.get("GALGAL3") == "GCF_000002315.5"


class TestApplyUCSCToGCFMapping:

    def _make_frame(self, assembly_values: list[str]) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "transcript_id": [f"TX{i}" for i in range(len(assembly_values))],
                "assembly_accession": assembly_values,
                "db_source": ["noncode_v4"] * len(assembly_values),
            }
        )

    def test_ucsc_name_replaced_with_gcf(self):
        df = self._make_frame(["ce10", "dm6"])
        result = apply_ucsc_to_gcf_mapping(df)
        assert result.loc[0, "assembly_accession"] == "GCF_000002985.6"
        assert result.loc[1, "assembly_accession"] == "GCF_000001215.4"

    def test_already_gcf_passthrough(self):
        df = self._make_frame(["GCF_000001405.40"])
        result = apply_ucsc_to_gcf_mapping(df)
        assert result.loc[0, "assembly_accession"] == "GCF_000001405.40"

    def test_unknown_name_passthrough(self):
        df = self._make_frame(["UnknownAssembly"])
        result = apply_ucsc_to_gcf_mapping(df)
        assert result.loc[0, "assembly_accession"] == "UnknownAssembly"

    def test_danrer7_replaced(self):
        df = self._make_frame(["danRer7"])
        result = apply_ucsc_to_gcf_mapping(df)
        assert result.loc[0, "assembly_accession"] == "GCF_000002035.5"

    def test_empty_frame_passthrough(self):
        df = pd.DataFrame(columns=["transcript_id", "assembly_accession", "db_source"])
        result = apply_ucsc_to_gcf_mapping(df)
        assert result.empty
