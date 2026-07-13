"""Tests for the static NC_→GCF exception map (ROI strategy #5).

Legacy RefSeq chromosome accessions (e.g. rice Build 4.0) carry no
``Assembly:`` dbxref and their nuccore→assembly elink is empty, so the
efetch-based mapping in ``resolve_ncbi_assembly_accessions`` returns None.
The static exception map recovers them without any network call.
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from ncbi_assembly_utils import NC_TO_ASSEMBLY_EXCEPTIONS, assembly_from_exceptions


class TestNCExceptionMap:
    def test_rice_chr1_maps_to_build_4_0(self):
        # NC_008394.4 = Oryza sativa Japonica chromosome 1, Build 4.0
        assert NC_TO_ASSEMBLY_EXCEPTIONS["NC_008394.4"] == "GCF_000005425.2"

    def test_rice_chr12_maps_to_build_4_0(self):
        assert NC_TO_ASSEMBLY_EXCEPTIONS["NC_008405.2"] == "GCF_000005425.2"

    def test_all_twelve_rice_chromosomes_present(self):
        chrom_accs = [
            f"NC_0083{n:02d}.{'4' if n == 94 else '2'}" for n in range(94, 100)
        ]
        chrom_accs += [f"NC_0084{n:02d}.2" for n in range(0, 6)]
        for acc in chrom_accs:
            assert acc in NC_TO_ASSEMBLY_EXCEPTIONS, acc

    def test_helper_returns_mapping_for_known_acc(self):
        assert assembly_from_exceptions("NC_008394.4") == "GCF_000005425.2"

    def test_helper_returns_none_for_unknown_acc(self):
        assert assembly_from_exceptions("NC_000001.11") is None
