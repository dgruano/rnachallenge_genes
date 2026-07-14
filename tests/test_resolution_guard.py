"""Tests for the namespace→assembly presence-verification guard."""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from resolution_guard import check_match_rates


def test_flags_zero_match_species():
    # brachypodium aimed at wrong assembly → 0/142 → flagged
    fails = check_match_rates(matched={}, attempted={"brachypodium_distachyon": 142})
    assert [f[0] for f in fails] == ["brachypodium_distachyon"]
    assert fails[0][1:] == (0, 142, 0.0)


def test_partial_match_passes():
    # maize mix: only one of two namespaces matches per resolver (~52%) → not flagged
    fails = check_match_rates(matched={"zea_mays": 561}, attempted={"zea_mays": 1073})
    assert fails == []


def test_below_min_attempts_ignored():
    # too few IDs to judge, even at 0 match
    fails = check_match_rates(matched={}, attempted={"tiny": 3}, min_attempts=5)
    assert fails == []


def test_threshold_boundary():
    # exactly at min_rate is OK; just under flags
    assert (
        check_match_rates(matched={"x": 2}, attempted={"x": 100}, min_rate=0.02) == []
    )
    fails = check_match_rates(matched={"x": 1}, attempted={"x": 100}, min_rate=0.02)
    assert [f[0] for f in fails] == ["x"]


def test_logs_each_failure():
    msgs = []

    class _L:
        def error(self, m):
            msgs.append(m)

    check_match_rates(matched={}, attempted={"a": 50, "b": 50}, log=_L())
    assert len(msgs) == 2
    assert all("presence-check FAIL" in m for m in msgs)
