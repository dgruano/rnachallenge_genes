"""
scripts/resolution_guard.py
Presence-verification guard for namespace-keyed plant resolution.
================================================================

We resolve plant transcripts by namespace → source (GTF/GFF3). The correctness
guarantee is not that the namespace→assembly mapping is *trusted*, but that it is
*verified*: if a species that has input IDs and a configured annotation matches
≈0 of those IDs, the source is almost certainly the wrong assembly (e.g. GRMZM
aimed at AGPv4, where 0/100 GRMZM cores exist). Turn that silent zero-coordinate
bug into a loud, named failure.

Pure function so it is unit-testable without Snakemake.
"""

from typing import Dict, List, Tuple


def check_match_rates(
    matched: Dict[str, int],
    attempted: Dict[str, int],
    min_rate: float = 0.02,
    min_attempts: int = 5,
    log=None,
) -> List[Tuple[str, int, int, float]]:
    """Return the species whose match rate is suspiciously low.

    matched / attempted are species -> count. A species is flagged when it had at
    least ``min_attempts`` IDs routed to a configured source but matched a fraction
    below ``min_rate`` (the signature of a wrong-assembly source). Species below
    ``min_attempts`` are ignored (too little signal to judge).

    Returns a list of (species, n_matched, n_attempted, rate). Empty == all good.
    The caller decides whether to exit non-zero (strict mode).
    """
    failures: List[Tuple[str, int, int, float]] = []
    for species, n_att in sorted(attempted.items()):
        if n_att < min_attempts:
            continue
        n_match = matched.get(species, 0)
        rate = n_match / n_att if n_att else 0.0
        if rate < min_rate:
            failures.append((species, n_match, n_att, rate))

    if log is not None:
        for species, n_match, n_att, rate in failures:
            log.error(
                f"presence-check FAIL: {species} matched {n_match}/{n_att} "
                f"({rate:.1%}) < {min_rate:.0%} — wrong assembly/source? "
                f"Verify with: zgrep -c -Ff <sample_ids> <that_source>.gff3.gz"
            )
    return failures


def _demo() -> None:
    """Runnable self-check: assert the guard flags a zero-match, spares a partial."""
    log_msgs = []

    class _L:
        def error(self, m):
            log_msgs.append(m)

    attempted = {
        "brachypodium_distachyon": 142,  # wrong assembly → 0 match, must flag
        "zea_mays": 1073,  # maize mix: one namespace matches (~52%), must NOT flag
        "tiny": 3,  # below min_attempts, ignore even at 0 match
    }
    matched = {"zea_mays": 561, "tiny": 0}

    fails = check_match_rates(attempted=attempted, matched=matched, log=_L())
    failed_species = {f[0] for f in fails}
    assert failed_species == {"brachypodium_distachyon"}, failed_species
    assert len(log_msgs) == 1 and "brachypodium_distachyon" in log_msgs[0]
    # partial maize passes; sub-min tiny ignored
    print("resolution_guard self-check OK:", fails)


if __name__ == "__main__":
    _demo()
