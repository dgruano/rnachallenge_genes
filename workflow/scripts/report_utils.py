"""Shared helpers for pipeline HTML reports."""

from __future__ import annotations

from html import escape

import pandas as pd


MAX_INLINE_ROWS = 25

FAILURE_ACTIONS = {
    "assembly_not_cached": (
        "Assembly not cached",
        "Re-run the assembly download stage or make sure the missing assembly is available locally.",
    ),
    "chrom_not_found": (
        "Chromosome not found",
        "Check chromosome naming and accession remapping before re-running extraction.",
    ),
    "missing_coordinates": (
        "Missing coordinates",
        "Inspect the upstream resolver output for rows without start/end coordinates.",
    ),
    "sequence_error": (
        "Sequence extraction error",
        "Inspect the failing row and the extraction log for a concrete parser or indexing issue.",
    ),
    "pattern_unmatched": (
        "Pattern unmatched",
        "Extend the ID classifier or add a source-specific parser if these IDs should be supported.",
    ),
}


def pct(num: int, denom: int) -> str:
    return f"{100 * num / denom:.0f}%" if denom else "—"


def stat_card(label: str, value: object, color: str = "#3b82f6") -> str:
    return f"""<div class="stat-card">
      <div class="stat-value" style="color:{color}">{value}</div>
      <div class="stat-label">{escape(str(label))}</div></div>"""


def df_to_html_table(
    df: pd.DataFrame,
    *,
    max_rows: int = MAX_INLINE_ROWS,
    full_href: str | None = None,
    full_label: str = "Full table",
) -> str:
    if df.empty:
        return "<p><em>No data.</em></p>"

    note_parts: list[str] = []
    if len(df) > max_rows:
        df = df.head(max_rows)
        note_parts.append(f"Showing first {max_rows} rows.")

    if full_href:
        note_parts.append(f"<a class=\"file-link\" href=\"{escape(full_href)}\">{escape(full_label)}</a>")

    footer = f"<p class=\"section-note\"><em>{' '.join(note_parts)}</em></p>" if note_parts else ""
    return df.to_html(index=False, border=0, classes="data-table") + footer


def report_intro_html(*, report_title: str, summary: str, related_href: str, related_label: str) -> str:
    return (
        '<div class="report-intro">'
        f"<strong>{escape(report_title)}</strong>"
        f"<p>{escape(summary)}</p>"
        f'<p class="section-note">Open the <a class="file-link" href="{escape(related_href)}">{escape(related_label)}</a> for the complementary view.</p>'
        "</div>"
    )


def next_actions_html(
    *,
    unresolved_count: int = 0,
    unclassified_count: int = 0,
    classified_unresolved_count: int = 0,
    ambiguous_count: int,
    unknown_prefix_count: int,
    extraction_failures: int = 0,
    chrom_failures: int = 0,
    assembly_failures: int = 0,
) -> str:
    actions: list[str] = []

    if unknown_prefix_count:
        actions.append(
            f"{unknown_prefix_count} unknown Ensembl prefixes were detected. Update ensembl_species_overrides and rerun classification."
        )
    if unclassified_count:
        actions.append(
            f"{unclassified_count} IDs were not classified. Inspect pattern_unmatched.tsv and expand parse patterns before rerunning."
        )
    if classified_unresolved_count:
        actions.append(
            f"{classified_unresolved_count} IDs were classified but not resolved. Inspect matched_not_found.tsv to target resolver-specific gaps."
        )
    if unresolved_count and not (unclassified_count or classified_unresolved_count):
        actions.append(
            f"{unresolved_count} unresolved IDs remain. Inspect unresolved.tsv to see which sources or patterns still need support."
        )
    if ambiguous_count:
        actions.append(
            f"{ambiguous_count} IDs mapped to multiple candidates. Review ambiguous.tsv if a deterministic choice matters."
        )
    if extraction_failures:
        if assembly_failures >= chrom_failures:
            actions.append(
                f"{extraction_failures} extraction failures were recorded, mostly from missing assemblies. Re-run the download stage for the affected accessions."
            )
        else:
            actions.append(
                f"{extraction_failures} extraction failures were recorded, mostly from chromosome lookup issues. Check chromosome naming and remapping before rerunning extraction."
            )

    if not actions:
        actions.append("No blocking issues were detected. The run is ready for downstream use.")

    items = "".join(f"<li>{escape(action)}</li>" for action in actions)
    return f"<div class=\"action-box\"><strong>What to do next</strong><ul>{items}</ul></div>"


def failure_summary_rows(df_failed: pd.DataFrame) -> tuple[str, dict[str, int]]:
    if df_failed.empty or "fail_reason" not in df_failed.columns:
        return "<p><em>No failure data.</em></p>", {}

    fail_counts = df_failed["fail_reason"].value_counts().to_dict()
    rows = []
    for reason, count in fail_counts.items():
        label, action = FAILURE_ACTIONS.get(
            reason,
            (
                reason.replace("_", " ").title(),
                "Inspect the corresponding record and pipeline log for details.",
            ),
        )
        rows.append({"Failure reason": label, "Count": count, "Action": action})
    return df_to_html_table(pd.DataFrame(rows), max_rows=MAX_INLINE_ROWS), fail_counts