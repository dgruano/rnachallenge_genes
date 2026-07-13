"""
Add tool-source annotations to unresolved transcript tables.

This joins an input TSV containing transcript_id against tool_source_map.tsv and
appends the most useful provenance columns for triage.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

log = get_logger("annotate_with_tools", snakemake.log[0])

in_records = snakemake.input.records
in_tool_map = snakemake.input.tool_map
out_annotated = snakemake.output.annotated

ANNOTATION_COLUMNS = ["transcript_id", "tools", "n_tools", "match_strategy"]


def main() -> None:
    records_df = pd.read_csv(in_records, sep="\t")
    tool_map_df = pd.read_csv(in_tool_map, sep="\t")

    missing_cols = [col for col in ANNOTATION_COLUMNS if col not in tool_map_df.columns]
    if missing_cols:
        raise ValueError(
            f"tool map is missing expected columns: {', '.join(missing_cols)}"
        )

    annotated_df = records_df.merge(
        tool_map_df[ANNOTATION_COLUMNS],
        on="transcript_id",
        how="left",
        validate="many_to_one",
    )
    annotated_df.to_csv(out_annotated, sep="\t", index=False)

    log.info("Annotated %s rows from %s", len(annotated_df), Path(in_records).name)
    log.info("Wrote %s", out_annotated)


main()
