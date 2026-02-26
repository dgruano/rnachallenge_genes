"""
utils/logging_utils.py
Shared structured logger for all pipeline scripts.
Each script gets both:
  - A Python logging.Logger writing structured lines to the Snakemake log file
  - Console output (stderr) for interactive runs
"""

import logging
import sys
from pathlib import Path


def get_logger(name: str, log_path: str | Path | None = None) -> logging.Logger:
    """
    Return a configured logger that writes to both stderr and an optional log file.

    Parameters
    ----------
    name     : Logger name (typically the script/module name)
    log_path : Path to the Snakemake log file (snakemake.log[0] or similar)
    """
    fmt = logging.Formatter(
        fmt="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()  # avoid duplicate handlers on re-import

    # Console handler (stderr)
    ch = logging.StreamHandler(sys.stderr)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    # File handler (Snakemake .log file)
    if log_path is not None:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(str(log_path), mode="w", encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)

    return logger
