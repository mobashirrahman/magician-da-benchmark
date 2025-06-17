#!/usr/bin/env python3
"""
@brief    Build a wide contig-level mean-depth matrix from per-sample coverage depth tables.
@param input   coverage/*.txt files produced by jgi_summarize_bam_contig_depths (one per sample).
@output matrix coverage/mean_depth.tsv – rows = contig IDs, columns = samples, values = mean depth.

The script is executed via Snakemake's `script:` directive.  All paths are supplied
through the global ``snakemake`` object injected by Snakemake at runtime.
"""

import pathlib
from typing import Sequence

import pandas as pd

# ------------------------------------------------------------------
# Access Snakemake-provided variables
# ------------------------------------------------------------------
INPUT_FILES: Sequence[str] = list(snakemake.input)  # type: ignore  # noqa: F821
OUT_MATRIX: str = snakemake.output["matrix"]  # type: ignore  # noqa: F821

# ------------------------------------------------------------------
# Constants – header spellings accepted for the ID & depth columns
# ------------------------------------------------------------------
ID_CANDIDATES = (
    "contigName",
    "contig",
    "name",
    "contig_id",
)
DEPTH_CANDIDATES = (
    "avg_cov",
    "avgDepth",
    "totalAvgDepth",
    "mean",
    "coverage",
)

# ------------------------------------------------------------------
# Parse each per-sample depth table into a Series → concatenate
# ------------------------------------------------------------------
series_list = []
for file_path in INPUT_FILES:
    sample = pathlib.Path(file_path).stem
    df = pd.read_csv(file_path, sep="\t")

    # ---- Identify contig ID column
    try:
        id_col = next(
            c for c in df.columns if c in ID_CANDIDATES or c.lower().startswith("contig")
        )
    except StopIteration as exc:
        raise ValueError(
            f"No contig column found in {file_path}. Headers: {list(df.columns)}"
        ) from exc

    # ---- Identify depth/coverage column
    try:
        depth_col = next(
            c
            for c in df.columns
            if c in DEPTH_CANDIDATES or "depth" in c.lower() or "cov" in c.lower()
        )
    except StopIteration as exc:
        raise ValueError(
            f"No depth/coverage column found in {file_path}. Headers: {list(df.columns)}"
        ) from exc

    series_list.append(df.set_index(id_col)[depth_col].rename(sample))

# ------------------------------------------------------------------
# Build wide matrix (fill missing contig-sample pairs with 0) & write
# ------------------------------------------------------------------
wide_df = pd.concat(series_list, axis=1).fillna(0)
wide_df.to_csv(OUT_MATRIX, sep="\t") 