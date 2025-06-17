#!/usr/bin/env python3
"""
@brief    Compare GTDB-Tk classifications before and after dRep and flag disagreements.
@param input.pre   gtdbtk.tsv file generated on MetaBAT2 bins (pre-dRep).
@param input.post  gtdbtk.tsv file generated on dRep representatives (post-dRep).
@output table      TSV with MAG-wise taxonomy comparison and agreement flag.
@output summary    Plain-text summary with counts of agreements/disagreements.

This script is executed by Snakemake via the `script:` directive.  All paths are
supplied through the global ``snakemake`` object injected by Snakemake.
"""

from __future__ import annotations

import pandas as pd
import textwrap

# ---------------------------------------------------------------------
# Access Snakemake IO
# ---------------------------------------------------------------------
pre_path: str = snakemake.input["pre"]  # type: ignore
post_path: str = snakemake.input["post"]  # type: ignore
out_table: str = snakemake.output["table"]  # type: ignore
out_summary: str = snakemake.output["summary"]  # type: ignore
sample: str = snakemake.wildcards["sample"]  # type: ignore

# ---------------------------------------------------------------------
# Load tables
# ---------------------------------------------------------------------
pre_df = pd.read_csv(pre_path, sep="\t")
post_df = pd.read_csv(post_path, sep="\t", names=["MAG", "Taxonomy"], header=0)

# ---------------------------------------------------------------------
# Merge + evaluate agreement
# ---------------------------------------------------------------------
merged = pre_df.merge(post_df, on="MAG", how="inner", suffixes=("_pre", "_post"))
merged["agree"] = merged["Taxonomy_pre"] == merged["Taxonomy_post"]
merged.to_csv(out_table, sep="\t", index=False)

# ---------------------------------------------------------------------
# Write summary
# ---------------------------------------------------------------------
_total = len(merged)
_agree = merged["agree"].sum()
_disagree = _total - _agree

with open(out_summary, "w") as fh:
    fh.write(
        textwrap.dedent(
            f"""\
            Sample: {sample}
            Representatives compared: {_total}
            Same taxonomy:            {_agree}
            Different taxonomy:       {_disagree}
            """
        )
    ) 