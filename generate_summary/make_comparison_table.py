import pathlib

from argparse import ArgumentParser

import pandas as pd


def create_comparison_table(base_table: pathlib.Path) -> pd.DataFrame:
    """Generate an overview table comparing bins to their closest reference genomes.
    Arguments:
        base_table: Path to the Excel file summarizing various QC parameters.
    Returns:
        A table directly comparing bins to closest reference genomes.
    """
    # Read in data
    bb_stats = pd.read_excel(base_table, sheet_name="BB_stats", index_col=0)
    checkm_stats = pd.read_excel(base_table, sheet_name="CheckM", index_col=0)
    drep_stats = pd.read_excel(base_table, sheet_name="dRep", index_col=0)

    summary_table = bb_stats[bb_stats["genome_type"] == "synthetic_MAG"][['bin_name']]

    # closest organism: SELECT drep_stats.query AS closest_genome, summary_table.bin_name, MAX(drep_stats.ani)
    # FROM summary table
    # LEFT JOIN drep_stats ON summary_table.bin_name = drep_stats.reference
    # GROUP BY summary_table.bin_name
    # closest organism: find highest ANI for each bin in drep_stats, filter down to only these
    # then left join (leaving all entries from the summary table intact) summary table on filtered drep_stats
    summary_table = pd.merge(summary_table, drep_stats[drep_stats.groupby("reference")['ani'].transform(max)
                                                       == drep_stats['ani']],
                            left_on='bin_name', right_on='reference', how='left')[["bin_name", "query", "ref_coverage",
                                                                                   "query_coverage", "ani"]]
    summary_table = summary_table.rename(columns={"query": "closest_genome", "ref_coverage": "bin_coverage",
                                                  "query_coverage": "source_coverage"})
    # scaffold etc. difference
    # create temp columns - one with scaffold stats, one with source genome stats
    summary_bb = pd.merge(summary_table, bb_stats, how="left", on="bin_name")
    summary_bb = pd.merge(summary_bb, bb_stats, left_on="closest_genome", right_on="bin_name", how="left",
                          suffixes=('_bin', '_ref'))
    summary_bb["scaffold_difference"] = summary_bb['n_scaffolds_bin'] - summary_bb["n_scaffolds_ref"]
    summary_bb["contig_difference"] = summary_bb['n_contigs_bin'] - summary_bb['n_contigs_ref']
    summary_bb["length_difference"] = summary_bb["scaf_bp_bin"] - summary_bb["scaf_bp_ref"]
    summary_bb["gc_difference"] = summary_bb["gc_avg_bin"] - summary_bb["gc_avg_ref"]

    summary_checkm = pd.merge(pd.merge(summary_table, checkm_stats, left_on="bin_name", right_on="Bin Id", how="left"),
                              checkm_stats, how="left", left_on="closest_genome", right_on="Bin Id",
                              suffixes=("_bin", "_ref"))
    summary_checkm["completeness_difference"] = summary_checkm["Completeness_bin"] - summary_checkm["Completeness_ref"]
    summary_table = summary_bb[["bin_name_bin", "closest_genome", "ani", "bin_coverage", "source_coverage",
                                "scaffold_difference", "contig_difference", "length_difference", "gc_difference"]]
    summary_table = summary_table.rename(columns={"bin_name_bin": "bin_name"})
    summary_table = pd.merge(summary_table, summary_checkm[["Bin Id_bin",
                                                            "completeness_difference"]],
                             how="left", left_on="bin_name", right_on="Bin Id_bin").drop("Bin Id_bin", axis=1)
    return summary_table



if __name__ == "__main__":
    parser = ArgumentParser(description="Generate a bin-focused overview from the general MAGICIAN summary file")
    parser.add_argument("infile", action="store", help="Path to summary .xlsx file")
    parser.add_argument("-o", "--outfile", action="store", help="Name of output file (default: bin_summary.xlsx)",
                        default="bin_summary.xlsx")
    args = parser.parse_args()

    infile = pathlib.Path(args.infile).resolve()
    outfile = pathlib.Path(args.outfile).resolve()

    summary = create_comparison_table(infile)

    explanations = pd.DataFrame.from_dict({"bin_name": ["Name of MetaBAT-assigned bin"],
                                           "closest_genome": ["Source genome with highest ANI for bin found by dRep"],
                                           "bin_coverage": ["Percent of bin covered by closest genome"],
                                           "source_coverage": ["Percent of closest genome covered by bin"],
                                           "scaffold_difference":
                                               ["Difference of the amount of scaffolds in the bin and the closest genome (bin - reference)"],
                                           "contig_difference":
                                               ["Difference of the amount of contigs in the bin and the closest genome (bin - reference)"],
                                           "length_difference":
                                               ["Length difference in base pairs between bin and closest genome (bin - reference)"],
                                           "gc_difference":
                                               ["Difference in GC content between bin and closest genome (bin - reference)"],
                                           "completeness_difference":
                                               ["Difference in completeness between bin and closest genome as calculated by CheckM (bin - reference). Note that completeness in bin and reference may be calculated with different marker genes; when in doubt, check the CheckM summary."]
                                           }, orient='index')

    with pd.ExcelWriter(outfile) as outfile_writer:
        summary.to_excel(outfile_writer, sheet_name="summary")
        explanations.to_excel(outfile_writer, sheet_name="explanations")
