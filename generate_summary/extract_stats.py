import pathlib

from argparse import ArgumentParser

import numpy as np
import pandas as pd


# TODO: style pass!
# TODO: fix for large dRep clusters!
def get_bb_stats(stats_file: pathlib.Path) -> pd.DataFrame:
    """Extract bin name, contig and scaffold counts, size, GC content, N50 and L50
    from a tab-separated file produced by bbtools's statswrapper.sh
    Arguments:
        stats_file: Path to the file containing assembly statistics
    Returns:
        A selection of the genome statistics.
    """
    stats = pd.read_csv(stats_file, sep="\t", usecols=['n_scaffolds', 'n_contigs', 'scaf_bp', 'scaf_N50', 'scaf_L50',
                                                       'gc_avg', 'filename'])[['filename', 'scaf_bp', 'gc_avg',
                                                                               'n_scaffolds', 'n_contigs', 'scaf_N50',
                                                                               'scaf_L50']]
    # N50 and L50 are switched around in BBTools - switch them back
    stats = stats.rename(columns={"filename": "bin_name", "scaf_N50": "scaffold_L50", "scaf_L50": "scaffold_N50"})
    # extract bin name from file name - TODO: make this more elegant than a double-replace?
    stats['bin_name'] = stats['bin_name'].apply(lambda x: (x.split('/')[-1].replace('.fa', '').replace('.', '_')))
    return stats


def get_checkm_stats(checkm_file: pathlib.Path) -> pd.DataFrame:
    """Extract taxonomy and marker genes from a tab-separated file produced by CheckM.
    Arguments:
        checkm_file: Path to the file containing CheckM results
    Returns:
        Taxonomy and absence/presence data for marker genes for all bins.
    """
    checkm_stats = pd.read_csv(checkm_file, sep="\t")
    # rename amount of marker genes for clarity
    checkm_stats = checkm_stats.rename(columns={'0': '0_markers', '1': '1_marker', '2': '2_markers', '3': '3_markers',
                                        '4': '4_markers', '5': '5_markers', '5+': '5_markers'})
    # unify bin naming for easier comparison
    checkm_stats["Bin Id"] = checkm_stats["Bin Id"].apply(lambda x: x.replace('.', '_'))
    return checkm_stats


def get_drep_stats(mummer_file: pathlib.Path) -> pd.DataFrame:
    """Extract ANI from comma-separated file produced by dRep's Mummer-based ANI step.
    Arguments:
        mummer_file:    Path to Ndb.csv
    Returns:
        Closest bin for each reference with ANI and coverage.
    """
    mummer_anis_original = pd.read_csv(mummer_file, sep=",")
    # fix typos
    mummer_anis_original = mummer_anis_original.rename(lambda x: x.replace('querry', 'query'), axis="columns")
    # unify genome/bin names for easier comparison with other stats
    mummer_anis_original['query'] = mummer_anis_original['query'].apply(lambda x:
                                                                        x.replace('.fa', '').replace('.', '_'))
    mummer_anis_original['reference'] = mummer_anis_original['reference'].apply(lambda x:
                                                                                x.replace('.fa', '').replace('.', '_'))

    # Goal: select one line per primary cluster, where bins are in reference and originals are in query
    # for self-comparisons, if the line doesn't contain "_bin_", delete all but cluster number and 'query'
    # if it does, delete all but cluster number and 'reference'

    # drop all self-comparisons initially
    mummer_anis = mummer_anis_original.drop(mummer_anis_original[mummer_anis_original['query']
                                                                 == mummer_anis_original['reference']].index)
    # still leaves multiple rows per primary cluster - select only these where bins are in reference
    mummer_anis = mummer_anis[mummer_anis['reference'].str.contains('_bin_', regex=False)][~mummer_anis['query'].str.contains('_bin_',
                                                                                                                              regex=False)]
    # re-add all primary clusters only consisting of one member (singletons)
    mummer_anis = pd.concat([mummer_anis, mummer_anis_original.groupby("primary_cluster").filter(lambda x: len(x) == 1)],
                                     ignore_index=True)
    print(mummer_anis)
    mummer_anis = mummer_anis[['query', 'reference', 'ref_coverage', 'query_coverage', 'ani',
                               'primary_cluster']].reset_index(drop=True)
    # self-comparisons are meaningless for our table, replace by NAN
    mummer_anis.loc[(mummer_anis['query'] == mummer_anis['reference'])
                    & (mummer_anis['reference'].str.contains('_bin_', regex=False)),
                    ['query', 'ref_coverage', 'query_coverage', 'ani']] = np.nan

    mummer_anis.loc[(mummer_anis['query'] == mummer_anis['reference'])
                    & (~mummer_anis['reference'].str.contains('_bin_', regex=False)),
                    ['reference', 'ref_coverage', 'query_coverage', 'ani']] = np.nan

    return mummer_anis


def merge_mag_and_ref_stats(mag_stats: pd.DataFrame, ref_stats: pd.DataFrame) -> pd.DataFrame:
    """Concatenate two tables with statistics for the synthetic MAGs and the reference genomes,
    marking which of the two each is.
     Arguments:
         mag_stats: table containing statistics of synthetic MAGs
         ref_stats: table containing statistics of reference genomes
     Returns:
         A merged table with sources of genomes marked.
    """
    mag_stats['genome_type'] = 'synthetic_MAG'
    ref_stats['genome_type'] = 'reference'
    merged_stats = pd.concat([mag_stats, ref_stats], ignore_index=True)
    return merged_stats


if __name__ == "__main__":
    parser = ArgumentParser(description="Extract genome statistics from CheckM, dRep and bbstats output.")
    parser.add_argument("stats", action="store", help="Path to BBstats file giving stats of bins")
    parser.add_argument("genome_stats", action="store", help="Path to BBstats file giving stats of reference genomes")
    parser.add_argument("checkm", action="store", help="Path to CheckM file for bins")
    parser.add_argument("genome_checkm", action="store", help="Path to CheckM file for reference genomes")
    parser.add_argument("drep_mummer", action="store", help="Path to dRep Mummer file (Ndb.csv)")
    parser.add_argument("-o", "--outfile", action="store",
                        help="Name of Excel file to write to (recommended extension: .xlsx) (default: samplestats.xlsx)",
                        default="samplestats.xlsx")
    args = parser.parse_args()

    # get absolute paths
    stats = pathlib.Path(args.stats).resolve()
    original_stats = pathlib.Path(args.genome_stats).resolve()
    checkm = pathlib.Path(args.checkm).resolve()
    original_checkm = pathlib.Path(args.genome_checkm).resolve()
    mummer = pathlib.Path(args.drep_mummer).resolve()
    outfile = pathlib.Path(args.outfile).resolve()

    drep_stats = get_drep_stats(mummer)
    # Extract stats for both MAGs and original genomes
    bb_stats = get_bb_stats(stats)
    reference_stats = get_bb_stats(original_stats)
    # Identify dRep primary clusters for each, selecting each bin only once (given a bin can match multiple genomes)
    bb_stats = pd.merge(bb_stats, drep_stats.drop_duplicates(subset=["reference"]), how="left", left_on="bin_name",
                        right_on="reference")[['bin_name', 'scaf_bp', 'gc_avg', 'n_scaffolds', 'n_contigs',
                                               'scaffold_L50', 'scaffold_N50', 'primary_cluster']]
    reference_stats = pd.merge(reference_stats, drep_stats.drop_duplicates(subset=["query"]), how="left", left_on="bin_name",
                               right_on="query")[['bin_name', 'scaf_bp', 'gc_avg', 'n_scaffolds', 'n_contigs',
                                                  'scaffold_L50', 'scaffold_N50', 'primary_cluster']]

    complete_stats = merge_mag_and_ref_stats(bb_stats, reference_stats)

    # Get CheckM stats for MAGs and original genomes
    checkm_stats = get_checkm_stats(checkm)
    reference_checkm = get_checkm_stats(original_checkm)
    complete_checkm = merge_mag_and_ref_stats(checkm_stats, reference_checkm)

    explanations = pd.DataFrame.from_dict({"genome_type": ["General information: type of genome (synthetic MAG or source genome)"],
                                           "BB_stats: bin_name": ["Name of MetaBAT-generated bin or genome"],
                                           "BB_stats: scaf_bp": ["Basepairs in scaffold(s)"],
                                           "BB_stats: gc_avg": ["Average GC content for bin/genome"],
                                           "BB_stats: n_scaffolds": ["Amount of scaffolds"],
                                           "BB_stats: n_contigs": ["Amount of contigs"],
                                           "BB_stats: scaffold_L50": ["scaffold L50: smallest number of scaffolds that cover half (or more) of the genome together"],
                                           "BB_stats: scaffold_N50": ["scaffold N50: length of the shortest scaffold from the smallest set of scaffolds that covers half (or more) of the genome."],
                                           "BB_stats: primary_cluster": ["dRep cluster based on estimate of average nucleotide identity; included here for ease of comparison."],
                                           "BB_stats: NOTE": ["In the raw BBstats output, N50 and L50 are reversed; here, the commonly used definitions are used instead."],
                                           "CheckM: Bin Id": ["Name of MetaBAT-generated bin or genome"],
                                           "CheckM: Marker lineage": ["Taxon for which specific marker genes could be found"],
                                           "CheckM: # genomes": ["Amount of genomes used to determine marker genes"],
                                           "CheckM: # markers": ["Amount of marker genes (single-copy genes occurring in more than 97 percent of the taxon's genomes) for the given taxon"],
                                           "CheckM: # marker sets": ["Amount of marker gene sets; marker genes are grouped into sets by combining all pairs of collocated marker genes - closer than 5 kb in 95 percent of genomes - which share a gene"],
                                           "CheckM: x_markers": ["Amount of markers occurring x times in the genome"],
                                           "CheckM: completeness": ["Completeness of genome, estimated by number of marker genes that are present"],
                                           "CheckM: contamination": ["Contamination of genome, estimated by number of marker genes occurring more than once"],
                                           "CheckM: Strain heterogeneity": ["Contamination specifically arising from closely related strains, identified via the amount of marker gene duplicates above a threshold of amino acid identity"],
                                           "dRep: query": ["Name of source genome"],
                                           "dRep: reference": ["Name of closest bin"],
                                           "dRep: ref_coverage": ["Percent of bin covered by source genome"],
                                           "dRep: query_coverage": ["Percent of source genome covered by bin"],
                                           "dRep: ani": ["Average nucleotide identity between source genome and bin, calculated by Nucmer alignment"],
                                           "dRep: primary_cluster": ["Cluster containing source genome and bin, determined by Mash-estimated ANI; source and bin need to have ANI of at least 90 percent to be in one cluster"],
                                           "dRep: NOTE": ["If a row contains only a source genome or only a bin, no bin/source genome had an estimated ANI of at least 90 percent."]},
                                          orient="index")

    with pd.ExcelWriter(outfile) as writer:
        complete_stats.to_excel(writer, sheet_name="BB_stats")
        complete_checkm.to_excel(writer, sheet_name="CheckM")
        drep_stats.to_excel(writer, sheet_name="dRep")
        explanations.to_excel(writer, sheet_name="explanations")
