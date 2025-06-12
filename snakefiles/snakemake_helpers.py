import pathlib

import pandas as pd

# hacky helper function for identifying whether plasmids are present
def check_plasmids(samples_file: pathlib.Path, sample: str) -> bool:
    """Identify whether a set of input sequences includes plasmids
    from information in sample file.
    Arguments:
        samples_file: path to file listing sample compositions
        sample:       name of the sample to examine
    Returns:
        True if sample contains plasmids, otherwise False
    """
    samples_table = pd.read_csv(samples_file,sep="\t",index_col=False)
    # check for any genomes with an abundance of 0 in the current sample, remove these
    samples_table = samples_table.loc[samples_table[sample] != 0]
    # does the sequence type column contain plasmids?
    plasmid_check = samples_table["seq_type"].str.contains("plasmid")
    return plasmid_check.any()