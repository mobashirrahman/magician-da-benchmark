import re
import shutil

from argparse import ArgumentParser
from pathlib import Path

import pandas as pd

from Bio import SeqIO


def get_camisim_per_sample(samples_file: Path, sample_col: str):
    """From a tab-separated table giving genbank files and their abundance in a given sample,
    create CAMISIM metadata, genome and abundance files.
    Arguments:
        samples_file:   Path to .tsv file
        sample_col:     column name for sample
    Returns:
        A tab-separated metadata file containing genome ID, OTU, NCBI taxid and novelty category,
        a tab-separated file listing genome ID and abundance,
        a directory with a fasta file with the sequence of each gbk in the table with nonzero abundance, and
        a tab-separated file listing the genome ID and path of each fasta file.
    """
    # TODO: make testable! Split out record parsing?
    samples_table = pd.read_csv(samples_file, sep="\t", index_col=False)
    # check for any genomes with an abundance of 0 in the current sample, remove these
    samples_table = samples_table.loc[samples_table[sample_col] != 0]
    # create metadata and fasta files
    records = [Path(record) for record in samples_table['genomes']]
    # set up metadata file with "genome_ID\tOTU\tNCBI_ID\tnovelty_category"
    metadata = ["genome_ID\tOTU\tNCBI_ID\tnovelty_category"]
    # initialize OTU count to 1 as we want all to be included
    otu_count = 1
    # initialize id_to_genome_file to blank
    id_to_genome = []
    # save record IDs separately as well
    record_ids = []
    # check if fasta dir exists
    fasta_dir = "camisim_fasta_{}".format(sample_col)
    if not Path(fasta_dir).is_dir():
        Path(fasta_dir).mkdir()
    # for each input file:
    for gene_file in records:
        taxon_id = ""
        strain_name = ""
        file_type = None
        # identify whether it's a fasta or genbank file
        with open(gene_file, "r") as infile:
            line = infile.readline()
            if line.startswith("LOCUS"):
                file_type = "genbank"
                # for genbanks, extract taxonomic ID and strain name here
                while line and not (taxon_id and strain_name):
                    if "taxon:" in line:
                        taxon_parts = line.strip().split(":")
                        taxon_id = taxon_parts[1].replace('"', "")
                    if "/strain" in line or "/isolate" in line:
                        strain_parts = line.strip().split("=")
                        strain_name = strain_parts[1].replace('"', "")
                    line = infile.readline()
            elif line.startswith(">"):
                file_type = "fasta"
            else:
                raise ValueError("Incorrect file type, only Genbank and Fasta files can be used")
        # for Genbank files, extract files
        if file_type == "genbank":
            record = SeqIO.read(gene_file, "genbank")  # TODO: this only works with no contigs - adapt for contigs later!
            #     extract ID: source and identifier
            # add strain identifier to source if it's not already there
            if strain_name in record.annotations["source"]:
                strain_source = record.annotations["source"]
            else:
                strain_source = record.annotations["source"] + " " + strain_name
            # clean up source: remove all characters that aren't alphanumeric, - or _
            sanitized_source = re.sub(r"[^A-Za-z0-9_\-]", "", strain_source.replace(" ", "_"))
            record_id = "{}_{}".format(sanitized_source, record.id.replace(".", "_"))
            # create and save FASTA file
            fasta_path = Path(fasta_dir,
                              '{}.fa'.format(record_id)).resolve()  # resolves as far as possible, appends the rest
            SeqIO.convert(gene_file, "genbank", fasta_path, "fasta")
        else:
            # for Fasta file, record ID is cleaned file name
            record_id = re.sub(r"[^A-Za-z0-9_\-]", "", gene_file.stem.replace(" ", "_"))
            # we cannot extract a taxon, so specify a placeholder - bacteria
            taxon_id = "2"
            # write a copy of the file under the cleaned name, ending in .fa
            fasta_path = Path(fasta_dir,
                              '{}.fa'.format(record_id)).resolve()  # resolves as far as possible, appends the rest
            shutil.copyfile(gene_file, fasta_path)
        record_ids.append(record_id)
        # create and append metadata line
        metadata.append("{}\t{}\t{}\tknown_strain".format(record_id, otu_count, taxon_id))
        otu_count += 1

        # id to file line
        id_to_genome.append("{}\t{}".format(record_id, fasta_path))

    # write metadata/id to fasta files
    # check if dir exists
    if not Path('camisim_configfiles').is_dir():
        Path('camisim_configfiles').mkdir()
    # write metadata file
    with open(Path("camisim_configfiles", "metadata_{}".format(sample_col)), "w") as meta_file:
        meta_file.write("\n".join(metadata))
    # write id_to_genome_file
    with open(Path("camisim_configfiles", "id_to_genome_file_{}".format(sample_col)), "w") as id_file:
        id_file.write("\n".join(id_to_genome))
    # extract genome IDs
    samples_table['genome_id'] = record_ids
    with open(Path("camisim_configfiles", "id_to_distributions_{}".format(sample_col)), "w") as abundance_file:
        abundance_file.write(samples_table.to_csv(sep="\t", header=False, index=False, columns=['genome_id', sample_col]))


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Prepare CAMISIM metadata and id to file mapping files and create FASTA files from gbks")
    parser.add_argument("sample_file", help="Tab-separated file with abundances for samples")
    parser.add_argument("sample_column", help="Column name of sample to extract")
    args = parser.parse_args()
    sample_file = args.sample_file
    sample_column = args.sample_column
    get_camisim_per_sample(sample_file, sample_column)
