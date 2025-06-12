from argparse import ArgumentParser
from pathlib import Path
from textwrap import dedent
from typing import Optional

from Bio import SeqIO


# utility function for getting amount of genomes
def get_file_length(infile: Path) -> int:
    """Count the number of lines in a file.
    Arguments:
        infile: Path to the file to be counted
    Returns:
        The amount of lines in a file
    """
    count = 0
    with open(infile, "r") as file_read:
        for count, line in enumerate(file_read, 1):
            pass
    return count


def get_sample_size(file_record: Path, coverage: Optional[float]=1) -> float:
    """Parse a id_to_genome file, sum up sizes of all fasta files given and,
    if desired, multiply with an average coverage factor.
    Arguments:
        file_record:    Path to id_to_genome file
        coverage:       desired average coverage for the sample
    Returns:
        The sample size for obtaining a given average coverage level, given
        the total size of all genomes examined.
    """
    # sanity check: is the value for coverage above 0?
    if coverage <= 0:
        raise ValueError("Coverage must be above 0")
    # get all paths
    with open(file_record, "r") as idfile:
        fasta_paths = [Path(line.strip().split("\t")[1]) for line in idfile]
    # for each file, get genome size and add up
    total_size = sum([len(record) for fasta_file in fasta_paths
                      for record in SeqIO.parse(fasta_file, "fasta")])
    return round(coverage*(total_size/1000000000), 2)

# TODO: this is becoming a giant almighty function, consider reworking


def generate_config_file(camisim_dir: Path, meta_file: Path, id_file: Path, output_dir: str,
                         readsim: str, readsim_path: Path, samtools_path: Path, sample: str,
                         amount_genomes: int, sample_size: float, error_profiles: Optional[Path] = "",
                         abundance_file: Optional[Path] = "", profile_name: Optional[str] = "mbarc",
                         own_error_basename: Optional[str] = "",
                         own_error_readlength: Optional[int] = "",
                         insert_size: Optional[int] = 270) -> str:
    """Generate a config file for CAMISIM and write it to a specified filename.
    Arguments:
        camisim_dir:            Path to the directory containing CAMISIM
        meta_file:              Path to the metadata file for CAMISIM
        id_file:                Path to the file linking genome IDs to fasta files
        output_dir:             output directory for CAMISIM
        readsim:                read simulator to use
        readsim_path:           Path to the read simulator to use
        samtools_path:          Path to samtools install to use
        sample:                 type of sample to use
        amount_genomes:         amount of genomes in the sample
        sample_size:            sample size in Gbp
        error_profiles:         Path to error profiles, can be blank if using wgsim
        abundance_file:         Path to file listing abundance for genomes, if given
        profile_name:           name of error profile to use; default "mbarc", options "mbarc",
                                "hi", "mi", "hi150", "own" (requires giving own profile & lengths)
        own_error_basename:     name of error profile files, without "[1/2].txt", if using own
        own_error_readlength:   length of reads to simulate with own error profile
        insert_size:            mean insert size (default: 270 bp)
    Returns:
        A CAMISIM config file with the chosen parameters.
    """
    # sanity checks
    # do we get a proper insert size?
    if insert_size <= 0:
        raise ValueError("Mean insert size needs to be above 0.")
    # is the read simulator a valid choice?
    if not readsim in {"art", "wgsim", "nanosim", "pbsim"}:
        raise ValueError("{} is not a valid read simulator. Valid options are art, wgsim, nanosim, pbsim.".format(readsim))
    # is the read simulator not wgsim, but there's no error profile?
    if not error_profiles and not readsim == "wgsim":
        raise ValueError("Error profile can only be omitted with wgsim")
    # is the sample type a valid choice?
    if not sample in {'replicates', 'timeseries_lognormal', 'timeseries_normal', 'differential'}:
        raise ValueError(
            """{} is not a valid sample type. Valid options are 'replicates', 'timeseries_lognormal',\
             'timeseries_normal', 'differential'.""".format(sample))
    # is the profile type a valid choice?
    if not profile_name in {"hi", "mi", "hi150", "mbarc", "own"}:
        raise ValueError("""{} is not a valid type of error profile. \
        Valid options are 'mbarc', 'hi', 'mi', 'hi150', 'own'.""".format(profile_name))
    # is it a custom profile...
    if profile_name == "own":
        # ...and are we using a read simulator that takes these?
        if readsim == "art":
            # ...but no values are given?
            if not own_error_basename:
                raise ValueError("Base profile name must be given when using custom error profile.")
            elif not own_error_readlength:
                raise ValueError("Read length for profile must be given when using custom error profile.")
        else:
            raise ValueError("{simulator} doesn't take custom profiles. Use ART or {simulator}'s builtin profiles.".format(simulator=readsim))
    # conversely, are we getting profile names or read lengths for a default profile?
    if (own_error_readlength or own_error_basename) and not profile_name == "own":
        raise ValueError("Custom error profile files and read lengths are only possible when specifying 'own' profiles.")

    # TODO: introduce more defaults?
    config_string = f'''\
    [Main]
    # maximum number of processes
    max_processors=8
    
    # 0: community design + read simulator,
    # 1: read simulator only
    phase=0
    
    # ouput directory, where the output will be stored (will be overwritten if set in from_profile)
    output_directory={output_dir}
    
    # temporary directory
    temp_directory=/tmp
    
    # gold standard assembly
    gsa=True
    
    # gold standard for all samples combined
    pooled_gsa=True
    
    # anonymize sequences?
    anonymous=False
    
    # compress data (levels 0-9, recommended is 1 the gain of higher levels is not too high)
    compress=1
    
    # id of dataset, used in foldernames and is prefix in anonymous sequences
    dataset_id=RL
    
    # Read Simulation settings, relevant also for from_profile
    [ReadSimulator]
    # which readsimulator to use:
    #           Choice of 'art', 'wgsim', 'nanosim', 'pbsim'
    type={readsim}
    
    # Samtools (http://www.htslib.org/) takes care of sam/bam files. Version 1.0 or higher required!
    # file path to executable
    samtools={samtools_path}
    
    # file path to read simulation executable
    readsim={readsim_path}
    
    #error profiles:
    #for ART:
    #HiSeq 150bp: hi150
    #MBARC-26 150bp: mbarc
    #for wgsim:
    #error rate as <float> (e.g. 0.05 for 5% error rate)
    #blank for nanosim and wgsim
    profile={profile_name}
    
    # Directory containing error profiles (can be blank for wgsim)
    error_profiles={error_profiles}
    
    # Custom error profile filenames if using own error profile
    base_profile_name={own_error_basename}
    
    # Read length for custom error profile if used
    profile_read_length={own_error_readlength}
    
    #paired end read, insert size (not applicable for nanosim)
    fragments_size_mean={insert_size}
    fragment_size_standard_deviation=27
    
    # Only relevant if not from_profile is run:
    [CommunityDesign]
    # optional: give abundance of genomes
    distribution_file_paths={abundance_file}
    # specify the samples size in Giga base pairs
    size={sample_size}
    
    # how many different samples?
    number_of_samples=1
    
    # how many communities
    num_communities=1
    
    # directory containing the taxdump of ncbi, version from 22.02.2017 is shipped
    # "nodes.dmp"
    # "merged.dmp"
    # "names.dmp"
    ncbi_taxdump={camisim_dir}/tools/ncbi-taxonomy_20170222.tar.gz
    
    # the strain simulator for de novo strain creation
    strain_simulation_template={camisim_dir}/scripts/StrainSimulationWrapper/sgEvolver/simulation_dir/
    
    # define communities: [community<integer>]
    [community0]
    # information about all included genomes:
    # can be used for multiple samples
    metadata={meta_file}
    id_to_genome_file={id_file}
    
    # how many genomes do you want to sample over all?
    genomes_total={amount_genomes}
    num_real_genomes={amount_genomes}
    
    # how many genomes per species taxon
    #   (species taxon will be replaced by OTU-cluster later on)
    max_strains_per_otu=1
    ratio=1
    
    # which kind of different samples do you need?
    #   replicates / timeseries_lognormal / timeseries_normal / differential
    mode={sample}
    
    # Part: community design
    # Set parameters of log-normal and normal distribution, number of samples
    # sigma > 0; influences shape (higher sigma -> smaller peak and longer tail),
    log_sigma=2
    
    # mu (real number) is a parameter for the log-scale
    log_mu=1
    
    # do you want to see a distribution before you decide to use it? yes/no
    view=no
    '''
    config_string = dedent(config_string)
    return config_string


if __name__ == "__main__":
    parser = ArgumentParser(description="Generate a config file for CAMISIM from parameters")
    parser.add_argument("camisim_dir", help="Directory containing CAMISIM")
    parser.add_argument("metadata", help="Path to CAMISIM metadata file")
    parser.add_argument("genome_file", help="Path to file containing genome ID to fasta path data")
    parser.add_argument("-f", "--filename", help="Name for the config file (default: camisim_config.ini)",
                        default="camisim_config.ini")
    parser.add_argument('-o', '--out_dir', action="store", help="Output directory for CAMISIM (default: camisim_out)",
                        default="camisim_out")
    parser.add_argument('-a', '--abundance_file', action="store",
                        help="Optional: file giving relative abundance of genomes", default="")
    #parser.add_argument('-c', '--coverage', action="store",
    #                    help="Desired average coverage for the sample (default: 1X)",
    #                    default=1)
    parser.add_argument("--samtools_path", action="store", help="Path to temp file containing samtools path (default: samtools_path.txt)",
                         default="samtools_path.txt")
    parser.add_argument('-s', '--sample_size', action="store",
                        help="Total size of sample in gigabasepairs (default: 1)", default=1)
    parser.add_argument('--insert_size', action="store", help="Mean insert size in bp (default: 270)",
                        default=270, type=int)
    parser.add_argument('--read_sim', action="store", help="Read simulator to use",
                        choices=["art", "wgsim", "nanosim", "pbsim"],
                        default="art")
    parser.add_argument('--read_sim_path', action="store",
                        help="Path to read simulator executable (default: ART shipped with CAMISIM)",
                        default=False)
    parser.add_argument('--sample_type', action="store",
                        help="Type of different samples to be simulated (default: replicates)",
                        choices=["replicates", "timeseries_lognormal", "timeseries_normal", "differential"],
                        default="replicates")
    parser.add_argument("--error_profile", action="store",
                        help="Path to error profiles (overridden by --errorfree). Default: ART error profiles",
                        default=False)
    parser.add_argument("--art_profile_type", action="store", default="mbarc",
                        choices=['mbarc', 'hi', 'mi', 'hi150', 'own'],
                        help="Type of ART error profile: mbarc, hi, mi, hi150, own (default: mbarc)")
    parser.add_argument("--profile_basename", action="store",
                        help="""Base name of custom error profile, if given (name of files without '[1/2].txt');\
                         required with 'own' error profile""")
    parser.add_argument("--profile_readlength", action="store",
                        help="Read length of custom error profile; required with 'own' error profile", type=int)
    parser.add_argument('--errorfree', action="store_true", help="Don't use an error profile (only works with wgsim)")
    args = parser.parse_args()

    # establish location of CAMISIM dir
    camisim_dir = Path(args.camisim_dir).resolve()

    # sanity check insert size
    if args.insert_size <= 0:
        parser.error("Mean insert size needs to be larger than 0 bp.")
    else:
        insert_size = args.insert_size

    # check if no errors is specified without using wgsim - fail early if so
    if args.errorfree:
        if not args.read_sim == "wgsim":
            parser.error("Omitting an error profile is only possible with wgsim.")
        else:
            error_profile = ""
    else:
        # was the path left blank? Then use CAMISIM default
        if not args.error_profile:
            error_profile = camisim_dir / "tools" / "art_illumina-2.3.6" / "profiles"
        else:
            error_profile = Path(args.error_profile).resolve()

    # custom profile arguments default to "" - handling both here for consistency
    profile_basename = ""
    profile_readlength = ""

    # sanity check arguments related to custom profile
    if args.art_profile_type == "own":
        if args.read_sim == "art":
            if not args.profile_basename:
                parser.error("Base name for custom error profile must be given with own error profiles.")
            elif not args.profile_readlength:
                parser.error("Read length for the custom error profile must be given with own error profiles.")
            else:
                profile_basename = args.profile_basename
                profile_readlength = args.profile_readlength
        else:
            parser.error("Supplying custom error profiles is only possible with ART.")

    if (args.profile_basename or args.profile_readlength) and not args.art_profile_type == "own":
        parser.error("To use custom error profiles, specify type of error profile as 'own' and use ART.")

    # establish remaining paths
    metadata = Path(args.metadata).resolve()
    genome_file = Path(args.genome_file).resolve()
    if args.abundance_file:
        abundance_file = Path(args.abundance_file).resolve()
    else:
        abundance_file = args.abundance_file
    # if no read_sim_path, it's the absolute path to CAMISIM's ART
    if not args.read_sim_path:
        read_sim_path = camisim_dir / "tools" / "art_illumina-2.3.6" / "art_illumina"
    else:
        read_sim_path = Path(args.read_sim_path).resolve()

    # get samtools path from file - TODO: move to a function?
    samtools_file = args.samtools_path
    with open(samtools_file, "r") as read_samtools:
        path_to_samtools = Path(read_samtools.readline().strip())
    if not path_to_samtools.exists():
        raise FileNotFoundError("Samtools not found at specified location.")
    # set remaining parameters
    filename = args.filename
    out_dir = args.out_dir
    read_sim = args.read_sim
    sample_type = args.sample_type
    #coverage = float(args.coverage)
    sample_size = float(args.sample_size)
    art_profile_type = args.art_profile_type

    #calculate amount of genomes and sample size
    genomes = get_file_length(genome_file)
    #size_total = get_sample_size(genome_file, coverage)

    config_str = generate_config_file(camisim_dir, metadata, genome_file, out_dir, read_sim, read_sim_path, path_to_samtools, sample_type,
                                      genomes, sample_size, error_profile, abundance_file, art_profile_type,
                                      profile_basename, profile_readlength, insert_size)
    with open(filename, "w") as outfile:
        outfile.write(config_str)
