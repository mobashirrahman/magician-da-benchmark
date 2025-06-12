import logging
import pathlib
import re
import subprocess
import sys

import pandas as pd
import yaml

import pipeline_config

from argparse import ArgumentParser
from typing import List, Optional

# import parameters
default_config_file = pipeline_config.default_config_file
workflow_config = pipeline_config.WORKFLOW_DEFAULT_CONF

logger = logging.getLogger("MAGICIAN")
logger.setLevel(logging.DEBUG)  # in case we need this later
plain_console_log = logging.StreamHandler()
plain_console_log.setLevel(logging.INFO)
plain_format = logging.Formatter("%(name)s: %(message)s")  # very basic format for this one
plain_console_log.setFormatter(plain_format)
logger.addHandler(plain_console_log)

# set defaults
DEFAULT_PROFILE = "mbarc"
DEFAULT_INSERT = 270
DEFAULT_CORES = 6


def make_demo_tempfile(tempfile: pathlib.Path) -> None:
    """Convert the sample distributions file given in test data to use absolute paths for portability
    and save to the current directory.

    Arguments:
        tempfile:   the file to which to save the input file for the demo run to

    """
    original_demo_file = pathlib.Path(__file__).parent / "test" / "data" / "test_genomes" / "sample_distributions.tsv"
    demo_files_to_abs = pd.read_csv(original_demo_file, sep="\t")
    demo_files_to_abs["genomes"] = demo_files_to_abs["genomes"].apply(lambda genome_path:
                                                                      pathlib.Path(genome_path).resolve())
    logger.info("Creating temporary input file %s", tempfile)
    demo_files_to_abs.to_csv(tempfile, sep="\t", index=False)


def get_snake_cmd(input_file, target: str, profile_type: Optional[str] = DEFAULT_PROFILE,
                  profile_base: Optional[str] = "", readlength: Optional[int] = None,
                  insert_size: Optional[int] = DEFAULT_INSERT, cluster_cmd: Optional[str] = "",
                  cores: Optional[int]=DEFAULT_CORES,
                  *snake_params, config_path: pathlib.Path = default_config_file) -> List[str]:
    """Get the Snakemake command with optional configuration parameters.
    Arguments:
        input_file:     File with paths to source genomes, sequence type (plasmid/chromosome) and desired relative
                        copy number in each community to simulate
        target:         rule or output file for Snakemake
        profile_type:   type of ART error profile to use for CAMISIM
        profile_base:   path to own ART error profile for CAMISIM
        readlength:     read length used for own ART error profile
        insert_size:    mean insert size to use for simulation (default: 270)
        cluster_cmd:    command to use for cluster mode
        cores:          the amount of cores Snakemake should use
        snake_params:   parameters to pass to the Snakefile
        config_path:    path to the config file to use with Snakemake

    Returns:
        The command for running Snakemake with the desired parameters.

    Raises:
        ValueError: if arguments contain invalid characters, or if a value that isn't a positive int
                    was given for read length, insert size or amount of cores

    """
    # check all elements of the command
    for input_param in [target, profile_type, profile_base, readlength, insert_size]:
        bad_chars = re.search(r"""[^a-zA-Z0-9"'./_\- ]""", str(input_param))
        if bad_chars:
            raise ValueError("Arguments can only consist of alphanumeric characters, quote marks, ., /, _, - and space.")

    # sanity check parameters
    if insert_size <= 0:
        raise ValueError("Insert size needs to be above 0.")

    core_error = "Amount of cores must be an integer above 0."
    # if we get a bad input for amount of cores we want to complain consistently
    try:
        cores = int(cores)
    except ValueError:
        raise ValueError(core_error)
    if cores <= 0:
        raise ValueError(core_error)

    # we only need to check read length when it's relevant - check explicitly for "not None"
    # so we can complain about read lengths <= 0 specifically
    if profile_type == "own":
        if not (profile_base and (readlength is not None)):
            raise ValueError("Both name of the custom error profile and read length of the error"
                             " profile must be given when using own profiles.")
        if readlength <= 0:
            raise ValueError("Read length needs to be above 0.")

    if (profile_base or readlength) and profile_type != "own":
        raise ValueError("Name of the error profile and read length can only be specified"
                         " when using own profiles.")

    # if custom parameters haven't been given, set to False for later processing in Snakefile
    if not profile_base:
        profile_base = "False"
    if not readlength:
        readlength = "False"

    # get path for snakefile
    snake_path = pathlib.Path(__file__).resolve().parent / "snakefiles" / "Snakefile"

    # load config
    with open(config_path, "r", encoding = "utf-8") as config_file:
        snake_config = yaml.safe_load(config_file)

    # get basic command
    snakemake_cmd = ["snakemake", target, "-s", snake_path]
    # add cluster commands if needed
    if cluster_cmd:
        snakemake_cmd += ["--cluster", cluster_cmd]
    snakemake_cmd += ["--config", 'profile_type="{}"'.format(profile_type),
                     'insert_size={}'.format(insert_size),
                     'samples_file={}'.format(input_file)]
    if profile_type == "own":
        snakemake_cmd += ['profile_name="{}"'.format(profile_base),
                          'readlength={}'.format(readlength)]
    snakemake_cmd += ['--use-conda', '--conda-frontend', snake_config["conda_frontend"],
                     "--configfile", str(config_path),
                     "--cores", str(cores),
                     *snake_params]

    return snakemake_cmd
    

if __name__ == "__main__":
    parser = ArgumentParser(description="Run MAGICIAN to simulate MAGs for a specified community"
                                        " or set of communities.\n"
                                        "Run without arguments for a test run.")
    parser.add_argument("community_file", action="store",
                        help="File with paths to source genomes, sequence type (plasmid/chromosome)"
                             " and their desired relative copy number in each community to simulate "
                             "(one column with organisms' copy numbers per community)")
    parser.add_argument("--target", action="store",
                        help="Desired output file or rule "
                             "(default: MAGs, statistics and summary files for all communities)",
                        default="all_bin_summaries")
    parser.add_argument("--profile_type", action="store",
                        help="Type of ART error profile to use for CAMISIM: mbarc, hi, mi, hi150, own "
                             f"(default: {DEFAULT_PROFILE})",
                        default=DEFAULT_PROFILE, choices=['mbarc', 'hi', 'mi', 'hi150', 'own'])
    parser.add_argument("--profile_name", action="store",
                        help="""Base name of custom error profile, if given (name of files without '[1/2].txt');\
                         required with 'own' error profile""", default="")
    parser.add_argument("--profile_readlength", action="store",
                        help="Read length of custom error profile; "
                             "required with 'own' error profile", default=None)
    parser.add_argument("--insert_size", action="store", type=int, default=DEFAULT_INSERT,
                        help=f"Mean insert size for read simulation (default: {DEFAULT_INSERT})")
    parser.add_argument("--cluster", action="store", default="",
                        help="""For use with snakemake's cluster mode; supply command for submitting jobs as you \
                        would with snakemake.""")
    parser.add_argument("--cores", action="store", default=DEFAULT_CORES,
                        help=f"Amount of cores Snakemake should use (default: {DEFAULT_CORES})")
    parser.add_argument("--config_file", help = "Config file for run")
    parser.add_argument("--snake_flags", nargs='*',
                        help="Flags to be passed to snakemake, enclosed in quotes")
    # if no arguments are given, show usage and offer demo;
    # in this case we're running with the default config
    if len(sys.argv) == 1:
        parser.print_usage()  # TODO: usage or full-on help?
        run_demo = input("Start local example run with sample genomes "
                         f"and output to {pathlib.Path.cwd()}? [y/n] \n"
                         f"Remember to edit {default_config_file}"
                         " to specify your CAMISIM installation.")
        if run_demo == "y":
            # make a copy of the input file with absolute paths so it can run anywhere
            local_tempfile = pathlib.Path.cwd() / "tmp_demo_sample_distributions.tsv"
            make_demo_tempfile(local_tempfile)
            snake_command = get_snake_cmd(local_tempfile, "all_bin_summaries", DEFAULT_PROFILE,
                                          insert_size = DEFAULT_INSERT, cluster_cmd = "",
                                          cores = DEFAULT_CORES)
            snake_run = subprocess.run(snake_command, check=True)
            sys.exit(snake_run.returncode)
        if run_demo == "n":
            sys.exit(0)
        # if we haven't stopped either of these places input was invalid
        logger.error("Invalid option %s", run_demo)
        sys.exit(1)

    # otherwise we have gotten args and need to handle them
    args = parser.parse_args()
    community_file = pathlib.Path(args.community_file)
    target_result = args.target
    profiletype = args.profile_type
    profilename = args.profile_name
    read_length = args.profile_readlength
    insert_size = args.insert_size
    cluster_cmd = args.cluster
    snake_cores = args.cores
    snake_flags = []
    if args.snake_flags:
        snake_flags = args.snake_flags[0].split()
    if args.config_file:
        default_config_file = pathlib.Path(args.config_file).resolve()

    snake_command = get_snake_cmd(community_file, target_result, profiletype, profilename,
                                  read_length, insert_size,
                                  cluster_cmd, snake_cores, *snake_flags)
    subprocess.run(snake_command, check=True)

