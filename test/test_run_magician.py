import pathlib
import unittest

import pandas as pd

import run_magician


class TestMakeDemoFile(unittest.TestCase):
    demo_tempfile = pathlib.Path(__file__).parent / "data" / "test_tmp_demo_distributions.tsv"
    @classmethod
    def setUpClass(cls):
        cls.demo_tempfile.unlink(missing_ok=True)

    @classmethod
    def tearDownClass(cls):
        cls.demo_tempfile.unlink(missing_ok=True)

    def test_write_demo_file(self):
        """Successfully convert relative paths in the example file to absolute paths."""
        expected_results = pd.DataFrame(data={"genomes": [str(pathlib.Path(__file__).resolve().parent
                                                          / "data"/ "test_genomes"
                                                          / "Bifidobacterium_choerinum_FMB-1_CP018044.1.gb"),
                                                          str(pathlib.Path(__file__).resolve().parent
                                                          / "data" / "test_genomes"
                                                          / "Enterococcus_faecium_Ef_aus00233_LT598663.1.gb"),
                                                          str(pathlib.Path(__file__).resolve().parent
                                                          / "data" / "test_genomes"
                                                          / "Mycoplasma_pneumoniae_C267_NZ_CP014267.gb")
                                                          ],
                                              "seq_type": ["chromosome", "chromosome", "chromosome"],
                                              "test_magician": [1, 2, 3]})
        log_msg = f'INFO:MAGICIAN:Creating temporary input file {self.demo_tempfile}'
        with self.assertLogs("MAGICIAN", level="INFO") as logged:
            run_magician.make_demo_tempfile(self.demo_tempfile)
            assert log_msg in logged.output
        test_result = pd.read_csv(self.demo_tempfile, sep="\t")
        pd.testing.assert_frame_equal(test_result, expected_results)
        

class TestRunMagician(unittest.TestCase):
    snake_path = pathlib.Path(__file__).resolve().parent.parent / "snakefiles" / "Snakefile"
    profile_type = "mbarc"
    profile_base = ""
    readlength = None
    insert_size = 270
    cluster_cmd = ""
    cores = 6
    distributions_file = pathlib.Path(__file__).parent / "data" / "test_distribution_file.tsv"

    def test_default_command_success(self):
        """Run Snakemake with default settings."""
        expected_command = ["snakemake", "all_bin_summaries", "-s", self.snake_path,
                            "--config", 'profile_type="mbarc"',
                            'insert_size=270', f"samples_file={self.distributions_file}",
                            "--use-conda",
                            "--conda-frontend", "conda",
                            "--configfile", str(run_magician.default_config_file),
                            "--cores", "6", "-n"]
        snake_flags = ["-n"]
        test_command = run_magician.get_snake_cmd(self.distributions_file, "all_bin_summaries",
                                                  self.profile_type,
                                                  self.profile_base, self.readlength, self.insert_size,
                                                  self.cluster_cmd, self.cores, *snake_flags)
        assert test_command == expected_command

    def test_use_config_success(self):
        """Run Snakemake and pass a config file."""
        config_file = pathlib.Path(__file__).parent / "data" / "test_wrapper" / "test_config.yml"

        expected_command = ["snakemake", "all_bin_summaries", "-s", self.snake_path,
                            "--config", 'profile_type="mbarc"',
                            'insert_size=270', f"samples_file={self.distributions_file}",
                            "--use-conda",
                            "--conda-frontend", "mamba",
                            "--configfile", str(config_file),
                            "--cores", "6",
                            "-n"]
        snake_flags = ["-n"]
        test_command = run_magician.get_snake_cmd(self.distributions_file, "all_bin_summaries",
                                                  self.profile_type, self.profile_base,
                                                  self.readlength, self.insert_size,
                                                  self.cluster_cmd, self.cores,
                                                  *snake_flags,
                                                  config_path = config_file)
        assert test_command == expected_command

    def test_run_cluster_command(self):
        """Run Snakemake in cluster mode"""
        expected_command = ["snakemake", "all_bin_summaries", "-s", self.snake_path,
                            "--cluster", "qsub -pe threaded {threads}",
                            "--config", 'profile_type="mbarc"',
                            'insert_size=270', f"samples_file={self.distributions_file}",
                            "--use-conda",
                            "--conda-frontend", "conda",
                            "--configfile", str(run_magician.default_config_file), "--cores", "6",
                            "-n"]
        snake_flags = ["-n"]
        cluster_cmd = "qsub -pe threaded {threads}"
        test_command = run_magician.get_snake_cmd(self.distributions_file, "all_bin_summaries",
                                                  self.profile_type, self.profile_base,
                                                  self.readlength, self.insert_size, cluster_cmd,
                                                  self.cores,
                                                  *snake_flags)
        assert test_command == expected_command

    def test_change_settings(self):
        """Use non-default settings for various parameters."""
        expected_command = ["snakemake", "all_bin_summaries", "-s", self.snake_path,
                            "--config", 'profile_type="mbarc"',
                            'insert_size=500',  f"samples_file={self.distributions_file}",
                            "--use-conda",
                            "--conda-frontend", "conda",
                            "--configfile", str(run_magician.default_config_file),
                            "--cores", "6", "-n"]
        snake_flags = ["-n"]
        insert_size = 500
        test_command = run_magician.get_snake_cmd(self.distributions_file, "all_bin_summaries",
                                                  self.profile_type, self.profile_base,
                                                  self.readlength, insert_size, self.cluster_cmd,
                                                  self.cores,
                                                  *snake_flags)
        assert test_command == expected_command

    def test_set_own_profile(self):
        """Successfully set a custom error profile."""
        expected_command = ["snakemake", "all_bin_summaries", "-s", self.snake_path,
                            "--config", 'profile_type="own"',
                            'insert_size=270', f"samples_file={self.distributions_file}",
                            'profile_name="test_R"', "readlength=150",
                            "--use-conda",
                            "--conda-frontend", "conda",
                            "--configfile", str(run_magician.default_config_file),
                            "--cores", "6", "-n"]
        snake_flags = ["-n"]
        profile_name = "own"
        profile_base = "test_R"
        readlength = 150
        test_command = run_magician.get_snake_cmd(self.distributions_file, "all_bin_summaries",
                                                  profile_name,
                                                  profile_base, readlength, self.insert_size,
                                                  self.cluster_cmd, self.cores, *snake_flags)
        assert test_command == expected_command

    def test_bad_readlength(self):
        """Catch bad read length."""
        snake_flags = ["-n"]
        profile_name = "own"
        profile_base = "test_R"
        readlength = 0
        error_msg = r"Read length needs to be above 0."
        with self.assertRaisesRegex(ValueError, error_msg):
            run_magician.get_snake_cmd(self.distributions_file, "all_bin_summaries",
                                       profile_name,
                                       profile_base, readlength, self.insert_size,
                                       self.cluster_cmd, self.cores, *snake_flags)


    def test_bad_insertsize(self):
        """Catch bad insert size."""
        snake_flags = ["-n"]
        insert_size = 0
        with self.assertRaisesRegex(ValueError, r"Insert size needs to be above 0."):
            run_magician.get_snake_cmd(self.distributions_file, "all_bin_summaries",
                                       self.profile_type,
                                       self.profile_base, self.readlength, insert_size, self.cluster_cmd, *snake_flags)

    def test_bad_file_name(self):
        bad_result = "test_ä.txt"
        with self.assertRaisesRegex(ValueError,
                                    r"""Arguments can only consist of alphanumeric characters, quote marks, \., \/, \_, \- and space."""):
            run_magician.get_snake_cmd(self.distributions_file, bad_result)

    def test_bad_cores(self):
        """Handle wrong dtype for cores."""
        error_msg = r"Amount of cores must be an integer above 0\."
        bad_cores = "six"
        with self.assertRaisesRegex(ValueError, error_msg):
            run_magician.get_snake_cmd(self.distributions_file, "all_bin_summaries", self.profile_type,
                                       self.profile_base, self.readlength, self.insert_size, self.cluster_cmd,
                                       bad_cores)
    def test_invalid_cores(self):
        """Handle invalid value for cores."""
        error_msg = r"Amount of cores must be an integer above 0\."
        bad_cores = 0
        with self.assertRaisesRegex(ValueError, error_msg):
            run_magician.get_snake_cmd(self.distributions_file, "all_bin_summaries", self.profile_type,
                                       self.profile_base, self.readlength, self.insert_size, self.cluster_cmd,
                                       bad_cores)
    def test_missing_profile_info(self):
        result = "test.txt"
        profile = "own"
        with self.assertRaisesRegex(ValueError,
                                    "Both name of the custom error profile and read length of the error profile must be given when using own profiles."):
            run_magician.get_snake_cmd(self.distributions_file, result, profile)

    def test_too_much_info(self):
        result = "test.txt"
        profile = "mbarc"
        extra_name = "TestR"
        with self.assertRaisesRegex(ValueError,
                                    "Name of the error profile and read length can only be specified when using own profiles."):
            run_magician.get_snake_cmd(self.distributions_file, result, profile, extra_name)