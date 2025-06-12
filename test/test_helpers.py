import pathlib
import unittest

import snakefiles.snakemake_helpers as snakehelper

class TestPlasmidCheck(unittest.TestCase):
    def test_plasmid_present(self):
        distribution_file = pathlib.Path(__file__).parent / "data" / "plasmid_distributions.tsv"
        no_plasmid_sample = "plasmidfree"
        yes_plasmid_sample = "plasmids"
        assert snakehelper.check_plasmids(distribution_file, yes_plasmid_sample)
        assert not snakehelper.check_plasmids(distribution_file, no_plasmid_sample)