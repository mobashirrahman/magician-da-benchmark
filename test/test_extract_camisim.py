import pathlib
import unittest

import camisim_setup.extract_camisim_data as extract_cami

class TestMetadataFromSample(unittest.TestCase):
    def test_fail_brokenfile(self):
        distribution_file = pathlib.Path(__file__).parent / "data" / "fail_distributions.tsv"
        fail_col = "fail_sample"
        with self.assertRaisesRegex(ValueError, "Incorrect file type, only Genbank and Fasta files can be used"):
            extract_cami.get_camisim_per_sample(distribution_file, fail_col)