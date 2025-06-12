import pathlib
import unittest

import numpy as np
import pandas as pd

import generate_summary.make_comparison_table as make_table


class TestMakeTable(unittest.TestCase):
    def test_table_success(self):
        test_file = pathlib.Path("test/data/summary_test_join_ids_new.xlsx")
        true_table = pd.DataFrame({"bin_name": ["test_hiseq_2500_bin_2", "test_hiseq_2500_bin_6",
                                                "test_hiseq_2500_bin_9", "test_hiseq_2500_bin_4"],
                                   "closest_genome": ["Streptomyces_collinus",
                                             "Salinispora_tropica_CNB-440_NC_009380_1",
                                             "Escherichia_coli_str_K-12_substr_MG1655_NC_000913_3",
                                             np.nan],
                                   "ani": [1, 0.9999351, 0.9998399, np.nan],
                                   "bin_coverage": [1, 0.94825554, 0.96552, np.nan],
                                   "source_coverage": [1, 0.98048997, 0.9568264, np.nan],
                                   "scaffold_difference": [3, 4, 3, np.nan],
                                   "contig_difference": [4, 5, 4, np.nan],
                                   "length_difference": [3255698, 4061655, 3255698, np.nan],
                                   "gc_difference": [0.24719, 0.06486, 0.24719, np.nan],
                                   "completeness_difference": [-1.72, 0, -6.39, np.nan]})
        test_table = make_table.create_comparison_table(test_file)
        pd.testing.assert_frame_equal(test_table, true_table)