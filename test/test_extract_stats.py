import unittest

import numpy as np
import pandas as pd

from pathlib import Path

import generate_summary.extract_stats as summary_stats


class TestBBStats(unittest.TestCase):
    def test_get_stats(self):
        true_stats = pd.DataFrame({"bin_name": ["sample14_bin_1", "sample14_bin_2"], "scaf_bp": [8667507, 8277261],
                                   "gc_avg": [0.6, 0.5],
                                    "n_scaffolds": [4, 5], "n_contigs": [5, 6], "scaffold_L50": [2, 3],
                                    "scaffold_N50": [1667507, 1277261]})
        stats_file = Path('test/data/fake_bb.csv')
        test_stats = summary_stats.get_bb_stats(stats_file)
        pd.testing.assert_frame_equal(true_stats, test_stats, check_exact=False)

    # Test parsing of reference genome names
    def test_get_stats_from_reference(self):
        true_reference = pd.DataFrame({"bin_name": ["Bacillus_cereus_ATCC_14579_NC_004722_1",
                                                    "Bacillus_subtilis_subsp_subtilis_str_168_NC_000964_3"],
                                       "scaf_bp": [5411809, 4215606],
                                       "gc_avg": [0.35281, 0.43514],
                                       "n_scaffolds": [1, 1],
                                       "n_contigs": [1, 1],
                                       "scaffold_L50": [1, 1],
                                       "scaffold_N50": [5411809, 4215606]})
        reference_file = Path('test/data/fake_refgenomes.csv')
        test_reference = summary_stats.get_bb_stats(reference_file)
        pd.testing.assert_frame_equal(true_reference, test_reference, check_exact=False)


class TestCheckMStats(unittest.TestCase):
    def test_get_checkm(self):
        true_checkm = pd.DataFrame({"Bin Id": ["sample13_bin_1", "sample13_bin_2"],
                                "Marker lineage": ["f__Enterobacteriaceae (UID5124)", "f__Streptomycetaceae (UID2048)"],
                                "# genomes": [134, 60], "# markers": [1173, 460], "# marker sets": [336, 233],
                                "0_markers": [79, 24], "1_marker": [1068, 429], "2_markers": [26, 7],
                                "3_markers": [0, 0], "4_markers": [0, 0], "5_markers": [0,0],
                                "Completeness": [93.13, 95.38], "Contamination": [3.11, 1.37],
                                "Strain heterogeneity": [38.46, 42.86]})
        checkm_file = Path('test/data/fake_checkm.txt')
        test_checkm = summary_stats.get_checkm_stats(checkm_file)
        pd.testing.assert_frame_equal(true_checkm, test_checkm, check_exact=False)


class TestDRepStats(unittest.TestCase):
    def test_get_drep(self):
        true_drep = pd.DataFrame({"query": ["Streptomyces_coelicolor_A32_NC_003888_3",
                                             "Salinispora_tropica_CNB-440_NC_009380_1",
                                             "Escherichia_coli_str_K-12_substr_MG1655_NC_000913_3",
                                             "Klebsiella_pneumoniae_subsp_pneumoniae_CP009208_1",
                                             np.nan],
                                  "reference": ["test_hiseq_2500_bin_2",
                                                "test_hiseq_2500_bin_6",
                                                "test_hiseq_2500_bin_9",
                                                np.nan,
                                                "test_hiseq_2500_bin_4"],
                                  "ref_coverage": [0.9987903, 0.94825554, 0.96552, np.nan, np.nan],
                                  "query_coverage": [0.95382077, 0.98048997, 0.9568264, np.nan, np.nan],
                                  "ani": [0.9999186, 0.9999351, 0.9998399, np.nan, np.nan],
                                  "primary_cluster": [1, 2, 3, 4, 5]})
        drep_mummer = Path('test/data/Ndb.csv')
        test_drep = summary_stats.get_drep_stats(drep_mummer)
        pd.testing.assert_frame_equal(test_drep, true_drep, check_exact=False)


class TestMergeStats(unittest.TestCase):
    def test_merge(self):
        stats = pd.DataFrame({"bin_name": ["sample14_bin_1", "sample14_bin_2"], "scaf_bp": [8667507, 8277261],
                              "gc_avg": [0.6, 0.5],
                              "n_scaffolds": [4, 5], "n_contigs": [5, 6], "scaffold_L50": [2, 3],
                              "scaffold_N50": [1667507, 1277261]})
        reference = pd.DataFrame({"bin_name": ["Bacillus_cereus_ATCC_14579_NC_004722_1",
                                               "Bacillus_subtilis_subsp_subtilis_str_168_NC_000964_3"],
                                  "scaf_bp": [5411809, 4215606],
                                  "gc_avg": [0.35281, 0.43514],
                                  "n_scaffolds": [1, 1],
                                  "n_contigs": [1, 1],
                                  "scaffold_L50": [1, 1],
                                  "scaffold_N50": [5411809, 4215606]})
        true_merged = pd.DataFrame({"bin_name": ["sample14_bin_1", "sample14_bin_2",
                                                 "Bacillus_cereus_ATCC_14579_NC_004722_1",
                                                "Bacillus_subtilis_subsp_subtilis_str_168_NC_000964_3"],
                                    "scaf_bp": [8667507, 8277261, 5411809, 4215606],
                                    "gc_avg": [0.6, 0.5, 0.35281, 0.43514],
                                    "n_scaffolds": [4, 5, 1, 1], "n_contigs": [5, 6, 1, 1],
                                    "scaffold_L50": [2, 3, 1, 1],
                                    "scaffold_N50": [1667507, 1277261, 5411809, 4215606],
                                    "genome_type": ["synthetic_MAG", "synthetic_MAG", "reference", "reference"]})
        test_merged = summary_stats.merge_mag_and_ref_stats(stats, reference)
        pd.testing.assert_frame_equal(test_merged, true_merged, check_exact=False)
