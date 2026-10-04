"""Regression tests for scientific correctness, rather than tool command spelling."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest
from magician.config import load, combinations
from magician.design import allocate, prepare_references, make_design
from magician.catalog import match_paf, union_length
from magician.matrices import build, normalized_metrics
from magician.evaluation import average_precision, confusion, evaluate_case, aggregate
from magician.io import save_table, table, write_json
from magician.reads import count_sam
from magician.report import preflight


@pytest.fixture
def cfg(tmp_path):
    return load(overrides={"output_dir": str(tmp_path / "out"), "cache_dir": str(tmp_path / "cache"),
                           "genomes": {"synthetic_length": 12000}})


def test_independent_reproducible_design_and_null(cfg, tmp_path):
    refs = tmp_path / "refs"
    prepare_references(cfg, refs)
    for case in ["null_s42", "spiked_s42"]:
        for repeat in [1, 2]:
            make_design(cfg, case, refs, tmp_path / f"{case}_{repeat}")
        a, b = tmp_path / f"{case}_1", tmp_path / f"{case}_2"
        assert (a / "allocations.tsv").read_bytes() == (b / "allocations.tsv").read_bytes()
        alloc = table(a / "allocations.tsv")
        meta = table(a / "samples.tsv")
        assert alloc.groupby("sample_id").read_pairs.sum().to_dict() == meta.set_index("sample_id").read_pairs.to_dict()
        proportions = alloc.pivot(index="genome_id", columns="sample_id", values="genomic_proportion")
        assert not np.allclose(proportions.C01, proportions.C02)
    null = table(tmp_path / "null_s42_1/truth.tsv")
    assert not null.is_da.any()
    truth = table(tmp_path / "spiked_s42_1/truth.tsv")
    unchanged = ~truth.is_da
    assert unchanged.any()
    assert np.allclose(truth.loc[unchanged, "baseline_proportion"], truth.loc[unchanged, "treatment_proportion"])
    assert truth.true_log2fc.min() < 0 < truth.true_log2fc.max()
    assert truth.baseline_proportion.sum() == pytest.approx(truth.treatment_proportion.sum())


def test_integer_pair_allocation():
    assert allocate(7, [0.5, 0.3, 0.2]).tolist() == [4, 2, 1]


def test_count_fragments_excludes_secondary_second_mate_and_low_mapq():
    def sam(flag, mapq=60):
        return f"read\t{flag}\tcontig\t1\t{mapq}\t150M\t=\t200\t350\tA\tI\n"
    values = count_sam([sam(99), sam(147), sam(355), sam(99, 0), sam(2147)], {"contig": "genome"}, 20)
    assert values == {"genome": 1}


def test_tied_ranks_do_not_create_artificial_ap():
    assert average_precision([True, False], [0, 0]) == pytest.approx(0.5)
    assert average_precision([False, True], [0, 0]) == pytest.approx(0.5)
    assert average_precision([True, False], [1, 0]) == 1
    assert confusion([True, False], [True, True])["fdr"] == 0.5


def test_matching_union_ambiguity_and_multiple_bins(cfg, tmp_path):
    assert union_length([(0, 60), (40, 100)]) == 100
    save_table(pd.DataFrame({"feature_id": ["MAG_1", "MAG_2", "MAG_3"], "length_bp": [100, 100, 100]}), tmp_path / "mags.tsv")
    save_table(pd.DataFrame({"genome_id": ["A", "B"], "length_bp": [200, 200]}), tmp_path / "refs.tsv")
    def paf(mag, source):
        return f"{mag}__c1\t100\t0\t100\t+\t{source}__c1\t200\t0\t100\t100\t100\t60\n"
    (tmp_path / "hits.paf").write_text(paf("MAG_1", "A") + paf("MAG_2", "A") + paf("MAG_3", "A") + paf("MAG_3", "B"))
    match_paf(cfg, tmp_path / "hits.paf", tmp_path / "mags.tsv", tmp_path / "refs.tsv", tmp_path / "matching.tsv")
    result = table(tmp_path / "matching.tsv").set_index("feature_id")
    assert result.loc["MAG_3", "match_status"] == "ambiguous"
    assert result.loc["MAG_3", "source_id"] == ""
    assert result.loc["MAG_1", "multiple_bins_per_source"]
    assert result.loc["MAG_1", "source_completeness"] == 0.5


def test_length_corrected_metrics_and_zero_library():
    counts = pd.DataFrame({"s1": [10, 20], "s2": [0, 0]}, index=["A", "B"])
    metrics = normalized_metrics(counts, pd.Series([100, 200], index=counts.index), pd.Series([30, 1], index=counts.columns))
    assert metrics["relative"].s1.tolist() == [0.5, 0.5]
    assert metrics["relative"].s2.sum() == 0
    assert metrics["tpm"].s1.sum() == pytest.approx(1e6)


def test_zero_mag_matrix_has_no_fabricated_features(cfg, tmp_path):
    samples = ["C01", "T01"]
    save_table(pd.DataFrame({"sample_id": samples, "group": ["Control", "Treatment"], "read_pairs": [10, 10]}), tmp_path / "samples.tsv")
    save_table(pd.DataFrame(columns=["feature_id", "length_bp"]), tmp_path / "features.tsv")
    save_table(pd.DataFrame(columns=["contig_id", "feature_id", "length_bp"]), tmp_path / "assign.tsv")
    paths = []
    for sample in samples:
        path = tmp_path / f"{sample}.tsv"
        save_table(pd.DataFrame({"feature_id": ["contig"], "count": [5], "length_bp": [100]}), path)
        paths.append(path)
    build(cfg, "mags", paths, tmp_path / "samples.tsv", tmp_path / "features.tsv", tmp_path / "assign.tsv", tmp_path / "matrices")
    matrix = pd.read_csv(tmp_path / "matrices/counts.tsv", sep="\t", index_col=0)
    assert matrix.shape == (0, 2)
    assert json.loads((tmp_path / "matrices/matrix.json").read_text())["status"] == "no_mags"


def test_count_methods_never_get_normalized_metrics(cfg):
    pairs = combinations(cfg)
    assert ("counts", "deseq2") in pairs
    assert ("tpm", "deseq2") not in pairs
    assert ("tpm", "wilcoxon_clr") in pairs


def test_preflight_rejects_over_budget(cfg, tmp_path):
    refs = tmp_path / "references.tsv"
    save_table(pd.DataFrame({"genome_id": ["A"], "length_bp": [1000000]}), refs)
    design = tmp_path / "design.json"
    write_json(design, {"actual_read_pairs": 10**10})
    with pytest.raises(ValueError, match="storage exceeds"):
        preflight(cfg, refs, [design], tmp_path / "preflight.json")


def test_unknown_configuration_and_duplicate_seeds_fail():
    with pytest.raises(ValueError, match="unique"):
        load(overrides={"experiment": {"seeds": [42, 42]}})
    with pytest.raises(ValueError, match="Unknown DA"):
        load(overrides={"analysis": {"methods": ["invented"]}})
