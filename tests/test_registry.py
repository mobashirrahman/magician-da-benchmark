"""Registry contracts, declared zero policies, truth lanes and matrix-only entry."""
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from magician.config import load, zero_policies, registry
from magician.registry import validate as validate_registry
from magician.truth import apply_zero_policy, centre_log_ratio, LANES
from magician.matrices import derive, OTHER
from magician import matrix
from magician.io import read_json, table


def test_registry_declares_a_compatible_job_grid():
    cfg = load()
    reg = registry(cfg)
    assert ("counts", "deseq2") in reg.select(cfg["analysis"]["methods"], cfg["analysis"]["metrics"])
    assert ("tpm", "deseq2") not in reg.select(cfg["analysis"]["methods"], cfg["analysis"]["metrics"])
    assert ("relative", "maaslin3_abundance") in reg.select(cfg["analysis"]["methods"] + ["maaslin3_abundance"], cfg["analysis"]["metrics"])


def test_default_registry_entry_points_keep_the_previous_grid_size():
    cfg = load()
    pairs = registry(cfg).select(cfg["analysis"]["methods"], cfg["analysis"]["metrics"])
    assert len(pairs) == 13


def test_gamma_variants_share_a_family_and_a_scale_model_does_not():
    reg = registry(load())
    assert reg.entry("aldex2")["method_family"] == reg.entry("aldex2_g0")["method_family"]
    assert len({reg.entry(name)["variant"] for name in ("aldex2", "aldex2_g0")}) == 2
    # ALDEx3 adds a scale model, so it is a different family and never a variant of it.
    assert reg.entry("aldex3")["method_family"] != reg.entry("aldex2")["method_family"]


def test_unknown_method_and_disabled_status_are_rejected():
    with pytest.raises(ValueError, match="Unknown DA method"):
        load(overrides={"analysis": {"methods": ["invented"]}})
    with pytest.raises(ValueError, match="does not allow scheduling"):
        load(overrides={"analysis": {"methods": ["metageenome_gee"]}})


def test_registry_endpoint_must_have_a_truth_lane(tmp_path):
    cfg = load()
    reg = registry(cfg)
    reg.endpoints["made_up"] = {"lane": "not_a_lane", "description": "x"}
    with pytest.raises(ValueError, match="unknown truth lane"):
        validate_registry(reg, cfg)


def test_registry_entry_needs_its_environment_and_adapter(tmp_path):
    cfg = load()
    reg = registry(cfg)
    reg.env_dir = tmp_path
    with pytest.raises(ValueError, match="missing environment"):
        validate_registry(reg, cfg)


def test_count_scale_pseudocount_is_unit_aware():
    values = pd.DataFrame([[10.0, 20.0], [30.0, 40.0]])
    library = np.array([100.0, 200.0])
    measurement = apply_zero_policy(values, library,
                                    {"zero_handling": "pseudocount", "pseudocount": 0.5,
                                     "pseudocount_scale": "measurement"})
    count_scale = apply_zero_policy(values, library,
                                    {"zero_handling": "pseudocount", "pseudocount": 0.5,
                                     "pseudocount_scale": "count"})
    # The same constant means different things in proportions and in expected counts.
    assert measurement.iloc[0, 0] == 10.5
    assert count_scale.iloc[0, 0] == 10.0 + 0.5 / 100
    assert count_scale.iloc[0, 1] == 20.0 + 0.5 / 200


def test_dirichlet_policy_replaces_only_zeros():
    # Samples by features, the orientation every Python truth helper uses.
    values = pd.DataFrame([[0.0, 4.0], [6.0, 0.0]])
    library = np.array([10.0, 10.0])
    replaced = apply_zero_policy(values, library, {"zero_handling": "dirichlet", "dirichlet_alpha": 1.0})
    # A zero becomes the posterior mean of a Dirichlet fitted to that sample's counts;
    # observed values are left alone, so only the imputed mass is added.
    assert replaced.iloc[0, 0] == pytest.approx(1 / 6)
    assert replaced.iloc[0, 1] == 4.0
    assert replaced.iloc[1, 0] == 6.0
    assert replaced.iloc[1, 1] == pytest.approx(1 / 8)


def test_centre_log_ratio_is_the_clr_of_the_contract():
    values = pd.DataFrame([[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]])
    clr = centre_log_ratio(values)
    # Each sample is centred, not each feature.
    assert np.allclose(clr.mean(axis=1).to_numpy(), 0.0, atol=1e-12)
    expected = np.log2(values.to_numpy())
    assert np.allclose(clr.to_numpy(), expected - expected.mean(axis=1, keepdims=True), atol=1e-12)


def test_declared_policies_all_validate():
    names = {policy["name"] for policy in zero_policies(load())}
    assert names == {"pseudocount_measurement", "pseudocount_count", "dirichlet"}


def test_truth_lanes_cover_every_declared_geometry(mini, cfg):
    lanes = mini["spiked_s42"]["lanes"]
    assert set(lanes.lane) == set(LANES)
    # The CLR lane exists once per input metric, zero policy and reserved-category variant.
    clr = lanes.loc[lanes.lane.eq("clr_log_ratio")]
    assert clr.input_metric.nunique() == len(cfg["analysis"]["metrics"])
    assert set(clr.zero_policy) >= {"none", "pseudocount_measurement", "pseudocount_count", "dirichlet"}
    assert clr.reserve_other.astype(str).str.lower().eq("false").any()


def test_absolute_lane_needs_load_information(mini, cfg):
    absolute = mini["spiked_s42"]["lanes"].query("lane == 'absolute_abundance'")
    assert not absolute.available.any()


def test_absolute_lane_is_available_with_load_information(cfg, tmp_path):
    cfg = load(overrides={"output_dir": str(tmp_path / "out"), "cache_dir": str(tmp_path / "cache"),
                          "genomes": {"synthetic_length": 12000},
                          "experiment": {"samples_per_group": 3, "read_pairs": 2000,
                                         "load_information": True},
                          "analysis": {"min_total_count": 0, "min_prevalence": 0.0,
                                       "metrics": ["counts"]}})
    from magician.design import prepare_references, make_design
    from magician.truth import build
    refs = tmp_path / "refs"
    prepare_references(cfg, refs)
    make_design(cfg, "spiked_s42", refs, tmp_path / "case")
    features = tmp_path / "features.tsv"
    truth = table(tmp_path / "case/truth.tsv")[["feature_id", "length_bp"]]
    truth["retained_for_da"] = True
    truth.to_csv(features, sep="\t", index=False)
    lanes = build(cfg, "spiked_s42", "sources", tmp_path / "case/expectations.tsv",
                  tmp_path / "case/samples.tsv", features, tmp_path / "lanes.tsv",
                  zero_policies=zero_policies(cfg))
    absolute = lanes.query("lane == 'absolute_abundance'")
    assert absolute.available.all()


def test_read_fraction_and_relative_lanes_differ_only_with_unequal_lengths(mini, cfg):
    lanes = mini["spiked_s42"]["lanes"]
    relative = lanes.query("lane == 'relative_genomic_abundance'").set_index("feature_id").true_effect
    read = lanes.query("lane == 'expected_read_fraction'").set_index("feature_id").true_effect
    # The synthetic fixture has equal genome lengths, so the two scales coincide here.
    assert np.allclose(relative.to_numpy(), read.to_numpy(), atol=1e-9)


def test_derive_is_cheap_and_repeatable_from_raw_counts(cfg, tmp_path, mini):
    state = mini["spiked_s42"]
    first = table(state["matrices"] / "counts.tsv")
    tighter = load(overrides=dict(output_dir=cfg["output_dir"], cache_dir=cfg["cache_dir"],
                                  genomes=cfg["genomes"], experiment=cfg["experiment"],
                                  analysis=dict(cfg["analysis"], min_total_count=10**9)))
    derive(tighter, "sources", state["raw"], tmp_path / "tighter")
    assert table(tmp_path / "tighter/counts.tsv").equals(first)
    features = table(tmp_path / "tighter/features.tsv")
    assert not features.retained_for_da.any()
    assert read_json(tmp_path / "tighter/matrix.json")["n_retained"] == 0


def test_reserved_other_category_preserves_totals_and_is_never_tested(cfg, tmp_path, mini):
    state = mini["spiked_s42"]
    plain = pd.read_csv(state["matrices"] / "counts.tsv", sep="\t", index_col=0)
    reserved = load(overrides=dict(output_dir=cfg["output_dir"], cache_dir=cfg["cache_dir"],
                                   genomes=cfg["genomes"], experiment=cfg["experiment"],
                                   analysis=dict(cfg["analysis"], reserve_other=True,
                                                 min_total_count=10**7)))
    derive(reserved, "sources", state["raw"], tmp_path / "reserved")
    counts = pd.read_csv(tmp_path / "reserved/counts.tsv", sep="\t", index_col=0)
    assert OTHER in counts.index
    # The filtered features leave as one category and the library totals still match.
    assert len(counts) < len(plain)
    assert np.array_equal(counts.sum(axis=0).to_numpy(), plain.sum(axis=0).to_numpy())
    features = table(tmp_path / "reserved/features.tsv").set_index("feature_id")
    assert not bool(features.loc[OTHER, "retained_for_da"])
    assert bool(features.loc[OTHER, "reserved_other"])


def _write_bundle(root, case, state, samples=None, counts=None):
    """A bundle case with every part the importer requires."""
    case_dir = root / case
    case_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"case": [case], "scenario": [case.split("_s")[0]], "seed": [int(case.split("_s")[1])]}
                 ).to_csv(root / "cases.tsv", sep="\t", index=False)
    ids = list(state["counts"].index)
    columns = list(state["counts"].columns)
    groups = ["Control" if s.startswith("C") else "Treatment" for s in (samples or columns)]
    pd.DataFrame({"sample_id": samples or columns, "group": groups,
                  "read_pairs": [2000] * len(samples or columns)}).to_csv(
        case_dir / "samples.tsv", sep="\t", index=False)
    pd.DataFrame({"feature_id": ids, "length_bp": [12000] * len(ids)}).to_csv(
        case_dir / "features.tsv", sep="\t", index=False)
    (counts if counts is not None else state["counts"]).reset_index().to_csv(
        case_dir / "counts.tsv", sep="\t", index=False)
    state["truth"].reset_index().to_csv(case_dir / "truth.tsv", sep="\t", index=False)
    table(state["design"] / "expectations.tsv").to_csv(case_dir / "expectations.tsv", sep="\t", index=False)
    pd.DataFrame({"feature_id": ids, "source_id": ids, "match_status": "matched"}).to_csv(
        case_dir / "matching.tsv", sep="\t", index=False)
    return case_dir


def test_matrix_bundle_import_needs_no_reads(cfg, tmp_path, mini):
    state = mini["spiked_s42"]
    bundle = tmp_path / "bundle"
    _write_bundle(bundle, "spiked_s42", state)
    matrix_cfg = load(overrides={"output_dir": str(tmp_path / "out"), "cache_dir": str(tmp_path / "cache"),
                                 "input_mode": "matrix", "matrix_input": {"bundle": "PLACEHOLDER"},
                                 "analysis": {"min_total_count": 0, "min_prevalence": 0.0,
                                              "metrics": ["counts", "relative"]}})
    matrix_cfg["matrix_input"]["bundle"] = str(bundle)
    assert matrix.cases(matrix_cfg) == ["spiked_s42"]
    matrix.import_case(matrix_cfg, "spiked_s42", tmp_path / "out")
    imported = pd.read_csv(tmp_path / "out/matrices/spiked_s42/raw_sources/counts.tsv", sep="\t", index_col=0)
    assert imported.equals(state["counts"])
    assert (tmp_path / "out/design/spiked_s42/expectations.tsv").exists()
    matrix.import_manifest(matrix_cfg, tmp_path / "out/provenance/matrix_bundle.json")
    manifest = read_json(tmp_path / "out/provenance/matrix_bundle.json")
    assert manifest["n_cases"] == 1
    assert "spiked_s42/counts.tsv" in manifest["files"]


def test_matrix_bundle_rejects_disagreeing_parts(cfg, tmp_path, mini):
    state = mini["spiked_s42"]
    bundle = tmp_path / "bad"
    _write_bundle(bundle, "spiked_s42", state, samples=["C01"])
    bad = load(overrides={"output_dir": str(tmp_path / "out"), "cache_dir": str(tmp_path / "cache"),
                          "input_mode": "matrix", "matrix_input": {"bundle": str(bundle)}})
    with pytest.raises(ValueError, match="Control and Treatment"):
        matrix.validate_case(bad, "spiked_s42")


def test_matrix_mode_requires_a_bundle():
    with pytest.raises(ValueError, match="requires matrix_input.bundle"):
        load(overrides={"input_mode": "matrix"})
    with pytest.raises(ValueError, match="only valid when input_mode is matrix"):
        load(overrides={"matrix_input": {"bundle": "somewhere"}})


def test_matrix_bundle_requires_integer_counts(cfg, tmp_path, mini):
    state = mini["spiked_s42"]
    bundle = tmp_path / "naive"
    _write_bundle(bundle, "spiked_s42", state, counts=state["counts"].astype(float) + 0.5)
    bad = load(overrides={"output_dir": str(tmp_path / "out"), "cache_dir": str(tmp_path / "cache"),
                          "input_mode": "matrix", "matrix_input": {"bundle": str(bundle)}})
    with pytest.raises(ValueError, match="integer counts"):
        matrix.import_case(bad, "spiked_s42", tmp_path / "out")

def test_mag_truth_inherits_source_labels_through_matching(cfg, tmp_path, mini):
    from magician.truth import build
    from magician.config import zero_policies
    state = mini["spiked_s42"]
    source_id = state["truth"].index[0]
    # Two MAGs for the first source, one ambiguous MAG with no source truth.
    matching = pd.DataFrame({"feature_id": ["MAG_1", "MAG_2", "MAG_3"],
                             "source_id": [source_id] * 2 + [""],
                             "match_status": ["matched", "matched", "ambiguous"]})
    matching.to_csv(tmp_path / "matching.tsv", sep="\t", index=False)
    features = pd.DataFrame({"feature_id": ["MAG_1", "MAG_2", "MAG_3"],
                             "length_bp": [6000, 6000, 6000],
                             "retained_for_da": [True, True, True],
                             "reserved_other": [False, False, False]})
    features.to_csv(tmp_path / "features.tsv", sep="\t", index=False)
    lanes = build(cfg, "spiked_s42", "mags", state["design"] / "expectations.tsv",
                  state["design"] / "samples.tsv", tmp_path / "features.tsv",
                  tmp_path / "lanes.tsv", zero_policies=zero_policies(cfg),
                  matching_file=tmp_path / "matching.tsv")
    relative = lanes.query("lane == 'relative_genomic_abundance'").set_index("feature_id")
    # Both MAGs of one source share its label and effect; the ambiguous MAG is absent.
    assert set(relative.index) == {"MAG_1", "MAG_2"}
    assert (relative.is_da == relative.is_da.iloc[0]).all()
    assert (relative.true_effect == relative.true_effect.iloc[0]).all()
    assert bool(relative.is_da.iloc[0]) == bool(state["truth"].loc[source_id, "is_da"])
    # The MAG effect is the source lane effect in the same coordinates: length
    # correction cancels under uniform recovery, so labels and effects agree.
    source_lanes = state["lanes"].query("lane == 'relative_genomic_abundance'").set_index("feature_id")
    assert relative.true_effect.iloc[0] == pytest.approx(
        float(source_lanes.true_effect.loc[source_id]))
