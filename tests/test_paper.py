"""Paper plan: Tier 1 generators, calibration, variance, replication, reference."""
import numpy as np
import pandas as pd
import pytest
from magician import generators, calibration, variance, replication, reference, truth as truth_module
from magician.io import table


def _spec(**over):
    base = dict(samples_per_group=5, n_features=20, depth=20000, fraction_da=0.2,
                log2fc=1.0, direction="balanced", total_load="none",
                zero_structure="sampling", confounder="none")
    base.update(over)
    return base


def test_g1_g2_g3_write_bundle_compatible_cases(tmp_path):
    for gen in ("g1", "g2", "g3"):
        out = tmp_path / gen / "null_s1"
        rec = generators.generate_case(gen, _spec(), "null", 1, out)
        assert rec["generator"] == gen
        for name in ("counts.tsv", "samples.tsv", "features.tsv", "truth.tsv",
                     "expectations.tsv", "matching.tsv"):
            assert (out / name).exists(), f"{gen}/{name} missing"
        truth = table(out / "truth.tsv")
        assert not truth.is_da.any()
        spiked = tmp_path / gen / "spiked_s2"
        generators.generate_case(gen, _spec(), "spiked", 2, spiked)
        assert table(spiked / "truth.tsv").is_da.any()


def test_g2_resamples_distinct_measured_profiles_and_never_invents_features(tmp_path):
    spec = _spec(samples_per_group=20, n_features=300, depth=500000)
    rec = generators.generate_case("g2", spec, "null", 21, tmp_path / "big")
    assert rec["provenance"]["study"] == "AsnicarF_2021" and rec["provenance"]["n_donors"] >= 900
    profiles = table(tmp_path / "big/expectations.tsv").pivot(
        index="feature_id", columns="sample_id", values="genomic_proportion")
    assert profiles.T.drop_duplicates().shape[0] == 40
    # Measured zeros are kept: a real profile is sparse and so are its counts.
    assert 0.5 < (profiles == 0).to_numpy().mean() < 0.9
    assert rec["zero_fraction"] >= (profiles == 0).to_numpy().mean()
    pool, _ = generators._empirical_pool(300)
    drawn = profiles.T.to_numpy()
    assert all(np.abs(pool - row).max(axis=1).min() < 1e-12 for row in drawn)
    available = generators._real_profile_data()[0].shape[1]
    with pytest.raises(ValueError, match="synthetic expansion is disabled"):
        generators.generate_case("g2", _spec(n_features=available + 1), "null", 22, tmp_path / "wide")


@pytest.mark.parametrize("gen", ["g1", "g3"])
def test_parametric_generators_take_spread_and_sparsity_from_the_donor_table(gen, tmp_path):
    spec = _spec(samples_per_group=20, n_features=300, depth=500000)
    rec = generators.generate_case(gen, spec, "null", 23, tmp_path / gen)
    assert rec["provenance"]["table_sha256"] == generators._real_profile_data()[1]["table_sha256"]
    counts = table(tmp_path / gen / "counts.tsv").set_index("feature_id").to_numpy()
    assert (counts == 0).mean() > (0.5 if gen == "g3" else 0.02)
    share = counts / counts.sum(axis=0)
    assert np.log(share.mean(axis=1)[share.mean(axis=1) > 0]).std() > 1.5


def test_generators_honour_load_and_zeros(tmp_path):
    rec = generators.generate_case("g1", _spec(direction="all-up", log2fc=2.0), "spiked", 3,
                                    tmp_path / "up")
    truth = table(tmp_path / "up" / "truth.tsv")
    assert truth.is_implanted.any()
    assert (truth.loc[truth.is_implanted, "true_log2fc"] > 0).all()
    generators.generate_case("g1", _spec(total_load="x2"), "spiked", 4, tmp_path / "load")
    exp = table(tmp_path / "load" / "expectations.tsv")
    assert "absolute_copies" in exp.columns
    generators.generate_case("g1", _spec(zero_structure="prevalence-dependent"), "spiked", 5,
                              tmp_path / "zero")
    counts = table(tmp_path / "zero" / "counts.tsv")
    assert (counts.set_index("feature_id") >= 0).all().all()


@pytest.mark.parametrize("gen", generators.GENERATORS)
@pytest.mark.parametrize("direction", ["balanced", "all-up", "all-down"])
def test_implants_survive_and_truth_includes_closure(gen, direction, tmp_path):
    out = tmp_path / "case"
    rec = generators.generate_case(gen, _spec(direction=direction), "spiked", 41, out)
    truth = table(out / "truth.tsv")
    chosen = truth.is_implanted
    assert chosen.sum() == 4
    assert rec["n_implanted"] == 4
    np.testing.assert_allclose(truth.treatment_proportion.sum(), 1)
    # Design truth is the implanted set; closure shifts are recorded beside it.
    assert truth.is_da.equals(truth.is_implanted)
    if direction == "balanced":
        assert truth.is_relative_shifted.sum() == 4
        np.testing.assert_allclose(truth.loc[~chosen, "true_log2fc"], 0, atol=1e-10)
    else:
        assert truth.is_relative_shifted.sum() == len(truth)
        sign = 1 if direction == "all-up" else -1
        assert (truth.loc[chosen, "true_log2fc"] * sign > 0).all()
        assert (truth.loc[~chosen, "true_log2fc"] * sign < 0).all()
    counts = table(out / "counts.tsv").set_index("feature_id")
    samples = table(out / "samples.tsv").set_index("sample_id")
    np.testing.assert_array_equal(counts.sum().to_numpy(), samples.loc[counts.columns, "read_pairs"].to_numpy())


@pytest.mark.parametrize("gen", generators.GENERATORS)
def test_zero_effect_uses_identical_null_and_spiked_noise(gen, tmp_path):
    spec = _spec(log2fc=0)
    for scenario in ("null", "spiked"):
        generators.generate_case(gen, spec, scenario, 18, tmp_path / scenario)
    pd.testing.assert_frame_equal(table(tmp_path / "null/counts.tsv"), table(tmp_path / "spiked/counts.tsv"))
    pd.testing.assert_frame_equal(table(tmp_path / "null/expectations.tsv"), table(tmp_path / "spiked/expectations.tsv"))


@pytest.mark.parametrize("gen", generators.GENERATORS)
def test_spike_changes_effect_without_shrinking_biological_noise(gen, tmp_path):
    for scenario in ("null", "spiked"):
        generators.generate_case(gen, _spec(direction="all-up"), scenario, 33, tmp_path / scenario)
    null = table(tmp_path / "null/expectations.tsv").pivot(index="sample_id", columns="feature_id", values="genomic_proportion")
    spike = table(tmp_path / "spiked/expectations.tsv").pivot(index="sample_id", columns="feature_id", values="genomic_proportion")
    truth = table(tmp_path / "spiked/truth.tsv").set_index("feature_id").reindex(null.columns)
    target_clr = truth.true_log2fc - truth.true_log2fc.mean()
    control = null.index.str.startswith("C")
    # A measured zero stays a zero; present features move by the designed effect.
    assert ((spike == 0) == (null == 0)).all().all()
    np.testing.assert_allclose(spike.loc[control], null.loc[control], atol=1e-15)
    for sample in null.index[~control]:
        present = null.loc[sample] > 0
        change = np.log2(spike.loc[sample, present] / null.loc[sample, present])
        offset = truth.true_log2fc[present]
        np.testing.assert_allclose(change - change.mean(), offset - offset.mean(), atol=1e-9)


@pytest.mark.parametrize("gen", generators.GENERATORS)
def test_unsupported_confounder_is_rejected(gen, tmp_path):
    with pytest.raises(ValueError, match="Batch confounding is excluded"):
        generators.generate_case(gen, _spec(confounder="batch_0.6"), "spiked", 6, tmp_path / "conf")
    assert not (tmp_path / "conf").exists()


def test_core_and_fractional_grids_cover_axes():
    core = generators.core_cells()
    assert len(core) == 24 and all(c["fraction_da"] > 0 for c in core)
    assert {c["generator"] for c in core} == {"g1", "g2", "g3"}
    frac = generators.fractional_grid()
    assert len(frac) > 10 and not any(c in core for c in frac)
    assert all(c["confounder"] == "none" for c in core + frac)
    for key, levels in generators.TIER1_AXES.items():
        assert {c[key] for c in core + frac} == set(levels)


def test_null_cases_are_generated_once_per_null_condition():
    plan = generators.tier1_plan()
    nulls = [e for e in plan if e["scenario"] == "null"]
    spiked = [e for e in plan if e["scenario"] == "spiked"]
    key = lambda cell: tuple(cell[f] for f in generators.NULL_FACTORS)
    assert len({key(e["cell"]) for e in nulls}) == len(nulls)
    assert {key(e["cell"]) for e in spiked} == {key(e["cell"]) for e in nulls}
    assert all(e["cell"]["fraction_da"] == 0 and e["cell"]["direction"] == "none" for e in nulls)
    assert sum(e["replicates"] for e in nulls if e["design"].startswith("core")) == 6 * 200
    assert sum(e["replicates"] for e in spiked if e["design"].startswith("core")) == 24 * 100
    replicates = {key(e["cell"]): e["replicates"] for e in nulls}
    # Every core spiked cell, and every fractional cell on a core null condition, has 200 nulls.
    assert all(replicates[key(e["cell"])] == 200 for e in spiked if e["design"].startswith("core"))


def test_bundle_preserves_scenario_factors_and_previous_data(tmp_path):
    out = tmp_path / "bundle"
    plan = generators.design_plan([dict(_spec(), generator="g1"), dict(_spec(fraction_da=0.4), generator="g1")],
                                  1, 1, "core")
    generators.build_bundle(plan, out)
    listing = pd.read_csv(out / "cases.tsv", sep="\t", keep_default_na=False).set_index("case")
    assert listing.scenario.tolist() == ["null", "spiked", "spiked"]
    assert set(listing.generator_version) == {generators.GENERATOR_VERSION}
    assert listing.design.tolist() == ["core_n000", "core_c000", "core_c001"]
    assert listing.direction.tolist() == ["none", "balanced", "balanced"]
    assert listing.fraction_da.tolist() == [0.0, 0.2, 0.4]
    assert not table(out / "null_s1000/truth.tsv").is_da.any()
    before = (out / "null_s1000/counts.tsv").read_bytes()
    with pytest.raises(ValueError, match="not empty"):
        generators.build_bundle(plan, out)
    assert (out / "null_s1000/counts.tsv").read_bytes() == before


def test_clr_truth_accounts_for_a_changed_geometric_reference(cfg):
    # Relative abundance changes only for a/b, but CLR changes for all four because
    # the geometric mean falls. This is an endpoint difference, not four implants.
    ids = ["a", "b", "c", "d"]
    control, treatment = np.full(4, 0.25), np.array([0.4, 0.1, 0.25, 0.25])
    records, samples = [], []
    for name, target in [("C01", control), ("C02", control), ("T01", treatment), ("T02", treatment)]:
        samples.append(dict(sample_id=name, read_pairs=1000))
        records.extend(dict(sample_id=name, feature_id=fid, genomic_proportion=p,
                            target_proportion=p, read_proportion=p, expected_count=p * 1000,
                            library_read_pairs=1000, length_bp=20000, absolute_copies=p)
                       for fid, p in zip(ids, target))
    expected = pd.DataFrame(records)
    lengths = pd.Series(20000, index=ids)
    retained = pd.Series(True, index=ids)
    lane = truth_module.clr_lane(cfg, "spiked_s1", "sources", expected, pd.DataFrame(samples),
                                 lengths, retained, "counts",
                                 {"name": "none", "zero_handling": "pseudocount", "pseudocount": 0}, False)
    assert lane.is_da.all()
    ratio = treatment / control
    theoretical = np.log2(ratio / np.prod(ratio) ** 0.25)
    np.testing.assert_allclose(lane.true_effect, theoretical, atol=1e-9)
    expected["target_proportion"] = 0.25
    null = truth_module.clr_lane(cfg, "null_s1", "sources", expected, pd.DataFrame(samples),
                                lengths, retained, "counts",
                                {"name": "none", "zero_handling": "pseudocount", "pseudocount": 0}, False)
    assert not null.is_da.any()


def test_calibration_pvalue_level():
    rng = np.random.default_rng(0)
    null_p = rng.uniform(size=1000)
    assert 0.03 <= calibration.fpr_at_alpha(null_p) <= 0.07
    assert calibration.ks_uniformity(null_p) < 0.05
    qq = calibration.qq_data(null_p)
    assert len(qq["expected"]) > 0
    # Cluster bootstrap keeps methods paired (resamples whole cases).
    frame = pd.DataFrame({"case": ["a", "b"] * 4, "method": ["m1", "m1", "m2", "m2"] * 2,
                          "fdr": [0.1, 0.2, 0.3, 0.4] * 2})
    ci = calibration.cluster_bootstrap_ci(frame, "fdr", ["method"])
    assert len(ci) == 2


def test_variance_contrasts_ranks():
    scores = pd.DataFrame({
        "case": ["a", "b", "c"] * 4, "method": (["m1"] * 3 + ["m2"] * 3) * 2,
        "fdr": [0.02, 0.03, 0.02, 0.20, 0.25, 0.22] * 2,
        "source_recall": [0.8, 0.7, 0.75, 0.9, 0.85, 0.9] * 2,
        "generator": ["g1"] * 12})
    var = variance.variance_decomposition(scores, ["method", "generator"], "fdr")
    assert var["method"] > var.get("generator", 0)
    ask = variance.pairwise_contrasts(scores, "method", "fdr", draws=200)
    assert len(ask) == 2
    ranks = variance.rank_uncertainty(scores, ("method",), "source_recall", draws=200)
    assert len(ranks) == 2


def test_rank_uncertainty_does_not_break_ties_by_method_name():
    frame = pd.DataFrame({"case": ["a", "b", "c"] * 2,
                          "method": ["first"] * 3 + ["second"] * 3,
                          "source_recall": [0.8] * 6})
    ranks = variance.rank_uncertainty(frame, draws=200)
    assert ranks[("first",)] == ranks[("second",)] == (1.5, 1.5, 1.5)


def test_replication_splits_implants_rate():
    splits = replication.random_splits([f"S{i:02d}" for i in range(10)], 4, seed=0)
    assert len(splits) == 4
    counts = pd.DataFrame({"A": [100, 50], "B": [100, 50]}, index=["g1", "g2"])
    implanted, truth = replication.implant_reads(counts, ["g1"], 2.0, seed=0)
    assert truth["g1"] and not truth["g2"]
    assert implanted.loc["g1", "A"] != counts.loc["g1", "A"]
    r = replication.replication_rate(["a", "b", "c"], ["b", "c", "d"])
    assert 0 <= r["jaccard"] <= 1


def test_reference_select_filter_catalogue(tmp_path):
    refs = pd.DataFrame({"genome_id": [f"g{i}" for i in range(10)],
                         "length_bp": [10000] * 10})
    refs.to_csv(tmp_path / "refs.tsv", sep="\t", index=False)
    frame = reference.select_present(tmp_path / "refs.tsv", 0.7, 7, tmp_path / "present.tsv")
    assert frame.present.sum() == 7
    reference.reference_matching(tmp_path / "present.tsv", tmp_path / "matching.tsv")
    assert len(table(tmp_path / "matching.tsv")) == 7
    reference.build_catalogue(tmp_path / "refs.tsv", tmp_path / "present.tsv",
                              tmp_path / "feats.tsv", tmp_path / "matching2.tsv")
    assert len(table(tmp_path / "feats.tsv")) == 7
