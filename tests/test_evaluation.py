"""Evaluation, aggregation and recommendation behaviour on miniature real cases."""
from pathlib import Path
import math
import numpy as np
import pandas as pd
import pytest
from magician.config import load
from magician.evaluation import evaluate_case, aggregate, binomial_interval
from magician.io import save_table, table, write_json, read_json
from conftest import write_result


def _status(path, case, kind, metric, method, endpoint, status="success", **extra):
    write_json(path, dict(case=case, kind=kind, metric=metric, method=method, status=status,
                          method_family="test_family", variant="test", endpoint=endpoint,
                          effect_scale="log2_clr", reference="declared", adjustment_family="BH",
                          zero_policy=extra.get("zero_policy", "none"), native_q=True,
                          package_version="test", raw_p_available=extra.get("raw_p_available", True),
                          elapsed_seconds=1.0, message=extra.get("message", "")))


def _inputs(mini, case, tmp_path):
    """Truth, merged lanes and an identity matching table for one miniature case."""
    state = mini[case]
    ids = list(state["counts"].index)
    state["truth"].reset_index()[["feature_id", "is_da", "true_log2fc"]].to_csv(
        tmp_path / "truth.tsv", sep="\t", index=False)
    state["lanes"].to_csv(tmp_path / "lanes.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame({"feature_id": ids, "source_id": ids, "match_status": "matched"}).to_csv(
        tmp_path / "matching.tsv", sep="\t", index=False)
    return state, ids


def _discoveries(truth, n):
    return np.where(truth.is_da.to_numpy(), 1e-6, 0.9), np.where(truth.is_da.to_numpy(), 1e-5, 0.9)


def test_count_method_is_scored_on_the_read_fraction_lane(mini, cfg, tmp_path):
    case = "spiked_s42"
    state, ids = _inputs(mini, case, tmp_path)
    read_truth = state["lanes"].query("lane == 'expected_read_fraction'").set_index("feature_id")
    results = tmp_path / "results"; results.mkdir()
    p, q = _discoveries(state["truth"], len(ids))
    write_result(results / "counts.tsv", ids, endpoint="read_fraction",
                 effect_scale="log2_read_fraction", log2fc=read_truth.true_effect.to_numpy(),
                 pvalue=p, qvalue=q)
    _status(results / "counts.status.json", case, "sources", "counts", "deseq2", "read_fraction")
    evaluate_case(cfg, case, tmp_path / "truth.tsv", tmp_path / "lanes.tsv", tmp_path / "matching.tsv",
                  [results / "counts.tsv"], [results / "counts.status.json"], tmp_path / "eval")
    score = pd.read_csv(tmp_path / "eval/scores.tsv", sep="\t").iloc[0]
    assert score.endpoint == "read_fraction"
    assert score.fdr == 0
    assert score.source_recall == 1
    # Perfect effects in matching units give zero fold-change error.
    assert score.log2fc_mae == pytest.approx(0, abs=1e-9)
    assert score.n_effect_comparable == len(ids)


def test_unmatched_units_report_no_fold_change_error(mini, cfg, tmp_path):
    case = "spiked_s42"
    state, ids = _inputs(mini, case, tmp_path)
    results = tmp_path / "results"; results.mkdir()
    p, q = _discoveries(state["truth"], len(ids))
    # ZicoSeq reports a slope on transformed ratios, not a log2 fold change.
    write_result(results / "z.tsv", ids, endpoint="clr", effect_scale="sqrt_ratio_slope",
                 log2fc=[0.5] * len(ids), pvalue=p, qvalue=q)
    _status(results / "z.status.json", case, "sources", "counts", "zicoseq", "clr")
    evaluate_case(cfg, case, tmp_path / "truth.tsv", tmp_path / "lanes.tsv", tmp_path / "matching.tsv",
                  [results / "z.tsv"], [results / "z.status.json"], tmp_path / "eval")
    score = pd.read_csv(tmp_path / "eval/scores.tsv", sep="\t").iloc[0]
    assert math.isnan(score.log2fc_mae)
    assert score.n_effect_comparable == 0
    assert score.average_precision_source == "raw_p"


def test_q_only_adapter_gets_a_labelled_q_ranked_metric(mini, cfg, tmp_path):
    case = "spiked_s42"
    state, ids = _inputs(mini, case, tmp_path)
    results = tmp_path / "results"; results.mkdir()
    p, q = _discoveries(state["truth"], len(ids))
    write_result(results / "r.tsv", ids, endpoint="read_fraction", effect_scale="log2_read_fraction",
                 log2fc=[0.0] * len(ids), pvalue=[np.nan] * len(ids), qvalue=q)
    _status(results / "r.status.json", case, "sources", "counts", "future", "read_fraction",
            raw_p_available=False)
    evaluate_case(cfg, case, tmp_path / "truth.tsv", tmp_path / "lanes.tsv", tmp_path / "matching.tsv",
                  [results / "r.tsv"], [results / "r.status.json"], tmp_path / "eval")
    score = pd.read_csv(tmp_path / "eval/scores.tsv", sep="\t").iloc[0]
    assert math.isnan(score.average_precision)
    assert score.average_precision_q == pytest.approx(1.0)
    assert score.average_precision_source == "q_ranked"


def test_false_discoveries_and_recovery_accounting_are_visible(mini, cfg, tmp_path):
    case = "spiked_s42"
    state, ids = _inputs(mini, case, tmp_path)
    results = tmp_path / "results"; results.mkdir()
    write_result(results / "r.tsv", ids, endpoint="read_fraction", effect_scale="log2_read_fraction",
                 log2fc=[0.0] * len(ids), pvalue=0.001, qvalue=0.002)
    _status(results / "r.status.json", case, "sources", "counts", "deseq2", "read_fraction")
    evaluate_case(cfg, case, tmp_path / "truth.tsv", tmp_path / "lanes.tsv", tmp_path / "matching.tsv",
                  [results / "r.tsv"], [results / "r.status.json"], tmp_path / "eval")
    score = pd.read_csv(tmp_path / "eval/scores.tsv", sep="\t").iloc[0]
    assert score.fdr > 0
    assert score.n_unassigned_significant == 0
    recovery = table(tmp_path / "eval/recovery.tsv").iloc[0]
    # Every source was recovered, so only the truly changed ones count as recovered DA.
    assert recovery.n_recovered_da_sources == int(state["truth"].is_da.sum())
    assert recovery.source_da_genomes == int(state["truth"].is_da.sum())


def test_failed_and_timeout_methods_are_ineligible(tmp_path):
    cfg = load()
    rows = [dict(kind="sources", metric="counts", endpoint="read_fraction", method="deseq2",
                 method_family="count_nb", variant="default", zero_policy="none", case=scenario+"_s42",
                 scenario=scenario, seed=42, status=status, n_evaluable=5, n_reserved_other=0,
                 n_unassigned_significant=0, fdr=0, recall=1, source_recall=1, average_precision=1,
                 average_precision_q=1, any_false_discovery=0, da_elapsed_seconds=1.0)
            for scenario, status in (("null", "success"), ("spiked", "timeout"))]
    save_table(pd.DataFrame(rows), tmp_path / "scores.tsv")
    save_table(pd.DataFrame({"case": ["null_s42"]}), tmp_path / "recovery.tsv")
    aggregate(cfg, [tmp_path / "scores.tsv"], [tmp_path / "recovery.tsv"], tmp_path / "benchmark")
    ranking = table(tmp_path / "benchmark/ranking.tsv")
    assert not ranking.eligible.any()
    assert ranking.n_failed.iloc[0] == 1
    assert ranking.n_timeout.iloc[0] == 1


def test_literal_null_scenario_retained_and_null_interval_is_exact(tmp_path):
    cfg = load()
    rows = [dict(kind="sources", metric="counts", endpoint="read_fraction", method="edger",
                 method_family="count_nb", variant="tmm", zero_policy="none", case=scenario+"_s42",
                 scenario=scenario, seed=42, status="success", n_evaluable=4, n_reserved_other=0,
                 n_unassigned_significant=0, fdr=0, recall=1, source_recall=1, average_precision=1,
                 average_precision_q=1, any_false_discovery=0, da_elapsed_seconds=1.0)
            for scenario in ["null", "spiked"]]
    rows += [dict(row, method="wilcoxon_clr", variant="pseudocount_measurement",
                  method_family="clr_wilcoxon", recall=0, source_recall=0) for row in rows.copy()]
    save_table(pd.DataFrame(rows), tmp_path / "scores.tsv")
    save_table(pd.DataFrame({"case": ["null_s42", "spiked_s42"]}), tmp_path / "recovery.tsv")
    aggregate(cfg, [tmp_path / "scores.tsv"], [tmp_path / "recovery.tsv"], tmp_path / "benchmark")
    ranks = table(tmp_path / "benchmark/ranking.tsv").set_index("method")
    ranking = ranks.loc["edger"]
    assert ranking.n_null == 1
    assert ranking.null_probability_any_false == 0
    assert ranking.eligible
    assert not ranks.loc["wilcoxon_clr", "eligible"]
    low, high = binomial_interval(0, 1)
    assert low == 0 and high > 0.95


def test_variants_of_one_family_are_reported_as_one_candidate(tmp_path):
    cfg = load()
    rows = [dict(kind="sources", metric="counts", endpoint="clr", method=method,
                 method_family="dirichlet_clr", variant=variant, zero_policy="none",
                 case=scenario+"_s42", scenario=scenario, seed=42, status="success", n_evaluable=4,
                 n_reserved_other=0, n_unassigned_significant=0, fdr=0, recall=1, source_recall=1,
                 average_precision=1, average_precision_q=1, any_false_discovery=0, da_elapsed_seconds=1.0)
            for method, variant in (("aldex2", "gamma_0.5"), ("aldex2_g0", "gamma_0"))
            for scenario in ["null", "spiked"]]
    save_table(pd.DataFrame(rows), tmp_path / "scores.tsv")
    save_table(pd.DataFrame({"case": ["null_s42", "spiked_s42"]}), tmp_path / "recovery.tsv")
    aggregate(cfg, [tmp_path / "scores.tsv"], [tmp_path / "recovery.tsv"], tmp_path / "benchmark")
    recommendations = read_json(tmp_path / "benchmark/recommendations.json")
    assert len(recommendations) == 1
    # Two equally scoring settings of one Dirichlet family are one candidate, not two.
    assert recommendations[0]["candidate_families"] == ["dirichlet_clr"]
    assert recommendations[0]["best_method"] is None
    assert len(recommendations[0]["candidate_methods"]) == 2