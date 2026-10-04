"""Truth-aware evaluation; missing MAGs, absent lanes and failed methods stay explicit.

A method is scored against the truth lane for the endpoint it declares. Fold-change
error is only computed when the reported units and the truth coordinates agree.
Native adjusted values are never re-adjusted; where a package reports no raw
p-values, average precision over raw p-values is reported as unavailable and a
separate q-ranked average precision is recorded instead.
"""
from pathlib import Path
import math
import warnings
import numpy as np
import pandas as pd
from .io import table, save_table, read_json, write_json
from .truth import ENDPOINT_LANE
from .matrices import OTHER


def average_precision(labels, scores):
    """Threshold-grouped AP: ties do not depend on feature ordering."""
    labels, scores = np.asarray(labels, bool), np.asarray(scores, float)
    if not labels.sum():
        return math.nan
    order = np.argsort(-scores, kind="stable")
    labels, scores = labels[order], scores[order]
    ends = np.r_[np.flatnonzero(np.diff(scores)) + 1, len(scores)]
    tp = np.cumsum(labels)[ends - 1]
    recall = tp / labels.sum()
    precision = tp / ends
    return float(np.sum(np.diff(np.r_[0, recall]) * precision))


def confusion(labels, discoveries):
    labels, discoveries = np.asarray(labels, bool), np.asarray(discoveries, bool)
    tp = int((labels & discoveries).sum())
    fp = int((~labels & discoveries).sum())
    fn = int((labels & ~discoveries).sum())
    tn = int((~labels & ~discoveries).sum())
    return dict(tp=tp, fp=fp, fn=fn, tn=tn,
                fdr=fp / max(tp + fp, 1),
                recall=tp / (tp + fn) if tp + fn else math.nan,
                precision=tp / (tp + fp) if tp + fp else math.nan,
                false_positive_rate=fp / (fp + tn) if fp + tn else math.nan,
                any_false_discovery=int(fp > 0))


def binomial_interval(successes, trials, level=0.95):
    """Exact Clopper-Pearson interval; no SciPy dependency and no normal approximation."""
    if trials <= 0:
        return math.nan, math.nan
    tail = (1 - level) / 2
    lower = 0.0 if successes == 0 else _binomial_quantile(successes, trials, tail, upper=False)
    upper = 1.0 if successes == trials else _binomial_quantile(successes, trials, 1 - tail, upper=True)
    return lower, upper


def _binomial_quantile(k, n, target, upper):
    low, high = 0.0, 1.0
    for _ in range(80):
        mid = (low + high) / 2
        cdf = sum(math.comb(n, i) * mid**i * (1 - mid)**(n - i)
                  for i in range(n + 1) if (i <= k if upper else i >= k))
        if cdf > target:
            high = mid
        else:
            low = mid
    return (low + high) / 2


def select_lane(lanes, endpoint, metric, zero_policy, reserve_other, kind):
    """The one truth lane matching a method's declared endpoint and input geometry."""
    lane = ENDPOINT_LANE.get(endpoint)
    if lane is None:
        raise ValueError(f"Result rows declare an unknown endpoint: {endpoint}")
    subset = lanes.loc[lanes.lane.eq(lane) & (lanes.kind == kind)]
    if not len(subset):
        raise ValueError(f"Truth lane {lane} is missing for endpoint {endpoint}")
    if lane == "clr_log_ratio":
        subset = subset.loc[subset.input_metric.eq(metric) & subset.zero_policy.eq(zero_policy)
                            & subset.reserve_other.astype(str).str.lower().eq(str(bool(reserve_other)).lower())]
    if subset.empty:
        raise ValueError(f"No {lane} truth lane matches metric={metric}, zero_policy={zero_policy}, "
                         f"reserve_other={reserve_other}")
    return subset.drop_duplicates("feature_id")


def evaluate_case(cfg, case, truth_file, lanes_file, matching_file, result_files, status_files, output):
    truth = table(truth_file).set_index("feature_id")
    truth.is_da = truth.is_da.astype(str).str.lower().eq("true")
    lanes = table(lanes_file)
    matching = table(matching_file)
    matched = matching.loc[matching.match_status == "matched"]
    known = set(truth.index)
    if set(matched.source_id) - known:
        raise ValueError("MAG assignments refer to unknown source genomes")
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    summaries, details, statuses = [], [], []
    total_positive = int(truth.is_da.sum())
    recovery = dict(case=case, source_genomes=len(truth), source_da_genomes=total_positive,
                    n_mags=len(matching), n_matched_mags=len(matched),
                    n_recovered_sources=matched.source_id.nunique(),
                    n_recovered_da_sources=int(truth.loc[list(set(matched.source_id)), "is_da"].sum()),
                    n_ambiguous=int(matching.match_status.eq("ambiguous").sum()),
                    n_unassigned=int(matching.match_status.isin(["unassigned", "low_coverage"]).sum()),
                    n_multiply_matched_sources=int((matched.source_id.value_counts() > 1).sum()))
    for result_file, status_file in zip(result_files, status_files):
        status = read_json(status_file)
        statuses.append(status)
        frame = table(result_file)
        if not frame.feature_id.is_unique:
            raise ValueError(f"Duplicate DA result IDs: {result_file}")
        kind = status["kind"]
        if kind == "sources":
            unexpected = set(frame.feature_id) - known - {OTHER}
            if unexpected:
                raise ValueError(f"Source DA results contain unknown feature IDs: {sorted(unexpected)[:3]}")
            frame["source_id"] = frame.feature_id
            frame["match_status"] = "matched"
        else:
            if set(frame.feature_id) != set(matching.feature_id):
                raise ValueError("MAG DA result IDs disagree with catalogue")
            frame = frame.merge(matching[["feature_id", "source_id", "match_status"]], on="feature_id", validate="one_to_one")
        endpoint = status["endpoint"]
        lane = select_lane(lanes, endpoint, status["metric"],
                           status.get("truth_lane_policy", status["zero_policy"]),
                           bool(status.get("reserve_other", False)), kind)
        frame = frame.merge(lane[["feature_id", "is_da", "true_effect", "available", "in_tested_set",
                                  "effect_scale", "reference"]].rename(
                                  columns={"is_da": "truth_is_da", "true_effect": "truth_effect",
                                           "available": "truth_available", "in_tested_set": "in_tested_set",
                                           "effect_scale": "truth_effect_scale", "reference": "truth_reference"}),
                          on="feature_id", how="left", validate="one_to_one")
        # An unrecovered or ambiguous bin is a recovery failure, not biological absence.
        frame["is_da"] = frame.truth_is_da
        for key in ("pvalue", "qvalue", "log2fc", "truth_effect"):
            frame[key] = pd.to_numeric(frame[key], errors="coerce")
        frame["tested"] = frame.tested.astype(str).str.lower().eq("true")
        frame["significant"] = frame.tested & frame.qvalue.le(cfg["analysis"]["alpha"])
        frame["evaluable"] = (frame.match_status.eq("matched") & frame.in_tested_set.fillna(False).astype(bool)
                              & frame.truth_available.fillna(False).astype(bool) & frame.tested)
        frame["case"] = case
        frame["method"] = status["method"]
        frame["method_family"] = status["method_family"]
        frame["variant"] = status["variant"]
        frame["method_status"] = status["status"]
        frame["metric"] = status["metric"]
        frame["effect_units_match"] = frame.effect_scale.eq(frame.truth_effect_scale)
        details.append(frame)
        evaluated = frame.loc[frame.evaluable]
        labels = evaluated.is_da.to_numpy(dtype=bool)
        discovery = evaluated.significant.to_numpy(dtype=bool)
        score = confusion(labels, discovery)
        raw_p = status.get("raw_p_available", False)
        pvalues = evaluated.pvalue.fillna(1).to_numpy(dtype=float)
        qvalues = evaluated.qvalue.fillna(1).to_numpy(dtype=float)
        score["average_precision"] = average_precision(labels, -pvalues) if len(labels) and raw_p else math.nan
        score["average_precision_q"] = average_precision(labels, -qvalues) if len(labels) else math.nan
        score["average_precision_source"] = "raw_p" if raw_p else ("q_ranked" if len(labels) else "unavailable")
        true_sources_called = set(evaluated.loc[evaluated.significant & evaluated.is_da.astype(bool), "source_id"])
        false_sources_called = set(evaluated.loc[evaluated.significant & ~evaluated.is_da.astype(bool), "source_id"])
        score["source_recall"] = len(true_sources_called) / total_positive if total_positive else math.nan
        score["source_fdr"] = len(false_sources_called) / max(len(true_sources_called | false_sources_called), 1)
        directional = evaluated.loc[evaluated.significant & evaluated.is_da.astype(bool) & evaluated.log2fc.notna()]
        score["direction_accuracy"] = float((np.sign(directional.log2fc) == np.sign(directional.truth_effect)).mean()) if len(directional) else math.nan
        # Fold-change error only where reported units and truth coordinates agree.
        comparable = evaluated.loc[evaluated.effect_units_match & evaluated.log2fc.notna()
                                   & evaluated.truth_effect.notna()]
        score["log2fc_mae"] = float((comparable.log2fc - comparable.truth_effect).abs().mean()) if len(comparable) else math.nan
        score["n_effect_comparable"] = int(len(comparable))
        score.update(case=case, scenario=case.split("_s")[0], seed=int(case.split("_s")[1]),
                     kind=kind, metric=status["metric"], method=status["method"],
                     method_family=status["method_family"], variant=status["variant"],
                     endpoint=endpoint, effect_scale=status["effect_scale"], reference=status["reference"],
                     adjustment_family=status["adjustment_family"], zero_policy=status["zero_policy"],
                     status=status["status"], method_message=status.get("message", ""),
                     package_version=status.get("package_version", "unknown"),
                     da_elapsed_seconds=status.get("elapsed_seconds", math.nan),
                     n_features=len(frame), n_evaluable=len(evaluated), n_tested=int(frame.tested.sum()),
                     n_reserved_other=int((frame.fit_status == "reserved_category").sum()),
                     n_unassigned_significant=int((~frame.evaluable & frame.significant).sum()))
        if status["status"] != "success" or not len(evaluated):
            for key in ("fdr", "recall", "precision", "false_positive_rate", "any_false_discovery",
                        "average_precision", "average_precision_q", "source_recall", "source_fdr",
                        "direction_accuracy", "log2fc_mae"):
                score[key] = math.nan
        summaries.append(score)
    save_table(pd.DataFrame(summaries), out / "scores.tsv")
    save_table(pd.concat(details, ignore_index=True), out / "feature_evaluation.tsv")
    save_table(pd.DataFrame([recovery]), out / "recovery.tsv")
    write_json(out / "method_statuses.json", statuses)


def paired_bootstrap(frame, statistic, group_columns, seed=214, draws=2000):
    """Resample whole independent seeds, keeping every method's comparison paired."""
    if frame.empty:
        return {}
    cases = sorted(frame.case.unique())
    if len(cases) < 2:
        return {}
    rng = np.random.default_rng(seed)
    picked = rng.integers(0, len(cases), (draws, len(cases)))
    out = {}
    for values, group in frame.groupby(group_columns, sort=True):
        means = []
        for case in cases:
            values_for_case = group.loc[group.case.eq(case), statistic].dropna().to_numpy()
            means.append(values_for_case.mean() if values_for_case.size else np.nan)
        means = np.asarray(means, dtype=float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            boot = np.nanmean(means[picked], axis=1)
        boot = boot[np.isfinite(boot)]
        if not boot.size:
            continue
        out[tuple(values)] = (float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975)))
    return out


def aggregate(cfg, score_files, recovery_files, output, settings_rows=None):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    # Pandas otherwise treats the literal scenario name 'null' as missing data.
    scores = pd.concat([pd.read_csv(p, sep="\t", keep_default_na=False, na_values=["NA"])
                        for p in score_files], ignore_index=True)
    recovery = pd.concat([pd.read_csv(p, sep="\t", keep_default_na=False, na_values=["NA"])
                          for p in recovery_files], ignore_index=True)
    save_table(scores, out / "scores.tsv")
    save_table(recovery, out / "recovery.tsv")
    declared = {row["method_id"]: row for row in (settings_rows or [])}
    grouping = ["kind", "metric", "endpoint", "method"]
    rows = []
    for values, group in scores.groupby(grouping, sort=True):
        successful = group.loc[group.status.eq("success") & group.n_evaluable.gt(0)]
        spike = successful.loc[successful.scenario.eq("spiked")]
        null = successful.loc[successful.scenario.eq("null")]
        first = lambda name: group[name].iloc[0] if name in group else ""
        family = first("method_family")
        row = dict(zip(grouping, values), method_family=family, variant=first("variant"),
                   effect_scale=first("effect_scale"), reference=first("reference"),
                   adjustment_family=first("adjustment_family"), zero_policy=first("zero_policy"),
                   n_runs=len(group),
                   n_success=len(successful), n_failed=int(group.status.isin(["failed", "unavailable", "timeout"]).sum()),
                   n_empty=int(group.status.eq("empty").sum()), n_timeout=int(group.status.eq("timeout").sum()),
                   n_spiked=len(spike), n_null=len(null),
                   fdr=float(spike.fdr.mean()), recall=float(spike.recall.mean()),
                   source_recall=float(spike.source_recall.mean()),
                   average_precision=float(spike.average_precision.mean()),
                   average_precision_q=float(spike.average_precision_q.mean()) if "average_precision_q" in spike else math.nan,
                   log2fc_mae=float(spike.log2fc_mae.mean()) if "log2fc_mae" in spike else math.nan,
                   n_spiked_fdr_passes=int((spike.fdr <= cfg["analysis"]["alpha"]).sum()),
                   mean_da_seconds=float(successful.da_elapsed_seconds.mean()) if "da_elapsed_seconds" in successful else math.nan,
                   unassigned_significant=int(group.n_unassigned_significant.sum()),
                   reserved_other_features=int(group.n_reserved_other.max()) if "n_reserved_other" in group else 0)
        for name in ("fdr", "source_recall", "average_precision", "average_precision_q"):
            vals = spike[name].dropna().to_numpy()
            if len(vals) >= 2:
                rng = np.random.default_rng(214)
                boot = rng.choice(vals, (2000, len(vals)), replace=True).mean(axis=1)
                row[f"{name}_ci_low"], row[f"{name}_ci_high"] = np.quantile(boot, [0.025, 0.975])
            else:
                row[f"{name}_ci_low"] = row[f"{name}_ci_high"] = math.nan
        intervals = paired_bootstrap(scores, "fdr", grouping)
        low, high = intervals.get(tuple(values), (math.nan, math.nan))
        row["fdr_paired_ci_low"], row["fdr_paired_ci_high"] = low, high
        # Exact binomial interval for the null rate: the probability of any false
        # discovery under the global null, from independent null runs only.
        n_null_runs = len(null)
        false_null = int(null.any_false_discovery.fillna(0).sum())
        row["null_any_false_discovery"] = false_null
        row["null_probability_any_false"] = false_null / n_null_runs if n_null_runs else math.nan
        low, high = binomial_interval(false_null, n_null_runs)
        row["null_probability_ci_low"], row["null_probability_ci_high"] = low, high
        # A failed method cannot win by omitting difficult runs.
        row["eligible"] = (len(successful) == len(group) and len(spike) > 0 and n_null_runs > 0
                           and row["fdr"] <= cfg["analysis"]["alpha"]
                           and row["null_probability_any_false"] <= cfg["analysis"]["alpha"]
                           and row["source_recall"] > 0
                           and not row["unassigned_significant"])
        rows.append(row)
    ranking = pd.DataFrame(rows).sort_values(["kind", "metric", "endpoint", "eligible",
                                             "source_recall", "average_precision_q", "fdr"],
                                            ascending=[True, True, True, False, False, False, True],
                                            na_position="last")
    ranking["rank"] = ranking.groupby(["kind", "metric", "endpoint"]).cumcount() + 1
    save_table(ranking, out / "ranking.tsv")
    if settings_rows:
        settings = pd.DataFrame(settings_rows)
        settings.to_csv(out / "method_settings.tsv", sep="\t", index=False, na_rep="NA")
    recommendations = []
    for values, group in ranking.groupby(["kind", "metric", "endpoint"]):
        eligible = group.loc[group.eligible]
        candidates = []
        if len(eligible):
            best = eligible.iloc[0]
            tied = eligible.loc[np.isclose(eligible.source_recall, best.source_recall)
                                & np.isclose(eligible.average_precision_q, best.average_precision_q)
                                & np.isclose(eligible.fdr, best.fdr)]
            # Variants of one family are one candidate, never several methods.
            candidates = sorted({f"{method} [{variant}]" for method, variant
                                 in zip(tied.method, tied.variant)})
            families = sorted(set(tied.method_family))
            if len(families) == 1:
                conclusion = "eligible candidate"
            elif len(families) > 1:
                conclusion = "equally scoring eligible candidates"
            else:
                conclusion = "eligible candidate"
            best_method = candidates[0] if len(candidates) == 1 and len(families) == 1 else None
        else:
            conclusion, best_method, families = "no method meets eligibility criteria", None, []
        enough = len(cfg["experiment"]["seeds"]) >= 5 and cfg["genomes"]["manifest"] is not None
        recommendations.append(dict(kind=values[0], metric=values[1], endpoint=values[2],
                                    best_method=best_method, candidate_methods=candidates,
                                    candidate_families=families,
                                    evidence="repeated_simulation" if enough else "smoke_only",
                                    conclusion=conclusion))
    write_json(out / "recommendations.json", recommendations)