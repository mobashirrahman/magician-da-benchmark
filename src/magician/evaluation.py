"""Truth-aware evaluation; missing MAGs, absent lanes and failed methods stay explicit.

Discoveries are scored against the design truth: the features whose abundance the
design changed. It is the same set for every method, so a method cannot pass by
declaring an endpoint on which nothing is null. The truth lane for the endpoint a
method declares is scored beside it (`endpoint_*`): a discovery outside the
design truth that is a real shift on the method's own scale is counted as
compositional spillover, and endpoint FDR is unavailable when that scale has no
null feature. Fold-change error is only computed when the reported units and the
truth coordinates agree.
Native adjusted values are never re-adjusted; where a package reports no raw
p-values, average precision over raw p-values is reported as unavailable and a
separate q-ranked average precision is recorded instead.
"""
from pathlib import Path
from functools import lru_cache
import math
import warnings
import numpy as np
import pandas as pd
from .io import table, save_table, read_json, write_json
from .truth import ENDPOINT_LANE
from .matrices import OTHER

CONDITION_COLUMNS = ("generator_version", "generator", "samples_per_group", "n_features",
                     "depth", "fraction_da", "log2fc", "direction", "total_load",
                     "zero_structure", "confounder")
# Null data depend on these factors only, so null runs are shared by every spiked
# condition that agrees on them.
NULL_CONDITION_COLUMNS = ("generator_version", "generator", "samples_per_group", "n_features",
                          "depth", "zero_structure", "confounder")
NULL_FPR_MAX, NULL_FPR_CONSERVATIVE = 0.07, 0.03


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
    # Use log probabilities: integer binomial coefficients overflow for thousands
    # of null replicates. The upper bound inverts a decreasing lower-tail CDF.
    indices = np.arange(n + 1, dtype=float)
    selected = indices <= k if upper else indices >= k
    indices = indices[selected]
    coefficients = _binomial_coefficients(n)[selected]
    if upper:
        target = 1 - target
    low, high = 0.0, 1.0
    for _ in range(80):
        mid = (low + high) / 2
        cdf = float(np.exp(coefficients + indices * math.log(mid)
                           + (n - indices) * math.log1p(-mid)).sum())
        if (cdf > target) != upper:
            high = mid
        else:
            low = mid
    return (low + high) / 2


@lru_cache(maxsize=64)
def _binomial_coefficients(n):
    return np.asarray([math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
                       for i in range(n + 1)])


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


def evaluate_case(cfg, case, truth_file, lanes_file, matching_file, result_files, status_files, output,
                  reference_matching_file=None):
    truth = table(truth_file).set_index("feature_id")
    # Design truth is the implanted set when a generator separates it from closure shifts.
    truth["is_da"] = truth["is_implanted" if "is_implanted" in truth else "is_da"].astype(str).str.lower().eq("true")
    lanes = table(lanes_file)
    matching = table(matching_file)
    matched = matching.loc[matching.match_status == "matched"]
    known = set(truth.index)
    if set(matched.source_id) - known:
        raise ValueError("MAG assignments refer to unknown source genomes")
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    summaries, details, statuses = [], [], []
    design_path = Path(truth_file).with_name("design.json")
    design_record = read_json(design_path) if design_path.exists() else {}
    condition = {k: v for k, v in design_record.get("cell", {}).items() if k in CONDITION_COLUMNS}
    for key in ("generator", "generator_version"):
        if design_record.get(key):
            condition[key] = design_record[key]
    total_positive = int(truth.is_da.sum())
    recovery = dict(case=case, source_genomes=len(truth), source_da_genomes=total_positive,
                    n_mags=len(matching), n_matched_mags=len(matched),
                    n_recovered_sources=matched.source_id.nunique(),
                    n_recovered_da_sources=int(truth.loc[list(set(matched.source_id)), "is_da"].sum()) if len(matched) else 0,
                    n_ambiguous=int(matching.match_status.eq("ambiguous").sum()),
                    n_unassigned=int(matching.match_status.isin(["unassigned", "low_coverage"]).sum()),
                    n_multiply_matched_sources=int((matched.source_id.value_counts() > 1).sum()) if len(matched) else 0)
    # Tier 2 MAG-quality covariates: share of DA genomes not recovered (MNAR signal),
    # plus reference-table recovery when the incomplete-DB table is enabled.
    recovery["share_da_not_recovered"] = (
        (total_positive - recovery["n_recovered_da_sources"]) / total_positive if total_positive else math.nan)
    if reference_matching_file is not None and Path(reference_matching_file).exists():
        ref_match = table(reference_matching_file)
        ref_matched = ref_match.loc[ref_match.match_status == "matched"]
        ref_ids = set(ref_matched.source_id.tolist()) | set(ref_matched.feature_id.tolist())
        ref_ids = {g for g in ref_ids if g in known}
        recovery["n_reference_features"] = int(len(ref_match))
        recovery["n_reference_recovered_sources"] = int(len(ref_ids))
        recovery["n_reference_recovered_da"] = int(truth.loc[list(ref_ids), "is_da"].sum()) if ref_ids else 0
        recovery["share_da_missing_from_reference"] = (
            (total_positive - recovery["n_reference_recovered_da"]) / total_positive if total_positive else math.nan)
    # Per-source recovery vs abundance/completeness for the MAG-cost figure
    # (power/FDR as a function of true abundance and coverage).
    try:
        match_detail = matching.loc[matching.match_status == "matched", ["source_id", "source_completeness"]]
        abund = truth[["is_da", "baseline_proportion"]].copy()
        abund["recovered"] = abund.index.isin(set(matched.source_id))
        abund["source_completeness"] = abund.index.map(
            match_detail.drop_duplicates("source_id").set_index("source_id").source_completeness)
        abund.index.name = "source_id"
        abund.reset_index().to_csv(out / "mag_quality.tsv", sep="\t", index=False, na_rep="NA")
    except Exception:
        pass
    for result_file, status_file in zip(result_files, status_files):
        status = read_json(status_file)
        statuses.append(status)
        frame = table(result_file)
        if not frame.feature_id.is_unique:
            raise ValueError(f"Duplicate DA result IDs: {result_file}")
        kind = status["kind"]
        if kind in ("sources", "reference"):
            unexpected = set(frame.feature_id) - known - {OTHER}
            if unexpected:
                raise ValueError(f"Source DA results contain unknown feature IDs: {sorted(unexpected)[:3]}")
            frame["source_id"] = frame.feature_id
            frame["match_status"] = "matched"
        else:
            if set(frame.feature_id) != set(matching.feature_id):
                raise ValueError("MAG DA result IDs disagree with catalogue")
            frame = frame.merge(matching[["feature_id", "source_id", "match_status"]], on="feature_id", validate="one_to_one")
        declared = status["endpoint"]
        lane = select_lane(lanes, declared, status["metric"],
                           status.get("truth_lane_policy", status["zero_policy"]),
                           bool(status.get("reserve_other", False)), kind)
        frame = frame.merge(lane[["feature_id", "is_da", "true_effect", "available", "in_tested_set",
                                  "effect_scale", "reference"]].rename(
                                  columns={"is_da": "truth_is_da", "true_effect": "truth_effect",
                                           "available": "truth_available", "in_tested_set": "in_tested_set",
                                           "effect_scale": "truth_effect_scale", "reference": "truth_reference"}),
                          on="feature_id", how="left", validate="one_to_one")
        # An unrecovered or ambiguous bin is a recovery failure, not biological absence.
        frame["endpoint_is_da"] = frame.truth_is_da.astype(str).str.lower().eq("true")
        frame["is_da"] = frame.source_id.map(truth.is_da).fillna(False).astype(bool)
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
        own = evaluated.endpoint_is_da.to_numpy(dtype=bool)
        endpoint = confusion(own, discovery)
        score["endpoint_n_positive"], score["endpoint_n_null"] = int(own.sum()), int((~own).sum())
        # Without a null feature on the method's own scale no discovery can be false there.
        score["endpoint_fdr"] = endpoint["fdr"] if (~own).any() else math.nan
        score["endpoint_recall"] = endpoint["recall"]
        score["n_discoveries"] = int(discovery.sum())
        score["n_spillover_discoveries"] = int((discovery & ~labels & own).sum())
        score["spillover_share"] = score["n_spillover_discoveries"] / max(score["n_discoveries"], 1)
        raw_p = status.get("raw_p_available", False)
        pvalues = evaluated.pvalue.fillna(1).to_numpy(dtype=float)
        qvalues = evaluated.qvalue.fillna(1).to_numpy(dtype=float)
        score["average_precision"] = average_precision(labels, -pvalues) if len(labels) and raw_p else math.nan
        score["average_precision_q"] = average_precision(labels, -qvalues) if len(labels) else math.nan
        score["average_precision_source"] = "raw_p" if raw_p else ("q_ranked" if len(labels) else "unavailable")
        true_sources_called = set(evaluated.loc[evaluated.significant & evaluated.is_da.astype(bool), "source_id"])
        false_sources_called = set(evaluated.loc[evaluated.significant & ~evaluated.is_da.astype(bool), "source_id"])
        # Every design-positive source counts, including those filtered out or not recovered.
        source_positive = total_positive
        score["n_source_positive"] = source_positive
        score["source_recall"] = len(true_sources_called) / source_positive if source_positive else math.nan
        null_p = pd.to_numeric(evaluated.loc[~evaluated.is_da.astype(bool), "pvalue"], errors="coerce").dropna()
        score["n_null_pvalues"] = int(len(null_p))
        score["null_pvalue_fpr"] = float((null_p <= cfg["analysis"]["alpha"]).mean()) if raw_p and len(null_p) else math.nan
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
                     endpoint=declared, effect_scale=status["effect_scale"], reference=status["reference"],
                     adjustment_family=status["adjustment_family"], zero_policy=status["zero_policy"],
                     status=status["status"], method_message=status.get("message", ""),
                     package_version=status.get("package_version", "unknown"),
                     da_elapsed_seconds=status.get("elapsed_seconds", math.nan),
                     n_features=len(frame), n_evaluable=len(evaluated), n_tested=int(frame.tested.sum()),
                     n_reserved_other=int((frame.fit_status == "reserved_category").sum()),
                     n_unassigned_significant=int((~frame.evaluable & frame.significant).sum()))
        score.update(condition)
        if status["status"] != "success" or not len(evaluated):
            for key in ("fdr", "recall", "precision", "false_positive_rate", "any_false_discovery",
                        "average_precision", "average_precision_q", "source_recall", "source_fdr",
                        "direction_accuracy", "log2fc_mae", "null_pvalue_fpr", "endpoint_fdr",
                        "endpoint_recall", "spillover_share"):
                score[key] = math.nan
        summaries.append(score)
    save_table(pd.DataFrame(summaries), out / "scores.tsv")
    save_table(pd.concat(details, ignore_index=True), out / "feature_evaluation.tsv")
    save_table(pd.DataFrame([recovery]), out / "recovery.tsv")
    write_json(out / "method_statuses.json", statuses)


def paired_bootstrap(frame, statistic, group_columns, seed=214, draws=2000):
    """Resample experiments within each condition, keeping methods paired."""
    if frame.empty:
        return {}
    strata = [c for c in group_columns if c != "method"]
    out = {}
    for values, group in frame.groupby(strata, sort=True, dropna=False):
        values = values if isinstance(values, tuple) else (values,)
        meta = dict(zip(strata, values))
        pivot = group.pivot_table(index="case", columns="method", values=statistic, aggfunc="mean")
        if len(pivot) < 2:
            continue
        picked = np.random.default_rng(seed).integers(0, len(pivot), (draws, len(pivot)))
        for method in pivot:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                boot = np.nanmean(pivot[method].to_numpy()[picked], axis=1)
            boot = boot[np.isfinite(boot)]
            if boot.size:
                key = tuple(method if c == "method" else meta[c] for c in group_columns)
                out[key] = tuple(float(v) for v in np.quantile(boot, [0.025, 0.975]))
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
    conditions = ([c for c in CONDITION_COLUMNS if c in scores.columns and scores[c].notna().all()]
                  if "generator_version" in scores.columns else [])
    if "generator_version" in scores and scores.generator_version.astype(str).str.startswith("tier1-v").any():
        if set(conditions) != set(CONDITION_COLUMNS):
            raise ValueError("Corrected Tier 1 scores require complete condition metadata; mixed or incomplete runs cannot be pooled")
    grouping = ["kind", "metric", "endpoint", "method", *conditions]
    comparison = ["kind", "metric", "endpoint", *conditions]
    intervals = paired_bootstrap(scores.loc[scores.scenario.eq("spiked") & scores.status.eq("success")], "fdr", grouping)
    # Null runs are matched on the factors that shape null data; effect factors are ignored.
    null_key = ["kind", "metric", "endpoint", "method", *[c for c in conditions if c in NULL_CONDITION_COLUMNS]]
    null_runs = {key if isinstance(key, tuple) else (key,): rows for key, rows
                 in scores.loc[scores.scenario.eq("null")].groupby(null_key, sort=False)}
    rows = []
    for values, group in scores.groupby(grouping, sort=True):
        shared = null_runs.get(tuple(dict(zip(grouping, values))[c] for c in null_key))
        group = pd.concat([group.loc[group.scenario.ne("null")]] + ([shared] if shared is not None else []))
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
                   endpoint_fdr=float(spike.endpoint_fdr.mean()) if "endpoint_fdr" in spike else math.nan,
                   n_endpoint_fdr_estimable=int(spike.endpoint_fdr.notna().sum()) if "endpoint_fdr" in spike else 0,
                   endpoint_recall=float(spike.endpoint_recall.mean()) if "endpoint_recall" in spike else math.nan,
                   spillover_share=float(spike.spillover_share.mean()) if "spillover_share" in spike else math.nan,
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
        # Pre-registered eligibility (plan.md Section 4), applied per method x
        # endpoint x condition in paper analysis; recorded here globally as well.
        # Upper FDR CI <=0.10 (twice nominal), null p-level FPR <=0.07 and
        # failure rate <5%. A null FPR below 0.03 is valid but flagged conservative.
        # Null any-FD with exact interval is reported above.
        n_fail = int(row["n_failed"] + row["n_empty"])
        row["failure_rate"] = n_fail / len(group) if len(group) else math.nan
        if "null_pvalue_fpr" in null.columns and len(null):
            fpr_vals = pd.to_numeric(null["null_pvalue_fpr"], errors="coerce").dropna().to_numpy()
            row["null_fpr"] = float(fpr_vals.mean()) if fpr_vals.size else math.nan
        else:
            row["null_fpr"] = math.nan
            fpr_vals = np.array([])
        if fpr_vals.size >= 2:
            rng = np.random.default_rng(215)
            boot = rng.choice(fpr_vals, (2000, len(fpr_vals)), replace=True).mean(axis=1)
            row["null_fpr_ci_low"], row["null_fpr_ci_high"] = (float(np.quantile(boot, 0.025)),
                                                               float(np.quantile(boot, 0.975)))
        else:
            row["null_fpr_ci_low"] = row["null_fpr_ci_high"] = math.nan
        fdr_ci_high = row.get("fdr_ci_high", math.nan)
        try:
            fdr_ok = bool(np.isfinite(fdr_ci_high) and fdr_ci_high <= 2 * cfg["analysis"]["alpha"])
        except Exception:
            fdr_ok = False
        null_ok = bool(np.isfinite(row["null_fpr"]) and row["null_fpr"] <= NULL_FPR_MAX)
        row["null_conservative"] = bool(np.isfinite(row["null_fpr"]) and row["null_fpr"] < NULL_FPR_CONSERVATIVE)
        row["prereg_eligible"] = bool(fdr_ok and null_ok and row["failure_rate"] < 0.05
                                      and row["source_recall"] > 0
                                      and not row["unassigned_significant"]
                                      and len(spike) > 0 and n_null_runs > 0)
        row["prereg_note"] = ("; ".join([
            f"fdr_ci_high={fdr_ci_high:.3f}" if np.isfinite(fdr_ci_high) else "fdr_ci_high=NA",
            f"null_fpr={row['null_fpr']:.3f}" if np.isfinite(row["null_fpr"]) else "null_fpr=NA",
            f"failure_rate={row['failure_rate']:.3f}" if np.isfinite(row["failure_rate"]) else "failure_rate=NA",
        ]))
        rows.append(row)
    eligibility_key = "prereg_eligible" if conditions else "eligible"
    ap_key = "average_precision" if conditions else "average_precision_q"
    ranking = pd.DataFrame(rows).sort_values([*comparison, eligibility_key,
                                             "source_recall", ap_key, "fdr"],
                                            ascending=[*[True] * len(comparison), False, False, False, True],
                                            na_position="last")
    ranking["rank"] = ranking.groupby(comparison).cumcount() + 1
    save_table(ranking, out / "ranking.tsv")
    if settings_rows:
        settings = pd.DataFrame(settings_rows)
        settings.to_csv(out / "method_settings.tsv", sep="\t", index=False, na_rep="NA")
    recommendations = []
    for values, group in ranking.groupby(comparison):
        enough = bool(conditions) and group.n_null.min() >= 100 and group.n_spiked.min() >= 50
        eligible = group.loc[group.prereg_eligible if conditions else group.eligible]
        if conditions and not enough:
            eligible = eligible.iloc[:0]
        candidates = []
        if len(eligible):
            best = eligible.iloc[0]
            tied = eligible.loc[np.isclose(eligible.source_recall, best.source_recall)
                                & np.isclose(eligible[ap_key], best[ap_key], equal_nan=True)
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
        enough = enough or (not conditions and len(cfg["experiment"]["seeds"]) >= 5 and cfg["genomes"]["manifest"] is not None)
        recommendations.append(dict(zip(comparison, [v.item() if isinstance(v, np.generic) else v for v in values]),
                                    best_method=best_method, candidate_methods=candidates,
                                    candidate_families=families,
                                    evidence="repeated_simulation" if enough else "smoke_only",
                                    conclusion=conclusion))
    write_json(out / "recommendations.json", recommendations)
