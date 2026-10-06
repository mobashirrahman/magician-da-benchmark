"""Tier 1 calibration at the p-value level + cluster bootstrap.

Null calibration is judged at the p-value level (QQ plots, false-positive rate
at alpha 0.05, KS on uniformity), which yields thousands of null tests per
cell, and by the any-false-discovery rate only on the core cells. Features
within an experiment are not independent, so intervals use the cluster
(per-experiment) bootstrap: resample whole experiments, keeping every method
paired.
"""
import math
import numpy as np
import pandas as pd


def ks_uniformity(pvalues):
    """Kolmogorov-Smirnov statistic against Uniform(0,1); NaN when empty."""
    p = np.asarray([v for v in np.ravel(pvalues) if np.isfinite(v)], dtype=float)
    if p.size == 0:
        return math.nan
    p = np.clip(p, 1e-12, 1 - 1e-12)
    ordered = np.sort(p)
    expected = (np.arange(1, len(p) + 1) - 0.5) / len(p)
    return float(np.max(np.abs(ordered - expected)))


def fpr_at_alpha(pvalues, alpha=0.05):
    """False-positive rate at alpha on null p-values."""
    p = np.asarray([v for v in np.ravel(pvalues) if np.isfinite(v)], dtype=float)
    if p.size == 0:
        return math.nan
    return float((p <= alpha).mean())


def qq_data(pvalues, points=200):
    """Observed vs expected -log10 p for a QQ plot (downsampled)."""
    p = np.asarray([v for v in np.ravel(pvalues) if np.isfinite(v)], dtype=float)
    if p.size == 0:
        return {"expected": [], "observed": []}
    p = np.clip(p, 1e-300, 1.0)
    ordered = np.sort(p)
    idx = np.unique(np.linspace(0, len(ordered) - 1, min(points, len(ordered))).astype(int))
    obs = -np.log10(ordered[idx])
    exp = -np.log10((idx + 0.5) / len(ordered))
    return {"expected": [float(v) for v in exp], "observed": [float(v) for v in obs]}


def cluster_bootstrap_ci(frame, statistic, group_columns, seed=214, draws=2000, level=0.95):
    """Cluster bootstrap: resample whole cases, keep method comparison paired.

    Returns {(group values): (low, high)} with quantile intervals.
    """
    if frame.empty or "case" not in frame.columns:
        return {}
    cases = sorted(frame["case"].unique())
    if len(cases) < 2:
        return {}
    rng = np.random.default_rng(seed)
    picked = rng.integers(0, len(cases), (draws, len(cases)))
    alpha = (1 - level) / 2
    out = {}
    for values, group in frame.groupby(list(group_columns), sort=True):
        means = []
        for case in cases:
            vals = group.loc[group["case"].eq(case), statistic].dropna().to_numpy()
            means.append(vals.mean() if vals.size else np.nan)
        means = np.asarray(means, dtype=float)
        boot = np.nanmean(means[picked], axis=1)
        boot = boot[np.isfinite(boot)]
        if not boot.size:
            continue
        key = tuple(values) if isinstance(values, tuple) else (values,)
        out[key] = (float(np.quantile(boot, alpha)), float(np.quantile(boot, 1 - alpha)))
    return out


def calibration_table(feature_evaluation, alpha=0.05):
    """Per method x case null calibration from feature-level p-values.

    Expects columns: case, method, method_family, pvalue, truth_is_da (or is_da).
    Null features are those with truth_is_da == False on null cases.
    """
    frame = feature_evaluation.copy()
    label_col = "truth_is_da" if "truth_is_da" in frame.columns else "is_da"
    rows = []
    for (case, method, family), group in frame.groupby(["case", "method", "method_family"], sort=True):
        null_p = group.loc[~group[label_col].astype(str).str.lower().eq("true"), "pvalue"]
        null_p = pd.to_numeric(null_p, errors="coerce").dropna().to_numpy()
        rows.append(dict(case=case, method=method, method_family=family,
                         n_null_p=int(null_p.size),
                         null_fpr=float((null_p <= alpha).mean()) if null_p.size else math.nan,
                         null_ks=ks_uniformity(null_p),
                         null_qq=qq_data(null_p)))
    return pd.DataFrame(rows)
