"""Variance decomposition, mixed-effects summaries and rank uncertainty.

Headline quantification of "what matters": which axis (method, depth, DA
fraction, direction, generator, tier) explains the most variance in FDP and
sensitivity. Rankings carry rank-uncertainty intervals from the bootstrap,
never bare ranks.

These are descriptive sequential variance shares and paired bootstrap summaries.
They do not fit a mixed-effects model. Contrasts against the observed best method
are exploratory; missing bootstrap values retain their draw indices.
"""
import itertools
import math
import numpy as np
import pandas as pd


def variance_decomposition(scores, factors, statistic="fdr"):
    """ANOVA-style variance share per factor (type-I sequential, in given order).

    Returns {factor: share} summing to <=1 (residual remainder).
    """
    frame = scores.dropna(subset=[statistic]).copy()
    if frame.empty:
        return {}
    y = frame[statistic].to_numpy(dtype=float)
    total = float(((y - y.mean()) ** 2).sum())
    if total <= 0:
        return {f: 0.0 for f in factors}
    out, explained = {}, 0.0
    # Sequential: residualise on previous factors via group means.
    resid = y.copy()
    for factor in factors:
        if factor not in frame.columns:
            out[factor] = 0.0
            continue
        means = frame.assign(_r=resid).groupby(factor)["_r"].transform("mean").to_numpy()
        ss = float(((means - resid.mean()) ** 2).sum())
        # Orthogonalise for the next factor.
        resid = resid - (means - resid.mean())
        share = max(0.0, ss / total)
        out[factor] = share
        explained += share
    out["residual"] = max(0.0, 1.0 - explained)
    return out


def holm(pvalues):
    """Holm-adjusted p-values (no SciPy dependency)."""
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p, kind="stable")
    m = len(p)
    adj = np.empty(m)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (m - rank) * p[idx])
        adj[idx] = min(1.0, running)
    return adj


def _bootstrap_groups(scores, group_columns, statistic, seed, draws):
    frame = scores.dropna(subset=[statistic, "case"])
    if frame.empty or frame.case.nunique() < 2:
        return {}, {}
    pivot = frame.pivot_table(index="case", columns=list(group_columns), values=statistic, aggfunc="mean")
    sampled = np.random.default_rng(seed).integers(0, len(pivot), (draws, len(pivot)))
    boot, observed = {}, {}
    for column in pivot:
        key = column if isinstance(column, tuple) else (column,)
        values = pivot[column].to_numpy()
        selected = values[sampled]
        count = np.isfinite(selected).sum(axis=1)
        means = np.divide(np.nansum(selected, axis=1), count,
                          out=np.full(draws, np.nan), where=count > 0)
        boot[key] = means
        observed[key] = float(np.nanmean(values))
    return boot, observed


def pairwise_contrasts(scores, group="method", statistic="fdr", seed=214, draws=2000):
    """Exploratory paired-bootstrap contrasts against the observed best group."""
    boot, observed = _bootstrap_groups(scores, (group,), statistic, seed, draws)
    if len(boot) < 2:
        return pd.DataFrame()
    best = min(observed, key=observed.get)
    rows = []
    for key in sorted(boot):
        diff = boot[key] - boot[best]
        diff = diff[np.isfinite(diff)]
        if key == best:
            row = dict(group=key[0], vs_best=best[0], mean_diff=0.0,
                       ci_low=math.nan, ci_high=math.nan, p_raw=1.0)
        elif not diff.size:
            row = dict(group=key[0], vs_best=best[0], mean_diff=math.nan,
                       ci_low=math.nan, ci_high=math.nan, p_raw=1.0)
        else:
            p = 2 * min((diff <= 0).mean(), (diff >= 0).mean())
            low, high = np.quantile(diff, [0.025, 0.975])
            row = dict(group=key[0], vs_best=best[0], mean_diff=float(diff.mean()),
                       ci_low=float(low), ci_high=float(high),
                       p_raw=max(float(p), 1 / len(diff)))
        rows.append(row)
    for row, adjusted in zip(rows, holm([r["p_raw"] for r in rows])):
        row["p_holm"] = float(adjusted)
        row["interpretation"] = "exploratory; comparator selected from observed data"
    return pd.DataFrame(rows)


def rank_uncertainty(scores, group_columns=("method",), statistic="source_recall",
                     higher_is_better=True, seed=214, draws=2000):
    """Paired experiment bootstrap ranks; equal scores receive equal mean ranks."""
    boot, _ = _bootstrap_groups(scores, group_columns, statistic, seed, draws)
    if not boot:
        return {}
    keys = list(boot)
    values = pd.DataFrame({i: boot[key] for i, key in enumerate(keys)})
    ranks = values.rank(axis=1, method="average", ascending=not higher_is_better)
    out = {}
    for i, key in enumerate(keys):
        finite = ranks[i].dropna().to_numpy()
        if finite.size:
            out[key] = tuple(float(v) for v in np.quantile(finite, [0.5, 0.025, 0.975]))
    return out
