"""Tier 3a real nulls, Tier 3b implants, mock anchor and Tier 4 replication.

Tier 3a: healthy-cohort metagenomes (>=100 samples, one study/site/similar
depth) subsampled to a common depth, assembled+binned once, then 200 random
label splits (balanced, and confounded with batch where metadata allows). Any
discovery is false: measures type I error with real structure at no assumption
cost.

Tier 3b: read-level implants in group-B samples (up/down-sample reads aligned
to chosen MAGs by a known factor, conserving library size). Truth exact by
construction; assembly/bins/noise/batch real. 4 effects x 3 fractions x 10 seeds.

Tier 3c: mock communities (ZymoBIOMICS-type log ratios, public CAMI data) as a
sanity anchor that quantification recovers known ratios; small, cheap, never
used for ranking.

Tier 4: IBD (IBDMDB/HMP2 + independent IBD) and CRC (>=3 public cohorts).
Without ground truth report only cross-cohort replication (rank agreement,
replicated-discovery rate per method) and between-method agreement.
Replication is reliability, not accuracy.

Nothing in Tier 3-4 tunes methods or thresholds: thresholds and analysis code
are frozen and committed before Tier 2 (see tools/preregister.py).
"""
import numpy as np
import pandas as pd


def random_splits(sample_ids, n_splits=200, seed=0, batch=None):
    """Balanced random label splits; optionally stratified/confounded with batch.

    Returns a list of {sample_id: label} dicts with Control/Treatment labels.
    When batch is given ({sample_id: batch}), half the splits are balanced
    within batch and half are confounded with batch (all of one batch in one
    group) where metadata allows.
    """
    rng = np.random.default_rng(seed)
    ids = list(sample_ids)
    n = len(ids)
    if n < 6 or n % 2:
        raise ValueError("Need an even number of >=6 samples for balanced splits")
    half = n // 2
    splits = []
    for i in range(n_splits):
        order = rng.permutation(n)
        labels = {}
        if batch is not None and i % 2 == 1:
            # Confounded split: sort by batch, then cut (batch predicts group).
            by_batch = sorted(ids, key=lambda s: (str(batch.get(s, "")), rng.random()))
            for s in by_batch[:half]:
                labels[s] = "Control"
            for s in by_batch[half:]:
                labels[s] = "Treatment"
        else:
            for k, s in enumerate([ids[j] for j in order[:half]]):
                labels[s] = "Control"
            for k, s in enumerate([ids[j] for j in order[half:2 * half]]):
                labels[s] = "Treatment"
        splits.append(labels)
    return splits


def implant_reads(counts, mag_ids, factor, seed=0, conserve_library=True):
    """Up/down-sample reads aligned to chosen MAGs by a known factor.

    counts: DataFrame features x samples (group-B columns implanted).
    factor >1 enriches, <1 depletes. Library size conserved by rescaling the
    non-implanted features (Nearing-style implants extended to MAG space).
    Returns (implanted counts, truth labels).
    """
    rng = np.random.default_rng(seed)
    out = counts.copy().astype(float)
    # Poisson resampling around the scaled expectation keeps integer counts.
    for mag in mag_ids:
        if mag not in out.index:
            raise ValueError(f"Unknown MAG for implant: {mag}")
        scaled = np.maximum(out.loc[mag].to_numpy() * factor, 0)
        out.loc[mag] = rng.poisson(scaled).astype(float)
    if conserve_library:
        orig_totals = counts.sum(axis=0)
        new_totals = out.sum(axis=0).replace(0, np.nan)
        scale = (orig_totals / new_totals).fillna(1.0)
        non_implant = [f for f in out.index if f not in set(mag_ids)]
        out.loc[non_implant] = (out.loc[non_implant] * scale).round()
        out.loc[list(mag_ids)] = out.loc[list(mag_ids)].round()
        # Fix rounding drift on the implanted rows.
        drift = orig_totals - out.sum(axis=0)
        for col in out.columns:
            if drift[col] != 0 and len(non_implant):
                pick = non_implant[int(rng.integers(len(non_implant)))]
                out.loc[pick, col] = max(0, out.loc[pick, col] + drift[col])
    truth = pd.Series(False, index=counts.index)
    truth.loc[[m for m in mag_ids if m in truth.index]] = True
    return out.astype(int), truth


def replication_rate(discoveries_a, discoveries_b):
    """Replicated-discovery rate: |A cap B| / |A union B| (Jaccard) + rank agreement.

    discoveries_* are ranked feature lists (ordered) or sets. Returns dict with
    jaccard, overlap, spearman rank agreement on the shared ranking.
    """
    a = list(discoveries_a)
    b = list(discoveries_b)
    sa, sb = set(a), set(b)
    union = sa | sb
    jaccard = len(sa & sb) / len(union) if union else float("nan")
    # Rank agreement: Spearman on ranks of shared features.
    shared = sorted(sa & sb)
    if len(shared) < 3:
        spearman = float("nan")
    else:
        ra = pd.Series({f: i for i, f in enumerate(a)})
        rb = pd.Series({f: i for i, f in enumerate(b)})
        spearman = float(ra.loc[shared].corr(rb.loc[shared], method="spearman"))
    return dict(jaccard=jaccard, n_a=len(sa), n_b=len(sb), n_shared=len(sa & sb),
                spearman=spearman)


def cross_cohort_table(replicated):
    """Aggregate replicated-discovery rows into a per-method summary.

    replicated: DataFrame with method, cohort_pair, jaccard, spearman.
    """
    rows = []
    for (method,), group in replicated.groupby(["method"], sort=True):
        rows.append(dict(method=method, n_pairs=len(group),
                         mean_jaccard=float(group.jaccard.mean()),
                         mean_spearman=float(group.spearman.mean())))
    return pd.DataFrame(rows).sort_values("mean_jaccard", ascending=False)
