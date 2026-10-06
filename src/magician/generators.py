"""Versioned Tier 1 simulations with shared null/spiked noise mechanisms.

All three generators take their abundance spread, between-sample variation and
sparsity from one empirical donor table (resources/real_profiles): species
profiles of 999 healthy adults from a single study and body site.

G1 is parametric: lognormal baseline abundances and independent lognormal
biological variation, with both standard deviations set to the donor table's.
G2 resamples whole donor profiles without replacement, zeros included, so real
between-feature dependence is retained. G3 is a synthetic hurdle model: each
feature has the prevalence and nonzero lognormal moments of a donor feature,
and presence and abundance follow prescribed block correlations. It is not a
fitted SparseDOSSA model. Sequencing counts are multinomial in every generator.

The implanted set is the design truth (`is_da`). One-directional implants also
shift the relative abundance of every other feature through closure; that
shift is recorded (`is_relative_shifted`, `true_log2fc`) but is not design truth.
Null data do not depend on effect factors, so null cases are generated once per
null condition and shared by every spiked cell with that condition.
Matrix features have equal lengths so read fractions and length-corrected
abundance use the same scale in Tier 1. Batch confounding is excluded until
covariate-aware adapters are implemented.
"""
from functools import lru_cache
from pathlib import Path
from statistics import NormalDist
import gzip
import hashlib
import numpy as np
import pandas as pd
from .io import save_table, write_json, read_json

GENERATOR_VERSION = "tier1-v3"
GENERATORS = ("g1", "g2", "g3")
GENERATOR_CHOICES = list(GENERATORS)
PROFILE_DIR = Path(__file__).resolve().parents[2] / "resources/real_profiles"
TIER1_AXES = {
    "samples_per_group": [5, 10, 20, 40],
    "n_features": [100, 300, 1000],
    "depth": [500_000, 2_000_000, 10_000_000],
    "fraction_da": [0.05, 0.10, 0.20, 0.40],
    "log2fc": [0.5, 1.0, 2.0, 4.0],
    "direction": ["balanced", "all-up", "all-down"],
    "total_load": ["none", "x2"],
    "zero_structure": ["sampling", "dropout", "prevalence-dependent"],
    "confounder": ["none"],
}
# Factors that change null data. Effect factors only act on spiked cases.
NULL_FACTORS = ("generator", "samples_per_group", "n_features", "depth", "zero_structure", "confounder")
NO_EFFECT = dict(fraction_da=0.0, log2fc=0.0, direction="none", total_load="none")


def _close(values):
    values = np.maximum(np.asarray(values, dtype=float), 0)
    if not np.isfinite(values).all() or values.sum() <= 0:
        raise ValueError("A community must have finite, positive total abundance")
    return values / values.sum()


def _effect_multipliers(n_changed, log2fc, direction, rng):
    fc = 2.0 ** float(log2fc)
    if direction == "balanced":
        mult = np.full(n_changed, 1.0 / fc)
        mult[:n_changed // 2] = fc
        rng.shuffle(mult)
        return mult
    if direction in ("all-up", "all-down"):
        return np.full(n_changed, fc if direction == "all-up" else 1.0 / fc)
    raise ValueError(f"Unknown direction: {direction}")


def _apply_zeros(counts, expected_props, zero_structure, rng, depth):
    """Detection loss after counting; sampling zeros come from counts alone."""
    if zero_structure in ("none", "sampling"):
        return counts.copy()
    if zero_structure == "dropout":
        order = np.argsort(expected_props, axis=0)
        third = max(1, counts.shape[0] // 3)
        mask = np.zeros_like(counts, dtype=bool)
        for col in range(counts.shape[1]):
            low, mid = order[:third, col], order[third:2 * third, col]
            mask[low, col] = rng.random(len(low)) < 0.25
            mask[mid, col] = rng.random(len(mid)) < 0.05
        return np.where(mask, 0, counts)
    if zero_structure == "prevalence-dependent":
        detect = 1.0 - np.exp(-np.maximum(expected_props * depth, 0) / 5.0)
        return np.where(rng.random(counts.shape) < detect, counts, 0)
    raise ValueError(f"Unknown zero_structure: {zero_structure}")


@lru_cache(maxsize=1)
def _real_profile_data():
    """Donor proportions (samples x features, most prevalent feature first) and provenance."""
    provenance = read_json(PROFILE_DIR / "provenance.json")
    raw = gzip.decompress((PROFILE_DIR / provenance["table"]).read_bytes())
    if hashlib.sha256(raw).hexdigest() != provenance["table_sha256"]:
        raise ValueError("The empirical donor table differs from its provenance record")
    frame = pd.read_csv(PROFILE_DIR / provenance["table"], sep="\t", index_col=0)
    matrix = frame.to_numpy(dtype=float).T
    if not np.isfinite(matrix).all() or (matrix < 0).any() or (matrix.sum(axis=1) <= 0).any():
        raise ValueError("The empirical donor table must hold finite, nonnegative, nonempty profiles")
    prevalence, mean = (matrix > 0).mean(axis=0), matrix.mean(axis=0)
    order = np.lexsort((np.arange(matrix.shape[1]), -mean, -prevalence))
    matrix = matrix[:, order]
    matrix.setflags(write=False)
    keep = ("study", "resource", "profiler", "profile_sha256", "table", "table_sha256", "n_samples", "n_features")
    return matrix, {k: provenance[k] for k in keep}


@lru_cache(maxsize=8)
def _empirical_pool(n_features):
    """Whole donor profiles over the most prevalent features, zeros as measured."""
    real, provenance = _real_profile_data()
    if n_features > real.shape[1]:
        raise ValueError(f"G2 cannot supply {n_features} features: the donor table has {real.shape[1]} "
                         "measured species and synthetic expansion is disabled")
    pool = real[:, :n_features]
    pool = pool[pool.sum(axis=1) > 0]
    pool = pool / pool.sum(axis=1, keepdims=True)
    pool.setflags(write=False)
    return pool, dict(provenance, n_donors=int(len(pool)), features="most prevalent species")


@lru_cache(maxsize=8)
def _calibration(n_features):
    """Donor moments: prevalence and nonzero log-abundance mean and SD per feature."""
    real, provenance = _real_profile_data()
    pool, _ = _empirical_pool(min(n_features, real.shape[1]))
    logged = np.where(pool > 0, np.log(np.where(pool > 0, pool, 1.0)), np.nan)
    present = (pool > 0).sum(axis=0)
    log_sd = np.where(present >= 3, np.nanstd(logged, axis=0), np.nan)
    typical_sd = float(np.nanmedian(log_sd))
    moments = dict(prevalence=(pool > 0).mean(axis=0), log_mean=np.nanmean(logged, axis=0),
                   log_sd=np.where(np.isfinite(log_sd), log_sd, typical_sd),
                   abundance_sd=float(np.log(pool.mean(axis=0)).std()), biological_sd=typical_sd)
    return moments, dict(provenance, calibrated_on=int(pool.shape[1]))


def _implant(baseline, spec, scenario, rng, stratified=False):
    """Preserve balanced-subset mass; close the whole one-directional community.

    One-directional implants change the relative abundances of unselected taxa
    through closure. Truth records those shifts separately from the implanted set.
    Total-load change is independent of these relative effects.
    """
    n_features = len(baseline)
    n_changed = min(n_features, max(2, round(n_features * spec["fraction_da"]))) if spec["fraction_da"] > 0 else 0
    factors = np.ones(n_features)
    implanted = np.zeros(n_features, dtype=bool)
    if not n_changed:
        return _close(baseline), factors, implanted
    if stratified:
        strata = [s for s in np.array_split(np.argsort(baseline), min(10, n_features)) if len(s)]
        selected = []
        while len(selected) < n_changed:
            for stratum in strata:
                candidates = np.setdiff1d(stratum, selected)
                if len(candidates) and len(selected) < n_changed:
                    selected.append(int(rng.choice(candidates)))
        chosen = np.asarray(selected)
    else:
        chosen = rng.choice(n_features, n_changed, replace=False)
    multipliers = _effect_multipliers(n_changed, spec["log2fc"], spec["direction"], rng)
    if scenario == "spiked":
        factors[chosen] = multipliers
        implanted[chosen] = True
        if spec["direction"] == "balanced":
            factors[chosen] *= baseline[chosen].sum() / (baseline[chosen] * factors[chosen]).sum()
    return _close(baseline * factors), factors, implanted


def _write_case(out, ids, lengths, samples, counts, changed, loads):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    names = [s["sample_id"] for s in samples]
    pd.DataFrame(counts, index=pd.Index(ids, name="feature_id"), columns=names).to_csv(out / "counts.tsv", sep="\t")
    clean = [{k: v for k, v in sample.items() if k not in ("composition", "target")} for sample in samples]
    save_table(pd.DataFrame(clean), out / "samples.tsv")
    save_table(pd.DataFrame({"feature_id": ids, "length_bp": lengths.astype(int)}), out / "features.tsv")
    save_table(changed, out / "truth.tsv")
    n = len(ids)
    props = np.concatenate([s["composition"] for s in samples])
    targets = np.concatenate([s["target"] for s in samples])
    pairs = np.repeat([int(s["read_pairs"]) for s in samples], n)
    load = np.repeat(np.asarray(loads, dtype=float), n)
    save_table(pd.DataFrame(dict(sample_id=np.repeat(names, n), feature_id=np.tile(ids, len(samples)),
                                 genomic_proportion=props, target_proportion=targets,
                                 read_proportion=props, library_read_pairs=pairs,
                                 length_bp=np.tile(lengths, len(samples)), expected_count=props * pairs,
                                 absolute_copies=props * load, target_absolute_copies=targets * load)),
               out / "expectations.tsv")
    save_table(pd.DataFrame({"feature_id": ids, "source_id": ids, "match_status": "matched"}), out / "matching.tsv")


def generate_case(generator, spec, scenario, seed, output):
    if generator not in GENERATORS or scenario not in ("null", "spiked"):
        raise ValueError(f"Unknown generator/scenario: {generator}/{scenario}")
    if spec.get("confounder", "none") != "none":
        raise ValueError("Batch confounding is excluded until simulated batch effects and adjusted adapters are supported")
    if not 0 <= float(spec["fraction_da"]) <= 1:
        raise ValueError("fraction_da must be within [0, 1]")
    if spec["total_load"] not in ("none", "x2"):
        raise ValueError("Unknown total_load")
    design_rng, noise_rng, count_rng = [np.random.default_rng(s) for s in np.random.SeedSequence(seed).spawn(3)]
    n_features, per_group, depth = int(spec["n_features"]), int(spec["samples_per_group"]), int(spec["depth"])
    if n_features < 2 or per_group < 3 or depth < 1:
        raise ValueError("Cases require >=2 features, >=3 samples/group and positive depth")
    if generator == "g2":
        pool, provenance = _empirical_pool(n_features)
        if len(pool) < 2 * per_group:
            raise ValueError("G2 needs one distinct donor profile per sample")
        baseline = _close(pool.mean(axis=0))
        donors = noise_rng.choice(len(pool), 2 * per_group, replace=False)
    else:
        moments, provenance = _calibration(n_features)
        if generator == "g1":
            abundance_sd = float(spec.get("abundance_sd", moments["abundance_sd"]))
            sd = float(spec.get("biological_sd", moments["biological_sd"]))
            baseline = _close(np.exp(design_rng.normal(0, abundance_sd, n_features)))
            provenance = dict(provenance, marginals="lognormal", abundance_sd=abundance_sd, biological_sd=sd)
        else:
            available = len(moments["prevalence"])
            picked = (np.arange(n_features) if n_features <= available
                      else np.r_[np.arange(available), design_rng.integers(0, available, n_features - available)])
            prevalence = np.clip(moments["prevalence"][picked], 1e-6, 1 - 1e-6)
            log_mean, log_sd = moments["log_mean"][picked], moments["log_sd"][picked]
            absent_below = np.asarray([NormalDist().inv_cdf(1 - p) for p in prevalence])
            baseline = _close(prevalence * np.exp(log_mean + log_sd ** 2 / 2))
            rho = float(spec.get("copula_rho", 0.3))
            if not 0 <= rho < 1:
                raise ValueError("copula_rho must be within [0, 1)")
            blocks = design_rng.permutation(np.repeat(np.arange(10), (n_features + 9) // 10)[:n_features])
            provenance = dict(provenance, marginals="synthetic hurdle lognormal with donor moments",
                              within_block_rho=rho, resampled_marginals=bool(n_features > available))
    treated, factors, implanted = _implant(baseline, spec, scenario, design_rng, generator == "g2")
    delta = np.log2(treated / baseline)
    ids = [f"feature_{i + 1:04d}" for i in range(n_features)]
    lengths = np.full(n_features, 20000.0)
    changed = pd.DataFrame(dict(feature_id=ids, is_da=implanted, is_implanted=implanted,
                                is_relative_shifted=np.abs(delta) > 1e-10, true_log2fc=delta,
                                baseline_proportion=baseline, treatment_proportion=treated,
                                length_bp=lengths.astype(int)))
    samples, compositions, loads = [], [], []
    counts = np.zeros((n_features, 2 * per_group), dtype=int)
    for col in range(2 * per_group):
        treatment = col >= per_group
        total = max(1, round(depth * noise_rng.lognormal(-float(spec.get("library_cv", 0.15)) ** 2 / 2,
                                                        float(spec.get("library_cv", 0.15)))))
        if generator == "g1":
            donor = _close(baseline * noise_rng.lognormal(-sd ** 2 / 2, sd, n_features))
        elif generator == "g2":
            donor = pool[donors[col]]
        else:
            while True:
                # Presence and abundance each share one latent factor per block.
                latent = [np.sqrt(rho) * noise_rng.normal(size=10)[blocks]
                          + np.sqrt(1 - rho) * noise_rng.normal(size=n_features) for _ in range(2)]
                donor = np.where(latent[0] > absent_below, np.exp(log_mean + log_sd * latent[1]), 0.0)
                if donor.sum() > 0:
                    break
            donor = _close(donor)
        props = _close(donor * (factors if treatment else 1))
        dispersion = float(spec.get("dispersion", 1.0))
        if dispersion < 1:
            raise ValueError("dispersion must be >=1")
        count_props = props
        if dispersion > 1:
            # Dirichlet-multinomial counting variation keeps the library total fixed.
            count_props = count_rng.dirichlet(np.maximum(props * total / (dispersion - 1), 1e-9))
        counts[:, col] = count_rng.multinomial(total, count_props)
        compositions.append(props)
        load = 2.0 if scenario == "spiked" and treatment and spec["total_load"] == "x2" else 1.0
        loads.append(load)
        samples.append(dict(sample_id=f"{'T' if treatment else 'C'}{col % per_group + 1:02d}",
                            group="Treatment" if treatment else "Control", read_pairs=total,
                            composition=props, target=treated if treatment else baseline))
    counts = _apply_zeros(counts, np.column_stack(compositions), spec["zero_structure"], count_rng, depth)
    _write_case(output, ids, lengths, samples, counts, changed, loads)
    rec = dict(generator=generator, generator_version=GENERATOR_VERSION, scenario=scenario, seed=seed,
               n_features=n_features, samples_per_group=per_group, n_changed=int(implanted.sum()),
               n_implanted=int(implanted.sum()), n_relative_shifted=int(changed.is_relative_shifted.sum()),
               zero_fraction=float((counts == 0).mean()), provenance=provenance)
    write_json(Path(output) / "generator.json", dict(rec, cell=dict(spec)))
    return rec


def g1_case(spec, scenario, seed, output):
    return generate_case("g1", spec, scenario, seed, output)


def g2_case(spec, scenario, seed, output):
    return generate_case("g2", spec, scenario, seed, output)


def g3_case(spec, scenario, seed, output):
    return generate_case("g3", spec, scenario, seed, output)


GENERATOR_FN = {"g1": g1_case, "g2": g2_case, "g3": g3_case}
# A fourfold change: with donor-level between-sample variation a twofold change at 20
# samples per group is detected by almost no method, which would leave power unranked.
BASE_CELL = dict(generator="g1", samples_per_group=20, n_features=300, depth=2_000_000,
                 fraction_da=0.10, log2fc=2.0, direction="balanced", total_load="none",
                 zero_structure="sampling", confounder="none")


def null_cell(cell):
    """The null condition a spiked cell is calibrated against."""
    return dict(cell, **NO_EFFECT)


def core_cells():
    """Spiked core cells: generator x depth x fraction x direction."""
    return [dict(BASE_CELL, generator=generator, depth=depth, fraction_da=fraction, direction=direction)
            for generator in GENERATORS for depth in (500_000, 2_000_000)
            for fraction in (0.10, 0.20) for direction in ("balanced", "all-up")]


def fractional_grid():
    """One factor at a time around the base cell; cells already in the core are not repeated."""
    core = core_cells()
    cells = []
    for key, levels in TIER1_AXES.items():
        for value in levels:
            cell = dict(BASE_CELL, **{key: value})
            if cell not in core and cell not in cells:
                cells.append(cell)
    return cells


def design_plan(cells, null_replicates, spiked_replicates, prefix, covered=()):
    """Spiked replicates per cell and null replicates once per distinct null condition."""
    nulls = []
    for cell in cells:
        null = null_cell(cell)
        if null not in nulls and null not in covered:
            nulls.append(null)
    plan = [dict(design=f"{prefix}_n{i:03d}", cell=c, scenario="null", replicates=null_replicates)
            for i, c in enumerate(nulls)]
    return plan + [dict(design=f"{prefix}_c{i:03d}", cell=c, scenario="spiked", replicates=spiked_replicates)
                   for i, c in enumerate(cells)]


def tier1_plan(null_replicates=200, spiked_replicates=100, fractional_replicates=20):
    core = design_plan(core_cells(), null_replicates, spiked_replicates, "core")
    covered = [entry["cell"] for entry in core if entry["scenario"] == "null"]
    return core + design_plan(fractional_grid(), fractional_replicates, fractional_replicates,
                              "fractional", covered)


def build_bundle(plan, output, first_seed=1000, progress=False):
    """Build a fresh immutable bundle; never overwrite earlier experiment data."""
    out = Path(output)
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"Bundle {out} is not empty; use a new versioned directory")
    if not plan or any(not 1 <= entry["replicates"] < 100000 for entry in plan):
        raise ValueError("Replicate counts must be within [1, 99999]")
    if len({entry["design"] for entry in plan}) != len(plan):
        raise ValueError("Design identifiers must be unique")
    out.mkdir(parents=True, exist_ok=True)
    listing, summary = [], []
    for pos, entry in enumerate(plan):
        cell, scenario = entry["cell"], entry["scenario"]
        if progress:
            print(f"building {entry['design']} ({scenario} x {entry['replicates']}): {cell}", flush=True)
        for rep in range(entry["replicates"]):
            seed = first_seed + 100000 * pos + rep
            case = f"{scenario}_s{seed}"
            rec = generate_case(cell["generator"], cell, scenario, seed, out / case)
            rec.update(case=case, design=entry["design"])
            listing.append(dict(case=case, scenario=scenario, seed=seed, design=entry["design"],
                                generator_version=GENERATOR_VERSION, **cell))
            summary.append({k: v for k, v in rec.items() if k != "provenance"})
    save_table(pd.DataFrame(listing), out / "cases.tsv")
    save_table(pd.DataFrame(summary), out / "designs.tsv")
    write_json(out / "designs.json", dict(plan=plan, first_seed=first_seed,
                                         generator_version=GENERATOR_VERSION, n_cases=len(listing),
                                         null_factors=list(NULL_FACTORS),
                                         donor_table=_real_profile_data()[1],
                                         excluded_conditions=["batch confounding"],
                                         code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    return summary


def _count_summary(counts):
    counts = np.asarray(counts, dtype=float)
    share = counts / counts.sum(axis=0, keepdims=True)
    prevalence = (counts > 0).mean(axis=1)
    seen = prevalence > 0
    logged = np.where(share > 0, np.log(np.where(share > 0, share, 1.0)), np.nan)
    with np.errstate(invalid="ignore"):
        spread = np.nanstd(np.where((counts > 0).sum(axis=1, keepdims=True) >= 3, logged, np.nan), axis=1)
    quantiles = [float(v) for v in np.quantile(prevalence, [0.1, 0.25, 0.5, 0.75, 0.9])]
    return dict(zero_fraction=float((counts == 0).mean()), prevalence_quantiles=quantiles,
                features_below_10pct_prevalence=float((prevalence < 0.1).mean()),
                median_nonzero_count=float(np.median(counts[counts > 0])),
                log_mean_abundance_sd=float(np.log(share.mean(axis=1)[seen]).std()),
                median_within_feature_log_sd=float(np.nanmedian(spread)) if np.isfinite(spread).any() else None)


def realism_check(bundle, output, cases_per_condition=20, seed=7):
    """Compare each null condition with donor profiles counted at the same depth.

    The reference is the donor table itself, multinomially counted to the same
    library sizes. G2 should match it closely; G1 and G3 show how far a
    parametric model departs. These summaries are descriptive, not a test.
    """
    bundle, out = Path(bundle), Path(output)
    out.mkdir(parents=True, exist_ok=True)
    listing = pd.read_csv(bundle / "cases.tsv", sep="\t", keep_default_na=False)
    listing = listing.loc[listing.scenario.eq("null")]
    real, provenance = _real_profile_data()
    rng = np.random.default_rng(seed)
    rows = []
    for design, group in listing.groupby("design", sort=True):
        first = group.iloc[0]
        simulated = [pd.read_csv(bundle / case / "counts.tsv", sep="\t", index_col=0).to_numpy()
                     for case in group.case.head(cases_per_condition)]
        row = dict(design=design, **{k: (first[k].item() if hasattr(first[k], "item") else first[k])
                                    for k in NULL_FACTORS}, n_cases=len(simulated))
        summaries = [_count_summary(c) for c in simulated]
        row["simulated"] = {k: (np.mean([s[k] for s in summaries], axis=0).tolist()
                                if all(s[k] is not None for s in summaries) else None) for k in summaries[0]}
        n = int(first.n_features)
        if n <= real.shape[1]:
            pool, _ = _empirical_pool(n)
            reference = []
            for counts in simulated:
                picked = rng.choice(len(pool), counts.shape[1], replace=False)
                reference.append(_count_summary(np.column_stack(
                    [rng.multinomial(int(total), pool[i]) for i, total in zip(picked, counts.sum(axis=0))])))
            row["donor_reference"] = {k: np.mean([s[k] for s in reference], axis=0).tolist() for k in reference[0]}
        else:
            row["donor_reference"] = None
        rows.append(row)
    record = dict(generator_version=GENERATOR_VERSION, donor_table=provenance, conditions=rows,
                  note="donor_reference counts donor profiles at the simulated library sizes; "
                       "zero-structure conditions add detection loss that the reference does not have")
    write_json(out / "realism.json", record)
    try:
        import os
        os.environ.setdefault("MPLCONFIGDIR", str(out / ".mpl"))
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return record
    shown = [r for r in rows if r["donor_reference"]]
    keys = ["zero_fraction", "features_below_10pct_prevalence", "log_mean_abundance_sd",
            "median_within_feature_log_sd"]
    fig, axes = plt.subplots(1, len(keys), figsize=(4.2 * len(keys), 3.6), constrained_layout=True)
    for ax, key in zip(axes, keys):
        ax.scatter([r["donor_reference"][key] for r in shown], [r["simulated"][key] for r in shown],
                   c=[GENERATORS.index(r["generator"]) for r in shown], cmap="viridis")
        low, high = ax.get_xlim()
        ax.plot([low, high], [low, high], color="grey", linewidth=1)
        ax.set_xlabel("donor profiles"), ax.set_ylabel("simulated"), ax.set_title(key, fontsize=9)
    fig.savefig(out / "realism.png", dpi=140)
    plt.close(fig)
    return record
