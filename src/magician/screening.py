"""Predeclared matrix screening designs and their replicate generator.

Matrix screening never needs reads, assembly or binning: a design is a latent
composition, and Dirichlet-multinomial counts are drawn directly from it. The
generated bundle is a versioned input to the matrix-only execution mode, so the same
adapter, registry and truth machinery run unchanged.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from .io import save_table, write_json

# A small balanced subset, not the full product of the design factors. Each entry is
# chosen to vary one factor at a time against a common baseline design.
DESIGNS = {
    "baseline": dict(n_features=200, samples_per_group=20, effect="2x",
                     changed_fraction=0.2, dispersion=1.0, library_cv=0.15,
                     read_pairs=20000),
    "small_features": dict(n_features=50, samples_per_group=20, effect="2x",
                           changed_fraction=0.2, dispersion=1.0, library_cv=0.15,
                           read_pairs=20000),
    "large_features": dict(n_features=500, samples_per_group=20, effect="2x",
                           changed_fraction=0.2, dispersion=1.0, library_cv=0.15,
                           read_pairs=20000),
    "few_samples": dict(n_features=200, samples_per_group=10, effect="2x",
                        changed_fraction=0.2, dispersion=1.0, library_cv=0.15,
                        read_pairs=20000),
    "many_samples": dict(n_features=200, samples_per_group=40, effect="2x",
                         changed_fraction=0.2, dispersion=1.0, library_cv=0.15,
                         read_pairs=20000),
    "strong_effect": dict(n_features=200, samples_per_group=20, effect="4x",
                          changed_fraction=0.2, dispersion=1.0, library_cv=0.15,
                          read_pairs=20000),
    "many_changes": dict(n_features=200, samples_per_group=20, effect="2x",
                         changed_fraction=0.4, dispersion=1.0, library_cv=0.15,
                         read_pairs=20000),
    "overdispersed": dict(n_features=200, samples_per_group=20, effect="2x",
                          changed_fraction=0.2, dispersion=6.0, library_cv=0.15,
                          read_pairs=20000),
    "low_count": dict(n_features=200, samples_per_group=20, effect="2x",
                      changed_fraction=0.2, dispersion=1.0, library_cv=0.15,
                      read_pairs=4000),
    "load_shift": dict(n_features=200, samples_per_group=20, effect="2x",
                       changed_fraction=0.2, dispersion=1.0, library_cv=0.60,
                       read_pairs=20000),
    "few_changes": dict(n_features=200, samples_per_group=20, effect="2x",
                        changed_fraction=0.1, dispersion=1.0, library_cv=0.15,
                        read_pairs=20000),
    "dropout": dict(n_features=200, samples_per_group=20, effect="2x",
                    changed_fraction=0.2, dispersion=6.0, library_cv=0.15,
                    read_pairs=4000, concentration=10.0),
}


def draw_counts(rng, expected, dispersion):
    """Poisson counts, or gamma-Poisson when overdispersion is requested."""
    expected = np.asarray(expected, dtype=float)
    if dispersion <= 1.0:
        return rng.poisson(np.maximum(expected, 0.0)).astype(int)
    alpha = max(float(dispersion), 1e-9)
    # A structurally absent feature has no count to overdisperse.
    absent = expected <= 0
    safe = np.where(absent, 1.0, expected)
    size = alpha / safe
    probability = np.clip(size / (size + safe), 1e-12, 1 - 1e-12)
    drawn = rng.negative_binomial(size, probability).astype(int)
    return np.where(absent, 0, drawn)


def sample_composition(rng, n_features, concentration=100.0):
    """A Dirichlet draw, renormalised: the latent relative genomic abundance."""
    raw = rng.gamma(concentration, 1.0, n_features)
    return raw / raw.sum()


def design_replicate(name, scenario, seed, output, concentration=None):
    """One matrix replicate: counts, lengths, samples, truth and latent expectations."""
    spec = dict(DESIGNS[name])
    concentration = spec.get("concentration", 100.0) if concentration is None else concentration
    n_features = spec["n_features"]
    per_group = spec["samples_per_group"]
    effect = float(str(spec["effect"]).rstrip("x"))
    read_pairs = int(spec["read_pairs"])
    ids = [f"feature_{i + 1:04d}" for i in range(n_features)]
    lengths = np.array([12000 + (i % 7) * 2500 for i in range(n_features)], dtype=float)
    rng = np.random.default_rng(seed)
    # Every genome is the same length here, so relative abundance and read fraction
    # coincide; real genomes with different lengths exercise the two lanes separately.
    baseline = sample_composition(rng, n_features, concentration)
    treated = baseline.copy()
    n_changed = max(2, int(round(n_features * spec["changed_fraction"])))
    chosen = rng.choice(n_features, n_changed, replace=False)
    changed = pd.DataFrame(index=ids)
    if scenario == "spiked":
        up, down = chosen[:n_changed // 2], chosen[n_changed // 2:]
        treated[up] *= effect
        treated[down] /= effect
        # Conserve the changed subset's relative mass, exactly as the read generator does.
        treated[chosen] *= baseline[chosen].sum() / treated[chosen].sum()
    delta = np.log2(treated / baseline)
    changed["feature_id"] = ids
    changed["is_da"] = np.abs(delta) > 1e-10
    changed["true_log2fc"] = delta
    changed["baseline_proportion"] = baseline
    changed["treatment_proportion"] = treated
    changed["length_bp"] = lengths.astype(int)
    samples, counts, expectations = [], {}, []
    for index in range(per_group):
        for group, target in (("Control", baseline), ("Treatment", treated)):
            sample = f"{'C' if group == 'Control' else 'T'}{index + 1:02d}"
            total = max(1, int(round(read_pairs * rng.lognormal(
                -spec["library_cv"]**2 / 2, spec["library_cv"]))))
            if scenario == "null":
                proportions = sample_composition(rng, n_features, concentration)
            else:
                proportions = np.maximum(target * rng.lognormal(-0.02, 0.2, n_features), 0)
                proportions = proportions / proportions.sum()
            drawn = draw_counts(rng, total * proportions, spec["dispersion"])
            for position, (feature, share, count) in enumerate(zip(ids, proportions, drawn)):
                counts.setdefault(feature, {})[sample] = int(count)
                intended = float(target[position]) if scenario == "spiked" else float(baseline[position])
                expectations.append(dict(sample_id=sample, feature_id=feature,
                                         genomic_proportion=share, target_proportion=intended,
                                         read_proportion=share, library_read_pairs=total,
                                         length_bp=float(lengths[position]),
                                         expected_count=float(count)))
            samples.append(dict(sample_id=sample, group=group, read_pairs=total))
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    # pandas builds a dict of dicts with the outer keys as columns, so transpose to
    # put features on the rows and samples on the columns.
    frame = pd.DataFrame(counts).T
    frame.index.name = "feature_id"
    frame = frame.reindex(index=ids, columns=[row["sample_id"] for row in samples])
    frame.to_csv(out / "counts.tsv", sep="\t")
    save_table(pd.DataFrame(samples), out / "samples.tsv")
    pd.DataFrame({"feature_id": ids, "length_bp": lengths.astype(int)}).to_csv(
        out / "features.tsv", sep="\t", index=False)
    changed.to_csv(out / "truth.tsv", sep="\t", index=False)
    save_table(pd.DataFrame(expectations), out / "expectations.tsv")
    pd.DataFrame({"feature_id": ids, "source_id": ids, "match_status": "matched"}).to_csv(
        out / "matching.tsv", sep="\t", index=False)
    return dict(design=name, case=f"{scenario}_s{seed}", scenario=scenario, seed=seed,
                n_features=n_features, samples_per_group=per_group,
                n_changed=int(changed.is_da.sum()), read_pairs=read_pairs,
                dispersion=spec["dispersion"], effect=str(spec["effect"]))


def build_bundle(designs, null_replicates, spiked_replicates, output, first_seed=1000):
    """Materialise a bundle of independent replicates as a matrix input mode bundle.

    Seeds are disjoint across designs and scenarios, so every replicate of every
    design is an independent experiment.
    """
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    listing, summary = [], []
    for position, name in enumerate(designs):
        for scenario, replicates in (("null", null_replicates), ("spiked", spiked_replicates)):
            for replicate in range(replicates):
                seed = first_seed + 1000 * position + replicate + (
                    0 if scenario == "null" else null_replicates)
                case = f"{scenario}_s{seed}"
                record = design_replicate(name=name, scenario=scenario, seed=seed, output=out / case)
                listing.append(dict(case=case, scenario=scenario, seed=seed, design=name))
                summary.append(record)
    save_table(pd.DataFrame(listing), out / "cases.tsv")
    save_table(pd.DataFrame(summary), out / "designs.tsv")
    write_json(out / "designs.json", {"designs": {name: DESIGNS[name] for name in designs},
                                      "null_replicates": null_replicates,
                                      "spiked_replicates": spiked_replicates,
                                      "first_seed": first_seed,
                                      "generator": "Dirichlet-multinomial matrix screening",
                                      "n_cases": len(listing)})
    return summary