"""Truth defined on the measurement scale being tested.

A correct read-fraction association is not a false positive when read fraction is
the tested scale, and the current pipeline measures relative genomic abundance.
Every lane here is computed from the latent design, including genome lengths and
the declared denominator, so a method is only ever scored against the hypothesis
it actually tests.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from .io import table, save_table

LANES = ("relative_genomic_abundance", "expected_read_fraction", "clr_log_ratio",
         "reference_relative_change", "absolute_abundance")
# Effect units are compared to the declared units before any fold-change error is
# computed; an unmatched pair is reported as unavailable rather than rescaled.
ENDPOINT_LANE = {"relative_abundance": "relative_genomic_abundance",
                 "read_fraction": "expected_read_fraction",
                 "clr": "clr_log_ratio",
                 "reference_relative": "reference_relative_change",
                 "absolute_abundance": "absolute_abundance"}
TOLERANCE = 1e-10


def read_expectations(path):
    """Latent expectations keyed by feature. Both generator and bundle spellings work."""
    frame = table(path)
    if "feature_id" in frame.columns and "genome_id" not in frame.columns:
        frame = frame.rename(columns={"feature_id": "genome_id"})
    if "genome_id" not in frame.columns:
        raise ValueError(f"{path}: expectations need a genome_id or feature_id column")
    return frame.rename(columns={"genome_id": "feature_id"})


def metric_matrix(expected, metric):
    """Expected values expressed in the units of one measured metric."""
    if metric not in {"counts", "relative", "tpm", "rpkm"}:
        raise ValueError(f"Unknown metric: {metric}")
    dna = expected.read_proportion * expected.length_bp
    denominator = dna.groupby(expected.sample_id).transform("sum")
    return {"counts": expected.expected_count,
            "relative": expected.genomic_proportion,
            "tpm": expected.genomic_proportion * 1e6,
            "rpkm": dna / denominator * 1e9}[metric]


def apply_zero_policy(matrix, library, policy):
    """Declared, unit-aware zero handling.

    `measurement` adds the pseudocount in the metric's own units, which is what a
    fixed 0.5 offset does to proportions and to counts alike. `count` adds it in
    expected-count units, so the same constant means the same thing whatever the
    input scale. `dirichlet` replaces zeros by the posterior mean of a Dirichlet
    fitted to the same sample, a multiplicative zero model.
    """
    kind = policy["zero_handling"]
    if kind == "pseudocount":
        offset = float(policy.get("pseudocount", 0.0))
        if offset < 0:
            raise ValueError("pseudocount must be nonnegative")
        scale = policy.get("pseudocount_scale", "measurement")
        if scale == "measurement":
            return matrix + offset
        if scale != "count":
            raise ValueError(f"Unknown pseudocount_scale: {scale}")
        # Expected counts are matrix * library in every current metric.
        return matrix + offset / library
    if kind == "dirichlet":
        alpha = float(policy.get("dirichlet_alpha", 0.5))
        totals = matrix.sum(axis=1)
        denominator = totals + alpha * matrix.shape[1]
        posterior = (matrix + alpha) / denominator.replace(0, np.nan)
        return matrix.where(matrix > 0, posterior).fillna(0)
    raise ValueError(f"Unknown zero_handling: {kind}")


def centre_log_ratio(matrix):
    """Explicit base-2 CLR: log2 of the transformed values, centred per sample."""
    values = np.log2(matrix.where(matrix > 0))
    return values.subtract(values.mean(axis=1), axis=0)


def _lane(ids, lane, endpoint, effect, tested, labels=None, **extra):
    values = pd.Series(np.asarray(effect, dtype=float), index=pd.Index(ids, name="feature_id"))
    finite = pd.Series(np.isfinite(values.to_numpy()), index=values.index) & pd.Series(tested, index=values.index)
    if labels is None:
        labels = finite & (values.abs() > TOLERANCE)
    else:
        labels = pd.Series(labels, index=values.index).astype(bool) & finite
    frame = pd.DataFrame({"feature_id": ids, "lane": lane, "endpoint": endpoint,
                          "in_tested_set": np.asarray(tested, dtype=bool),
                          "available": finite.to_numpy(),
                          "is_da": labels.to_numpy(),
                          "true_effect": values.to_numpy()})
    for key, value in extra.items():
        frame[key] = value
    return frame


def target_labels(expected, treated, ids):
    """A feature is truly differential when the design's group composition differs.

    Labels come from the target composition, not from the noisy per-sample draws around
    it: under the null the targets are identical and every measured difference is
    finite-sample variation.
    """
    targets = expected.pivot(index="sample_id", columns="feature_id",
                             values="target_proportion").reindex(columns=ids)
    if targets.isna().all().any():
        raise ValueError("expectations.tsv needs target_proportion for truth labels")
    control = targets[~treated].mean(axis=0)
    case_group = targets[treated].mean(axis=0)
    genomic = np.log2(case_group / control)
    dna_control = control * expected.groupby("feature_id").length_bp.first().reindex(ids)
    dna_treated = case_group * expected.groupby("feature_id").length_bp.first().reindex(ids)
    read_fraction = np.log2(dna_treated / dna_control)
    return (genomic.abs() > TOLERANCE), (read_fraction.abs() > TOLERANCE), genomic


def lanes(cfg, case, kind, expected, samples, lengths, retained, reference="median_feature",
          reference_feature=None):
    """Latent truth lanes that do not depend on the input metric or zero policy.

    Relative abundance, expected read fraction and absolute abundance are properties
    of the latent design. The CLR and reference-relative lanes are defined on the
    tested composition and are produced by `clr_lane`.
    """
    expected = expected.copy()
    expected["length_bp"] = expected.feature_id.map(lengths).astype(float)
    ids = sorted(expected.feature_id.unique())
    wide = expected.pivot(index="sample_id", columns="feature_id",
                          values="genomic_proportion").reindex(columns=ids)
    treated = wide.index.to_series().str.startswith("T").to_numpy()
    if not treated.any() or treated.all():
        raise ValueError("Both Control and Treatment samples are required for truth")
    keep = pd.Series(retained, index=pd.Index(ids, name="feature_id")).reindex(ids).fillna(False).astype(bool)
    tested = keep.to_numpy()
    tested_ids = [i for i in ids if keep.get(i, False)]
    if not tested_ids:
        raise ValueError("No features pass filtering, so no truth lane can be evaluated")
    control = wide[~treated].mean(axis=0)
    case_group = wide[treated].mean(axis=0)
    genomic = np.log2(case_group / control)
    da_genomic, da_read, target_genomic = target_labels(expected, treated, ids)
    # The reference-relative hypothesis is a hypothesis about change against the typical
    # feature, so its labels come from the median-centred target differences.
    shifted = target_genomic - target_genomic.median()
    da_reference = (shifted.reindex(ids).abs() > TOLERANCE).to_numpy()
    frames = []

    # Lane 1: relative genomic abundance. The length-corrected share of mapped DNA.
    frames.append(_lane(ids, "relative_genomic_abundance", "relative_abundance",
                        genomic.reindex(ids).to_numpy(), np.ones(len(ids), bool),
                        labels=da_genomic.to_numpy(),
                        effect_scale="log2_relative_abundance", reference="length_corrected_total",
                        denominator="relative genomic abundance of the source or MAG in the sample",
                        input_metric="any", zero_policy="none", reserve_other=False))

    # Lane 2: expected read fraction. The raw count scale is DNA fraction.
    dna = wide * expected.groupby("feature_id").length_bp.first().reindex(ids)
    read_effect = np.log2(dna[treated].mean(axis=0) / dna[~treated].mean(axis=0))
    frames.append(_lane(ids, "expected_read_fraction", "read_fraction",
                        read_effect.reindex(ids).to_numpy(), np.ones(len(ids), bool),
                        labels=da_read.to_numpy(),
                        effect_scale="log2_read_fraction", reference="total_dna_or_read_mass",
                        denominator="expected read or DNA fraction of the source or MAG in the sample",
                        input_metric="any", zero_policy="none", reserve_other=False))

    # Lane 4: change relative to a declared reference or the typical feature.
    typical = float(np.median(genomic.reindex(tested_ids).to_numpy()))
    if reference in {"median_feature", "learned_set"}:
        denominator = "median log2 change over the tested features"
        if reference == "learned_set":
            # ADAPT learns its reference set from the data. The truth counterpart is
            # the typical change over the tested catalogue; a wrongly learned set is a
            # property of the method, not a change of the truth definition.
            denominator += " (counterpart of a learned set)"
        frames.append(_lane(ids, "reference_relative_change", "reference_relative",
                            (genomic - typical).reindex(ids).to_numpy(), tested,
                            labels=da_reference,
                            effect_scale="log2_reference_relative",
                            reference="median_feature" if reference == "median_feature" else "learned_set",
                            denominator=denominator, input_metric="any", zero_policy="none",
                            reserve_other=False))
    elif reference == "declared_reference_feature":
        if reference_feature not in tested_ids:
            raise ValueError(f"Declared reference feature is not tested: {reference_feature}")
        shift = float(genomic[reference_feature])
        frames.append(_lane(ids, "reference_relative_change", "reference_relative",
                            (genomic - shift).reindex(ids).to_numpy(), tested,
                            labels=da_reference,
                            effect_scale="log_ratio_reference", reference=reference_feature,
                            denominator="log ratio to the declared reference feature",
                            input_metric="any", zero_policy="none", reserve_other=False))
    else:
        raise ValueError(f"Unknown truth reference: {reference}")

    # Lane 5: absolute abundance. Only defined when load information exists.
    copies = pd.to_numeric(expected.absolute_copies, errors="coerce")
    copies = pd.DataFrame({"sample_id": expected.sample_id, "feature_id": expected.feature_id,
                           "absolute_copies": copies}).pivot(
        index="sample_id", columns="feature_id", values="absolute_copies")
    copies = copies.reindex(index=wide.index, columns=ids)
    if copies.notna().all().all() and (copies.sum(axis=1) > 0).all():
        absolute = np.log2(copies[treated].mean(axis=0) / copies[~treated].mean(axis=0))
        frames.append(_lane(ids, "absolute_abundance", "absolute_abundance",
                            absolute.reindex(ids).to_numpy(), np.ones(len(ids), bool),
                            labels=da_genomic.to_numpy(),
                            effect_scale="log2_absolute_abundance", reference="total_genome_copies",
                            denominator="expected genome copies per sample",
                            input_metric="any", zero_policy="none", reserve_other=False))
    else:
        frames.append(_lane(ids, "absolute_abundance", "absolute_abundance",
                            np.full(len(ids), np.nan), np.ones(len(ids), bool),
                            effect_scale="log2_absolute_abundance", reference="total_genome_copies",
                            denominator="requires simulated or measured load information",
                            input_metric="any", zero_policy="none", reserve_other=False))
    out = pd.concat(frames, ignore_index=True)
    out.insert(0, "case", case)
    out.insert(1, "kind", kind)
    return out


def clr_lane(cfg, case, kind, expected, samples, lengths, retained, input_metric,
             zero_policy, reserve_other):  # noqa: D401
    """CLR truth in exactly the coordinates the adapter's declared transform produces."""
    expected = expected.copy()
    expected["length_bp"] = expected.feature_id.map(lengths).astype(float)
    ids = sorted(expected.feature_id.unique())
    wide = expected.pivot(index="sample_id", columns="feature_id",
                          values="genomic_proportion").reindex(columns=ids)
    library = samples.set_index("sample_id").read_pairs.reindex(wide.index).astype(float)
    treated = wide.index.to_series().str.startswith("T").to_numpy()
    keep = pd.Series(retained, index=pd.Index(ids, name="feature_id")).reindex(ids).fillna(False).astype(bool)
    tested_ids = [i for i in ids if keep.get(i, False)]
    if not tested_ids:
        raise ValueError("No features pass filtering, so no CLR lane can be evaluated")
    columns = tested_ids + (["other"] if reserve_other else [])
    series = metric_matrix(expected, input_metric)
    series.index = pd.MultiIndex.from_arrays([expected.sample_id.to_numpy(),
                                              expected.feature_id.to_numpy()])
    values = series.unstack().reindex(index=wide.index, columns=columns)
    if reserve_other:
        totals = expected.groupby("sample_id").expected_count.sum()
        tested_total = expected[expected.feature_id.isin(tested_ids)].groupby("sample_id").expected_count.sum()
        values["other"] = (totals - tested_total).reindex(values.index).fillna(0)
    transformed = centre_log_ratio(apply_zero_policy(values, library.reindex(values.index), zero_policy))
    effect = (transformed[treated].mean(axis=0) - transformed[~treated].mean(axis=0)).reindex(ids)
    da_genomic, _, _ = target_labels(expected, treated, ids)
    frame = _lane(ids, "clr_log_ratio", "clr", effect.to_numpy(), keep.to_numpy(),
                  labels=da_genomic.to_numpy(),
                  effect_scale="log2_clr", reference="geometric_mean_of_tested_features",
                  denominator=("geometric mean of log2 values over the tested features per sample"
                               + (", including the reserved other category" if reserve_other else "")),
                  input_metric=input_metric, zero_policy=zero_policy.get("name", "none"),
                  reserve_other=bool(reserve_other))
    frame.insert(0, "case", case)
    frame.insert(1, "kind", kind)
    return frame


def mag_expectations(expected, samples, features, matching):
    """Latent expectations for MAGs, inherited from their assigned sources.

    A MAG recovers a fraction of its source's length at roughly uniform coverage, so
    its length-corrected relative abundance equals its source's. Its expected counts
    scale with its recovered length, and the composition it participates in is the
    tested MAG catalogue, not the source list. Unmatched and ambiguous MAGs have no
    source expectations and are left out; evaluation already excludes them.
    """
    matched = matching.loc[matching.match_status.eq("matched"),
                           ["feature_id", "source_id"]].drop_duplicates()
    if matched.empty:
        raise ValueError("No matched MAGs, so no MAG truth lane can be evaluated")
    frame = expected.merge(matched.rename(columns={"feature_id": "mag_id", "source_id": "feature_id"}),
                           on="feature_id", how="inner", validate="many_to_many")
    mag_lengths = features.length_bp.reindex(matched.feature_id)
    source_lengths = frame.groupby("feature_id").length_bp.first()
    share = (mag_lengths.reindex(frame.mag_id).to_numpy(dtype=float)
             / source_lengths.reindex(frame.feature_id).to_numpy(dtype=float))
    frame["mag_id"] = frame.mag_id
    frame["feature_id"] = frame.mag_id
    frame["length_bp"] = mag_lengths.reindex(frame.feature_id).to_numpy(dtype=float)
    frame["expected_count"] = frame.expected_count * share
    frame["read_proportion"] = frame.read_proportion * share
    totals = frame.groupby("sample_id").read_proportion.transform("sum")
    frame["read_proportion"] = (frame.read_proportion / totals.where(totals.ne(0))).fillna(0)
    frame["absolute_copies"] = (pd.to_numeric(frame.absolute_copies, errors="coerce")
                                * share)
    return frame.drop(columns=["mag_id"])


def build(cfg, case, kind, expectations_file, samples_file, features_file, output,
          zero_policies=None, reference="median_feature", reference_feature=None,
          matching_file=None):
    """Write the lane table for one catalogue, from latent expectations and filtering.

    Latent lanes are emitted once. CLR lanes are emitted for every declared input
    metric and zero policy, plus the reserved-category variant, because each of those
    is a different geometry and therefore a different truth.
    """
    expected = read_expectations(expectations_file)
    samples = table(samples_file)
    features = table(features_file).set_index("feature_id")
    if not features.index.is_unique:
        raise ValueError("Feature IDs are not unique")
    expected["library_read_pairs"] = expected.sample_id.map(samples.set_index("sample_id").read_pairs)
    expected["expected_count"] = expected.read_proportion * expected.library_read_pairs
    if "absolute_copies" not in expected:
        expected["absolute_copies"] = np.nan
    if matching_file is not None and kind == "mags":
        expected = mag_expectations(expected, samples, features, table(matching_file))
    ids = expected.feature_id.unique()
    retained = features.reindex(ids).retained_for_da.astype(str).str.lower().eq("true")
    # Lengths come from the measured catalogue; the design supplies the rest, so a
    # filtered feature still has the length its read fraction needs.
    lengths = features.length_bp.reindex(ids)
    if "length_bp" in expected:
        lengths = lengths.fillna(expected.groupby("feature_id").length_bp.first().reindex(ids))
    # The untransformed geometry is always emitted: it is what a package that centres
    # raw counts itself, such as ALDEx with gamma zero, is scored against.
    policies = [{"name": "none", "zero_handling": "pseudocount", "pseudocount": 0.0}] + list(zero_policies or [])
    frames = [lanes(cfg, case, kind, expected, samples, lengths, retained,
                    reference=reference, reference_feature=reference_feature)]
    for metric in cfg["analysis"]["metrics"]:
        for policy in policies:
            for other in (False, True):
                try:
                    frames.append(clr_lane(cfg, case, kind, expected, samples, lengths, retained,
                                           metric, policy, other))
                except ValueError as exc:
                    if other:
                        # No reserved category exists when nothing was filtered out.
                        continue
                    raise ValueError(f"Cannot compute CLR truth for {metric}/{policy['name']}: {exc}") from None
    out = pd.concat(frames, ignore_index=True)
    save_table(out, Path(output))
    return out