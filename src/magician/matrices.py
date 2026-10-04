"""Durable raw quantification, then cheap derivation of filtered and normalized inputs.

Raw counts, feature lengths and library totals are written once and kept. Filtering
and normalization are derived from them in a rule that costs nothing to repeat, so
changing a threshold no longer requires deleted upstream counts.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from .io import table, save_table, write_json

OTHER = "other__filtered_reads"


def normalized_metrics(counts, lengths, libraries):
    per_base = counts.divide(lengths, axis=0)
    # Relative genomic abundance is length corrected, rather than read fraction.
    denom = per_base.sum(axis=0).astype(float)
    denom = denom.where(denom.ne(0))
    return {"counts": counts,
            "tpm": per_base.divide(denom, axis=1).fillna(0) * 1e6,
            "relative": per_base.divide(denom, axis=1).fillna(0),
            "rpkm": per_base.divide(libraries, axis=1).fillna(0) * 1e9}


def _read_counts(path):
    df = table(path)
    if not df.feature_id.is_unique:
        raise ValueError(f"Duplicate feature IDs in {path}")
    return df.set_index("feature_id")["count"]


def quantify(cfg, kind, count_files, samples_file, features_file, assignments_file, output):
    """Persist raw counts, lengths and library totals. Nothing here filters."""
    meta = table(samples_file).set_index("sample_id")
    samples = list(meta.index)
    features = table(features_file).rename(columns={"genome_id": "feature_id"}).set_index("feature_id")
    if not features.index.is_unique:
        raise ValueError("Feature IDs are not unique")
    counts = pd.DataFrame(0, index=features.index, columns=samples, dtype=int)
    assignment = table(assignments_file).set_index("contig_id").feature_id.to_dict() if kind == "mags" else None
    if len(count_files) != len(samples):
        raise ValueError("Missing sample count files")
    library = pd.DataFrame({"sample_id": samples,
                            "read_pairs": meta.read_pairs.astype(int),
                            "mapped_pairs": 0, "catalogue_pairs": 0,
                            "unassigned_pairs": 0, "unmapped_pairs": 0})
    library = library.set_index("sample_id")
    for sample, path in zip(samples, count_files):
        values = _read_counts(path)
        if kind == "mags":
            contig_total = int(values.sum())
            mapped = values.to_frame("count")
            mapped["mag"] = mapped.index.map(assignment)
            totals = mapped.loc[mapped["mag"].notna()].groupby("mag")["count"].sum()
            library.loc[sample, "catalogue_pairs"] = int(totals.sum())
            library.loc[sample, "mapped_pairs"] = contig_total
            library.loc[sample, "unassigned_pairs"] = contig_total - int(totals.sum())
        else:
            totals = values
            library.loc[sample, "catalogue_pairs"] = int(totals.sum())
            library.loc[sample, "mapped_pairs"] = int(totals.sum())
        if set(totals.index) - set(counts.index):
            raise ValueError("Counts contain unknown feature IDs")
        counts.loc[totals.index, sample] = totals.astype(int)
    library["unmapped_pairs"] = library.read_pairs - library.mapped_pairs
    if (library[["unassigned_pairs", "unmapped_pairs"]] < 0).any().any():
        raise ValueError("Library accounting is negative; mapped pairs exceed simulated pairs")
    counts.index.name = "feature_id"
    features.index.name = "feature_id"
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    save_table(counts.reset_index(), out / "counts.tsv")
    save_table(features.reset_index(), out / "features.tsv")
    save_table(library.reset_index(), out / "library.tsv")
    write_json(out / "raw.json", {"kind": kind, "n_features": int(len(features)),
                                  "n_samples": len(samples),
                                  "count_unit": "properly paired primary first-mate fragments",
                                  "raw_retained": True})


def derive(cfg, kind, raw, output):
    """Filtering, reserved-category bookkeeping and normalization, from raw counts only."""
    raw = Path(raw)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    counts = table(raw / "counts.tsv").set_index("feature_id")
    features = table(raw / "features.tsv").set_index("feature_id")
    library = table(raw / "library.tsv").set_index("sample_id")
    if not counts.index.is_unique or set(counts.columns) != set(library.index):
        raise ValueError("Raw quantification disagrees with itself")
    counts = counts.reindex(features.index).fillna(0)
    counts.index.name = "feature_id"
    analysis = cfg["analysis"]
    retained = ((counts.sum(axis=1) >= analysis["min_total_count"])
                & ((counts > 0).mean(axis=1) >= analysis["min_prevalence"]))
    features["retained_for_da"] = retained.astype(bool)
    features["reserved_other"] = False
    reserve = bool(analysis["reserve_other"])
    if reserve and (~retained).any():
        # Filtered counts leave the composition as one `other` category so that library
        # totals are preserved. This is a declared variant: it changes the reference
        # geometry of every CLR method, and the truth lanes are computed in the matching
        # coordinates. It is never a tested biological discovery.
        totals = counts.loc[~retained].sum(axis=0)
        counts = counts.loc[retained].copy()
        features = features.loc[retained].copy()
        row = pd.DataFrame({column: [np.nan] for column in features.columns},
                           index=pd.Index([OTHER]))
        row["length_bp"] = 0
        row["retained_for_da"] = False
        row["reserved_other"] = True
        features = pd.concat([features, row])
        counts.loc[OTHER] = totals
    metrics = normalized_metrics(counts, features.length_bp.astype(float).replace(0, np.nan),
                                library.read_pairs)
    if reserve and (~retained).any():
        # The reserved category is a mass, not a genome: keep the totals of the
        # features it collects rather than dividing by its zero length.
        totals = counts.loc[OTHER]
        for name, scale in (("tpm", 1e6), ("rpkm", 1e9)):
            denominator = metrics[name].drop(index=OTHER).sum(axis=0).astype(float)
            denominator = denominator.where(denominator.ne(0))
            metrics[name].loc[OTHER] = totals.divide(denominator).fillna(0) * scale
        grand = counts.sum(axis=0).astype(float)
        metrics["relative"].loc[OTHER] = (totals / grand.where(grand.ne(0))).fillna(0)
    features.index.name = "feature_id"
    features = features.reset_index()
    save_table(features, out / "features.tsv")
    for metric in analysis["metrics"]:
        metrics[metric].to_csv(out / f"{metric}.tsv", sep="\t", na_rep="NA")
    save_table(library.reset_index(), out / "library.tsv")
    write_json(out / "matrix.json", {"kind": kind, "n_features": int(len(counts)),
                                     "n_retained": int(retained.sum()),
                                     "reserved_other": OTHER if reserve else None,
                                     "n_samples": int(counts.shape[1]),
                                     "status": "ready" if len(counts) else "no_mags",
                                     "count_unit": "properly paired primary first-mate fragments",
                                     "relative_unit": "length-corrected fraction of mapped genomic abundance",
                                     "rpkm_denominator": "total simulated read pairs"})


def build(cfg, kind, count_files, samples_file, features_file, assignments_file, output):
    """Quantify and derive in one call, for callers that do not split the two rules."""
    raw = Path(output) / "raw"
    quantify(cfg, kind, count_files, samples_file, features_file, assignments_file, raw)
    derive(cfg, kind, raw, output)