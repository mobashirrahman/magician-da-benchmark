"""Reference-based table: mapping to a deliberately incomplete reference.

Emulates the usual database-gap situation with 70% of genomes present.
Reads from absent genomes become unmapped (library unmapped increases);
counts for present genomes are unchanged. The present set is fixed per run
(from the manifest + seed) so every case shares one incomplete database.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from .io import table, save_table, write_json


def select_present(references_file, fraction_present, seed, output):
    """Choose the present subset; writes reference_genomes.tsv + json."""
    refs = table(references_file)
    ids = sorted(refs.genome_id.astype(str).tolist()) if len(refs) else []
    n_present = int(round(len(ids) * float(fraction_present)))
    rng = np.random.default_rng(int(seed))
    present = set(rng.choice(ids, n_present, replace=False).tolist()) if ids else set()
    frame = pd.DataFrame({"genome_id": ids, "present": [g in present for g in ids]})
    save_table(frame, output)
    write_json(Path(str(output)).with_suffix(".json"),
               {"fraction_present": float(fraction_present), "seed": int(seed),
                "n_total": len(ids), "n_present": len(present)})
    return frame


def build_raw_reference(raw_sources, reference_genomes_file, output):
    """Derive raw_reference from raw_sources by dropping absent genomes."""
    raw_sources, out = Path(raw_sources), Path(output)
    out.mkdir(parents=True, exist_ok=True)
    counts = table(raw_sources / "counts.tsv").set_index("feature_id")
    features = table(raw_sources / "features.tsv").set_index("feature_id")
    library = table(raw_sources / "library.tsv").set_index("sample_id")
    present = table(reference_genomes_file)
    keep = present.loc[present.present.astype(str).str.lower().eq("true"), "genome_id"].tolist()
    keep = [g for g in keep if g in counts.index]
    dropped = counts.loc[[i for i in counts.index if i not in set(keep)]].sum(axis=0) if len(counts) else 0
    counts = counts.reindex(keep) if keep else counts.iloc[0:0]
    counts.index.name = "feature_id"
    features = features.reindex(counts.index)
    features.index.name = "feature_id"
    # Reads from absent genomes leave the catalogue: mapped/catalogue shrink.
    lib = library.copy()
    lib["catalogue_pairs"] = counts.sum(axis=0).reindex(lib.index).fillna(0).astype(int)
    lib["mapped_pairs"] = lib["catalogue_pairs"]
    lib["unassigned_pairs"] = 0
    lib["unmapped_pairs"] = (lib["read_pairs"] - lib["mapped_pairs"]).clip(lower=0)
    save_table(counts.reset_index(), out / "counts.tsv")
    save_table(features.reset_index(), out / "features.tsv")
    save_table(lib.reset_index(), out / "library.tsv")
    write_json(out / "raw.json", {"kind": "reference", "n_features": int(len(counts)),
                                  "n_samples": int(len(lib)),
                                  "count_unit": "properly paired primary first-mate fragments",
                                  "raw_retained": True, "derived_from": "raw_sources",
                                  "dropped_pairs_per_sample": {str(k): int(v) for k, v in dropped.items()}
                                  if hasattr(dropped, "items") else {}})
    return out


def build_catalogue(references_file, reference_genomes_file, features_output, matching_output):
    """Reference catalogue features + matching from the present subset."""
    refs = table(references_file)
    present = table(reference_genomes_file)
    keep = present.loc[present.present.astype(str).str.lower().eq("true"), "genome_id"].tolist()
    # references.tsv has genome_id,length_bp (from design.prepare_references genomes.tsv).
    if "length_bp" not in refs.columns:
        raise ValueError("references file needs length_bp")
    feats = refs.loc[refs.genome_id.isin(keep), ["genome_id", "length_bp"]].copy()
    feats = feats.rename(columns={"genome_id": "feature_id"})
    save_table(feats, features_output)
    reference_matching(reference_genomes_file, matching_output)
    return features_output


def filter_counts(source_counts_file, reference_genomes_file, output):
    """Filter one per-sample source count file to the present subset."""
    counts = table(source_counts_file)
    present = table(reference_genomes_file)
    keep = set(present.loc[present.present.astype(str).str.lower().eq("true"), "genome_id"].tolist())
    # Source count files have feature_id,count (and length_bp for contigs files).
    filtered = counts.loc[counts.feature_id.isin(keep)].copy()
    save_table(filtered, output)
    return output


def reference_matching(reference_genomes_file, output):
    """Matching table for the reference catalogue: present genomes, 1:1 matched.

    Only present genomes are listed (as both feature and source); absent DA
    genomes count as misses via source recall, exactly like unrecovered MAGs.
    """
    present = table(reference_genomes_file)
    rows = []
    for row in present.itertuples():
        if str(row.present).lower() != "true":
            continue
        rows.append(dict(feature_id=row.genome_id, source_id=row.genome_id,
                         match_status="matched",
                         aligned_fraction=1.0, identity=1.0,
                         source_completeness=1.0))
    save_table(pd.DataFrame(rows, columns=["feature_id", "source_id", "match_status",
                                           "aligned_fraction", "identity",
                                           "source_completeness"]), output)
    return output
