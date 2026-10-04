"""Shared fixtures: a complete miniature case built through the real code path."""
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from magician.config import load, samples as sample_names
from magician.design import prepare_references, make_design, allocate
from magician.matrices import derive
from magician.truth import build as build_truth
from magician.io import save_table, table


@pytest.fixture
def cfg(tmp_path):
    return load(overrides={"output_dir": str(tmp_path / "out"), "cache_dir": str(tmp_path / "cache"),
                          "genomes": {"synthetic_length": 12000},
                          "experiment": {"samples_per_group": 3, "read_pairs": 2000},
                          "analysis": {"min_total_count": 0, "min_prevalence": 0.0,
                                       "metrics": ["counts", "relative"]}})


@pytest.fixture
def mini(cfg, tmp_path):
    """References, two cases, durable raw counts, derived matrices and truth lanes."""
    refs = tmp_path / "refs"
    prepare_references(cfg, refs)
    cases = {}
    for case in ("null_s42", "spiked_s42"):
        design_dir = tmp_path / case
        make_design(cfg, case, refs, design_dir)
        meta = table(design_dir / "samples.tsv").set_index("sample_id")
        allocations = table(design_dir / "allocations.tsv")
        counts = pd.DataFrame(0, index=sorted(allocations.genome_id.unique()),
                              columns=list(meta.index), dtype=int)
        for sample in meta.index:
            rows = allocations.loc[allocations.sample_id.eq(sample)].sort_values("genome_id")
            values = allocate(int(meta.loc[sample, "read_pairs"]), rows.read_proportion.to_numpy())
            counts[sample] = values
        counts.index.name = "feature_id"
        raw = tmp_path / "raw" / case
        raw.mkdir(parents=True)
        counts.to_csv(raw / "counts.tsv", sep="\t")
        features = table(design_dir / "truth.tsv")[["feature_id", "length_bp"]]
        features.to_csv(raw / "features.tsv", sep="\t", index=False)
        library = pd.DataFrame({"sample_id": list(meta.index),
                                "read_pairs": meta.read_pairs.to_numpy(),
                                "mapped_pairs": counts.sum(axis=0).to_numpy(),
                                "catalogue_pairs": counts.sum(axis=0).to_numpy(),
                                "unassigned_pairs": 0,
                                "unmapped_pairs": 0})
        library.to_csv(raw / "library.tsv", sep="\t", index=False)
        derive(cfg, "sources", raw, tmp_path / "matrices" / case / "sources")
        from magician.config import zero_policies
        lanes = build_truth(cfg, case, "sources", design_dir / "expectations.tsv",
                            design_dir / "samples.tsv", tmp_path / "matrices" / case / "sources" / "features.tsv",
                            tmp_path / "matrices" / case / "sources" / "lanes.tsv",
                            zero_policies=zero_policies(cfg),
                            reference=cfg["analysis"]["truth_reference"],
                            reference_feature=cfg["analysis"]["truth_reference_feature"])
        cases[case] = dict(design=design_dir, raw=raw, matrices=tmp_path / "matrices" / case / "sources",
                           lanes=lanes, truth=table(design_dir / "truth.tsv").set_index("feature_id"),
                           counts=counts)
    return cases


def write_result(path, ids, *, endpoint, effect_scale, log2fc, pvalue, qvalue, zero_policy="none",
                 tested=None, fit_status=None):
    n = len(ids)
    frame = pd.DataFrame({"feature_id": ids, "endpoint": endpoint, "contrast": "Treatment / Control",
                          "effect_scale": effect_scale, "reference": "declared",
                          "adjustment_family": "BH", "zero_policy": zero_policy,
                          "log2fc": log2fc, "pvalue": pvalue, "qvalue": qvalue,
                          "tested": True if tested is None else tested,
                          "fit_status": "fit" if fit_status is None else fit_status,
                          "exclusion_reason": None})
    save_table(frame, path)
    return frame