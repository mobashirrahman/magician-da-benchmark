"""Analysis preserves condition keys and blocks recommendations from smoke runs."""
from pathlib import Path
import importlib.util
import json
import pandas as pd
import pytest

spec = importlib.util.spec_from_file_location("paper_analysis", Path(__file__).resolve().parents[1] / "tools/paper_analysis.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def test_factor_join_uses_case_identifiers_and_requires_complete_metadata(tmp_path):
    metadata = pd.DataFrame(dict(case=["null_s1", "spiked_s2"], design=["c000", "c000"],
                                 generator=["g1", "g2"], depth=[500000, 2000000]))
    metadata.to_csv(tmp_path / "cases.tsv", sep="\t", index=False)
    scores = pd.DataFrame(dict(case=["spiked_s2", "null_s1", "spiked_s2"], method=["a", "a", "b"]))
    enriched = analysis.enrich_conditions(scores, tmp_path)
    assert enriched.generator.tolist() == ["g2", "g1", "g2"]
    assert enriched.depth.tolist() == [2000000, 500000, 2000000]
    with pytest.raises(ValueError, match="does not cover"):
        analysis.enrich_conditions(pd.DataFrame(dict(case=["null_s9"])), tmp_path)


def test_decision_guide_requires_sufficient_versioned_evidence():
    base = dict(kind="sources", metric="counts", endpoint="read_fraction", method="m",
                generator_version="tier1-v3", generator="g1", depth=500000, n_null=200,
                samples_per_group=20, n_features=300, fraction_da=0.1, log2fc=1.0,
                direction="balanced", total_load="none", zero_structure="sampling", confounder="none",
                n_spiked=100, prereg_eligible=True, fdr=0.02, fdr_ci_high=0.04,
                fdr_ci_low=0.01, source_recall=0.7, null_fpr=0.05,
                null_probability_any_false=0.04, failure_rate=0)
    ranks = pd.DataFrame([base, dict(base, depth=2000000, n_null=20),
                          dict(base, generator_version="original"),
                          dict(base, prereg_eligible=False)])
    guide = analysis.decision_guide(ranks, "tier1")
    assert [g["recommendation"] for g in guide] == ["eligible candidate", "not established", "not established", "not established"]
    assert [g["depth"] for g in guide[:2]] == [500000, 2000000]
    json.dumps(guide, allow_nan=False)
