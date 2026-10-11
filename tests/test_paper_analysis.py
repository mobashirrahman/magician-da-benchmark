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


def _ranking():
    """One stratum across two generators, each method with a different eligibility profile."""
    profiles = {'always': (True, True), 'partial': (True, False), 'never': (False, False)}
    rows = []
    for index, generator in enumerate(('g1', 'g2')):
        for method, eligible in profiles.items():
            rows.append(dict(kind='sources', metric='counts', endpoint='clr', method=method,
                             generator=generator, fdr=0.01 + 0.04 * len(rows),
                             source_recall=0.2 + 0.1 * len(rows), null_fpr=0.02 + 0.01 * len(rows),
                             spillover_share=0.01, prereg_eligible=eligible[index]))
    return pd.DataFrame(rows)


def test_cell_summary_keeps_methods_that_are_never_eligible():
    summary = analysis._cell_summary(_ranking(), ['fdr', 'null_fpr'])
    counts = summary.set_index('method')[['n_cells', 'n_eligible']]
    assert counts.to_dict('index') == {'always': {'n_cells': 2, 'n_eligible': 2},
                                       'partial': {'n_cells': 2, 'n_eligible': 1},
                                       'never': {'n_cells': 2, 'n_eligible': 0}}


def test_write_figures_renders_every_panel_and_leaves_no_empty_directory(tmp_path):
    pytest.importorskip('matplotlib')
    ranking = _ranking()
    variance = pd.DataFrame([dict(kind='sources', metric='counts', endpoint='clr',
                                  statistic='fdr', method=0.3, generator=0.2, depth=0.1,
                                  fraction_da=0.1, direction=0.1, residual=0.2)])
    ranks = pd.DataFrame([dict(kind='sources', metric='counts', endpoint='clr', method=m,
                               median_rank=r, low=r - 0.5, high=r + 0.5)
                          for m, r in (('strict', 1.0), ('partial', 2.0))])
    written = analysis.write_figures(ranking, variance, ranks, tmp_path)
    assert written == ['f1_calibration_vs_sensitivity.png', 'f2_null_calibration.png',
                       'f3_variance_shares.png', 'f4_rank_uncertainty.png',
                       'f5_closure_spillover.png']
    assert sorted(p.name for p in tmp_path.glob('*.png')) == written
    assert all((tmp_path / name.replace('.png', '.svg')).exists() for name in written)
    assert all((tmp_path / name).stat().st_size > 0 for name in written)


def test_write_figures_honours_the_requested_formats(tmp_path):
    pytest.importorskip('matplotlib')
    written = analysis.write_figures(_ranking(), pd.DataFrame(), pd.DataFrame(), tmp_path,
                                     formats=('png',))
    assert written
    assert not list(tmp_path.glob('*.svg'))
    assert all((tmp_path / name).exists() for name in written)


def test_write_figures_survives_an_empty_ranking(tmp_path):
    pytest.importorskip('matplotlib')
    empty = pd.DataFrame(columns=['kind', 'metric', 'endpoint', 'method', 'fdr', 'source_recall',
                                  'null_fpr', 'spillover_share', 'prereg_eligible'])
    assert analysis.write_figures(empty, pd.DataFrame(), pd.DataFrame(), tmp_path) == []
