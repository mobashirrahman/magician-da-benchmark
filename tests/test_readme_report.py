"""The README's measured section is generated, and the generator is honest about it."""
from pathlib import Path
import importlib.util
import json

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results_paper_tier1_v3'
BENCHMARK = RESULTS / 'benchmark'


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'tools' / f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


report = _load('readme_report')


def _ranking():
    """Three conditions for two methods so medians, eligibility and ranks are all defined."""
    rows = []
    for generator in ('g1', 'g2', 'g3'):
        for method, eligible in (('limma_voom', True), ('deseq2', False)):
            index = len(rows)
            rows.append(dict(kind='sources', metric='counts', endpoint='read_fraction',
                             method=method, generator=generator, n_failed=0, n_empty=0,
                             n_timeout=0, null_fpr=0.05 + 0.01 * index, fdr=0.04 * index,
                             source_recall=0.1 * index, spillover_share=0.01,
                             prereg_eligible=eligible, rank='1' if index else '2'))
    return pd.DataFrame(rows)


def test_splice_replaces_only_the_generated_block():
    readme = f"# Title\n\nhand written\n\n{report.BEGIN}\nold\n{report.END}\n\ntail\n"
    updated = report.splice(readme, f'{report.BEGIN}\nnew\n{report.END}')
    assert updated == f"# Title\n\nhand written\n\n{report.BEGIN}\nnew\n{report.END}\n\ntail\n"
    assert report.splice(updated, f'{report.BEGIN}\nnew\n{report.END}') == updated


def test_splice_refuses_a_readme_without_markers():
    with pytest.raises(SystemExit, match='missing the generated-block markers'):
        report.splice('# Title only\n', 'block')


def test_label_and_span_keep_registry_ids_readable():
    assert report.label('deseq2') == 'DESeq2'
    assert report.label('aldex2_g0') == report.label('aldex3_g0') == 'ALDEx2/3 (gamma 0)'
    assert report.label('unregistered_method') == 'unregistered_method'
    assert report.span([0.5]) == '0.500'
    assert report.span([0.5, 0.25]) == '0.250-0.500'


def test_verdict_numbers_are_derived_not_written():
    frame = _ranking()
    counts = frame[frame.metric.eq('counts')]
    deseq2 = counts[counts.method.eq('deseq2')]
    facts = dict(ranking_rows=frame, spillover_max=float(frame.spillover_share.max()),
                 winners={('read_fraction', 'counts'): (3, pd.Series({'limma_voom': 2}))})
    bullets = report.verdict(facts)
    assert all(isinstance(b, str) and b.startswith('- ') for b in bullets)
    assert any(f'{deseq2.null_fpr.median():.2f}' in b for b in bullets)
    assert any(f'{deseq2.source_recall.median():.2f}' in b for b in bullets)
    assert any(f'eligible in 0 of {len(deseq2)}' in b for b in bullets)
    assert any(f'{frame.spillover_share.max():.3f}' in b for b in bullets)


def test_method_table_merges_registry_variants_into_one_row(tmp_path):
    frame = _ranking()
    facts = dict(per_method=pd.DataFrame([dict(metric='counts', endpoint='clr', method=m,
                                              cells=2, eligible=1, null_fpr=0.02,
                                              fdr=0.01, recall=0.3)
                                         for m in ('aldex2', 'aldex3')]))
    rows = report.render_method_table(facts)
    assert len(rows) == 3
    assert rows[2].count('|') == 8
    assert 'ALDEx2/3 (gamma 0.5)' in rows[2]


def test_rendered_block_is_deterministic_and_self_describing(tmp_path):
    facts = dict(results='results_x', datasets=10, methods=3, rows=9, rows_eligible=4,
                 failed=0, empty=0, timeout=0, steps='9 of 9', job_seconds=3600.0,
                 retained_gb=1.0, budget_gb=200.0, spillover_max=0.01,
                 per_method=pd.DataFrame([dict(metric='counts', endpoint='clr', method='linda',
                                              cells=2, eligible=1, null_fpr=0.02, fdr=0.01,
                                              recall=0.3)]),
                 winners={('clr', 'counts'): (1, pd.Series({'linda': 1}))},
                 ranking_rows=_ranking())
    first = report.render_block(facts, tmp_path)
    assert first == report.render_block(facts, tmp_path)
    assert first.startswith(report.BEGIN) and first.rstrip().endswith(report.END)
    assert 'tools/readme_report.py' in first
    assert '| Datasets evaluated | 10 |' in first


@pytest.mark.skipif(not (BENCHMARK / 'ranking.tsv').exists(),
                    reason='no completed Tier 1 run on this machine')
def test_committed_readme_matches_the_committed_run():
    readme = (ROOT / 'README.md').read_text()
    block = report.render_block(report.facts(RESULTS), ROOT / 'docs/figures/tier1_v3')
    assert report.splice(readme, block) == readme


@pytest.mark.skipif(not (BENCHMARK / 'ranking.tsv').exists(),
                    reason='no completed Tier 1 run on this machine')
def test_readme_figures_exist_and_are_not_empty():
    text = (ROOT / 'README.md').read_text()
    block = text.split(report.BEGIN)[1].split(report.END)[0]
    referenced = [line.split('](')[1].rstrip(')') for line in block.splitlines()
                  if line.startswith('![')]
    assert referenced
    for rel in referenced:
        path = ROOT / rel
        assert path.exists(), rel
        assert path.stat().st_size > 0, rel
    json.dumps(referenced)