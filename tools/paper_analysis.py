#!/usr/bin/env python3
"""Condition-specific summaries from benchmark archives.

  python tools/paper_analysis.py --tier1 results_paper_tier1_v2/benchmark --out analysis_v2

Rank uncertainty and contrasts compare methods only within matching conditions,
endpoints, feature-table kinds and input metrics. Variance shares are descriptive,
order-dependent summaries, not a fitted mixed-effects model. Consensus is a
specified rule only; this script does not compute or validate consensus calls.
"""
from pathlib import Path
import argparse
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from magician.variance import variance_decomposition, pairwise_contrasts, rank_uncertainty
from magician.evaluation import CONDITION_COLUMNS
from magician.generators import GENERATOR_VERSION
from magician.io import save_table, write_json, read_json

CONSENSUS_K = 3


def load_scores(benchmark):
    return pd.read_csv(Path(benchmark) / 'scores.tsv', sep='\t', keep_default_na=False, na_values=['NA'])


def enrich_conditions(scores, bundle):
    """Recover factors using an explicit case key; never join on reused cell ids."""
    metadata = pd.read_csv(Path(bundle) / 'cases.tsv', sep='\t', keep_default_na=False)
    if not metadata.case.is_unique:
        raise ValueError('Bundle metadata contains duplicate cases')
    missing = set(scores.case) - set(metadata.case)
    if missing:
        raise ValueError(f'Bundle metadata does not cover {len(missing)} score cases')
    columns = [c for c in CONDITION_COLUMNS if c in metadata and c not in scores]
    return scores.merge(metadata[['case', *columns]], on='case', how='left', validate='many_to_one')


def decision_guide(ranking, tier):
    """Use the benchmark's single implementation of predeclared eligibility."""
    columns = ['kind', 'metric', 'endpoint', 'method', *CONDITION_COLUMNS]
    out = []
    for _, row in ranking.iterrows():
        condition = {c: (None if pd.isna(row[c]) else row[c].item() if hasattr(row[c], 'item') else row[c])
                     for c in columns if c in ranking}
        enough = int(row.n_null) >= 100 and int(row.n_spiked) >= 50
        versioned = (condition.get('generator_version') == GENERATOR_VERSION
                     and all(condition.get(c) is not None for c in CONDITION_COLUMNS))
        eligible = str(row.get('prereg_eligible', False)).lower() == 'true'
        record = dict(condition, tier=tier, n_null=int(row.n_null), n_spiked=int(row.n_spiked),
                      evidence='repeated_simulation' if enough and versioned else 'exploratory',
                      recommendation='eligible candidate' if enough and versioned and eligible else 'not established')
        record['null_conservative'] = str(row.get('null_conservative', False)).lower() == 'true'
        for c in ('fdr', 'fdr_ci_low', 'fdr_ci_high', 'source_recall', 'null_fpr',
                  'null_probability_any_false', 'failure_rate', 'endpoint_fdr',
                  'endpoint_recall', 'spillover_share'):
            v = row.get(c)
            record[c] = float(v) if v is not None and pd.notna(v) else None
        out.append(record)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tier1')
    p.add_argument('--tier1-bundle')
    p.add_argument('--tier2')
    p.add_argument('--tier3a')
    p.add_argument('--out', default='analysis_v2')
    p.add_argument('--bootstrap-draws', type=int, default=2000)
    args = p.parse_args()
    if args.bootstrap_draws < 100:
        p.error('bootstrap-draws must be >=100')
    output = Path(args.out)
    (output / 'figures').mkdir(parents=True, exist_ok=True)
    frames, guide = [], []
    for tier, path in [('tier1', args.tier1), ('tier2', args.tier2), ('tier3a', args.tier3a)]:
        if not path:
            continue
        scores = load_scores(path)
        if tier == 'tier1':
            bundle = args.tier1_bundle
            provenance = Path(path).parent / 'provenance/matrix_bundle.json'
            if bundle is None and provenance.exists():
                bundle = read_json(provenance)['bundle']
            if bundle:
                scores = enrich_conditions(scores, bundle)
        scores['tier'] = tier
        frames.append(scores)
        rank_path = Path(path) / 'ranking.tsv'
        if rank_path.exists():
            ranking = pd.read_csv(rank_path, sep='\t', keep_default_na=False, na_values=['NA'])
            guide.extend(decision_guide(ranking, tier))
    if not frames:
        p.error('Supply a completed benchmark directory')
    scores = pd.concat(frames, ignore_index=True)
    write_json(output / 'decision_guide.json', guide)
    scope = ['tier', 'kind', 'metric', 'endpoint']
    scope += [c for c in CONDITION_COLUMNS if c in scores]
    # No comparison of different endpoints, input metrics or design conditions.
    spiked = scores.loc[scores.scenario.eq('spiked') & scores.status.eq('success')].copy()
    contrasts, ranks, variance_rows = [], [], []
    for values, group in spiked.groupby(scope, sort=True, dropna=False):
        condition = dict(zip(scope, values))
        c = pairwise_contrasts(group, 'method', 'fdr', draws=args.bootstrap_draws)
        if len(c):
            for key, value in condition.items():
                c[key] = value
            contrasts.append(c)
        for key, interval in rank_uncertainty(group, ('method',), 'source_recall', True,
                                             draws=args.bootstrap_draws).items():
            ranks.append(dict(condition, method=key[0], median_rank=interval[0],
                              low=interval[1], high=interval[2]))
    for values, group in spiked.groupby(['tier', 'kind', 'metric', 'endpoint'], sort=True):
        context = dict(zip(['tier', 'kind', 'metric', 'endpoint'], values))
        factors = [c for c in ['method', 'generator', 'depth', 'fraction_da', 'direction'] if c in group]
        for statistic in ('fdr', 'source_recall', 'spillover_share'):
            variance_rows.append(dict(context, statistic=statistic,
                interpretation='descriptive; sequential factor order matters',
                **variance_decomposition(group, factors, statistic)))
    save_table(pd.DataFrame(variance_rows), output / 'variance.tsv')
    save_table(pd.concat(contrasts, ignore_index=True) if contrasts else pd.DataFrame(), output / 'contrasts.tsv')
    save_table(pd.DataFrame(ranks), output / 'rank_uncertainty.tsv')
    write_json(output / 'consensus.json', dict(k=CONSENSUS_K, status='specified_only_not_evaluated'))
    (output / 'report.md').write_text(
        '# Conditional benchmark analysis\n\n'
        f'Cases: {scores.case.nunique()}; method configurations: {scores.method.nunique()}.\n\n'
        'Discoveries are scored against the implanted design truth, the same set for every method.\n'
        'Endpoint-lane scores are secondary; endpoint FDR is unavailable where a scale has no null feature,\n'
        'and spillover_share is the share of discoveries that are closure shifts outside the design.\n'
        'Method comparisons are separated by endpoint, metric, table kind and recorded design conditions.\n'
        'Eligibility needs an upper FDR bound <=0.10, null p-value rejection rate <=0.07, failure rate <5%,\n'
        'and >=100 null and >=50 spiked runs per condition; null runs are shared across effect factors.\n'
        'Thin validation/fractional cells remain exploratory. G2 resamples measured profiles; G1 and G3 are\n'
        'synthetic models calibrated on the same donor table.\n'
        'Variance shares depend on factor order and are not a fitted mixed-effects model. Contrasts against the observed best method are exploratory.\n'
        'Batch confounding and consensus performance have not been evaluated.\n')
    print(f'Wrote conditional analysis to {output}', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
