#!/usr/bin/env python3
"""Condition-specific summaries from benchmark archives.

  python tools/paper_analysis.py --tier1 results_paper_tier1_v2/benchmark --out analysis_v2

Rank uncertainty and contrasts compare methods only within matching conditions,
endpoints, feature-table kinds and input metrics. Variance shares are descriptive,
order-dependent summaries, not a fitted mixed-effects model. Consensus is a
specified rule only; this script does not compute or validate consensus calls.

Figures pool the cells of one endpoint and input metric so that a panel never mixes
endpoints, input metrics or table kinds. Pooling across design conditions is
descriptive: a marker is filled when a method is eligible in every condition it was
evaluated in, hollow when that eligibility is partial and grey when it is absent.
"""
from pathlib import Path
import argparse
import os
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from magician.variance import variance_decomposition, pairwise_contrasts, rank_uncertainty
from magician.evaluation import CONDITION_COLUMNS
from magician.generators import GENERATOR_VERSION
from magician.io import save_table, write_json, read_json

CONSENSUS_K = 3
NULL_LIMIT = 0.07
NULL_CONSERVATIVE = 0.03
STRATUM = ['kind', 'metric', 'endpoint']


def _pyplot(config_dir):
    """Import pyplot headless with a writable config directory outside the tracked tree."""
    os.environ.setdefault('MPLCONFIGDIR', str(Path(config_dir)))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _save_factory(plt, out, formats=('png', 'svg')):
    """Return a write-then-close callback so figures and formats stay in one place."""
    def save(fig, stem):
        for extension in formats:
            fig.savefig(Path(out) / f'{stem}.{extension}',
                        **({'dpi': 160} if extension == 'png' else {}))
        plt.close(fig)
        return f'{stem}.png'
    save.plt = plt
    return save


def _palette(methods):
    """One colour per method, stable across every figure of a run."""
    return {name: f'C{i % 20}' for i, name in enumerate(sorted(set(methods)))}


def _mark(ax, x, y, colour, eligible, cells, size=48):
    """Filled when eligible everywhere, hollow when partly, grey when never."""
    if cells and eligible >= cells:
        ax.scatter(x, y, color=colour, s=size, zorder=3)
    elif eligible:
        ax.scatter(x, y, facecolors='white', edgecolors=colour, linewidths=1.6, s=size, zorder=3)
    else:
        ax.scatter(x, y, color='0.78', s=size * 0.7, zorder=2)


def _cell_summary(frame, statistics):
    """One row per method and stratum: pooled cell medians, spread and eligibility."""
    keys = [*STRATUM, 'method']
    grouped = frame.groupby(keys)[statistics]
    summary = grouped.median().add_suffix('_median')
    for quantile in (0.25, 0.75):
        summary = summary.join(grouped.quantile(quantile).add_suffix(f'_q{int(quantile * 100)}'))
    summary = summary.join(grouped.size().rename('n_cells'))
    summary = summary.join(frame.groupby(keys).prereg_eligible.sum().rename('n_eligible'))
    return summary.reset_index()


def _strata(frame):
    return sorted(set(zip(frame.kind, frame.metric, frame.endpoint)))


def _panels(plt, rows, width, unit=0.36, floor=2.4):
    """One stacked panel per stratum, each sized to the number of methods it lists."""
    ratios = [max(floor, unit * n + 1.4) for n in rows]
    fig, axes = plt.subplots(len(rows), 1, figsize=(width, sum(ratios)), squeeze=False,
                             gridspec_kw=dict(height_ratios=ratios), constrained_layout=True)
    return fig, axes.ravel()


def _ceiling(values, floor):
    """Upper axis limit that keeps the data legible without clipping."""
    finite = [v for v in values if v == v and v is not None]
    return max(floor, max(finite, default=floor) * 1.18)


def figure_calibration(ranking, save, colour, alpha):
    """Design-truth FDR against recall: the error-control versus sensitivity trade."""
    if ranking.empty:
        return None
    strata = _strata(ranking)
    summary = _cell_summary(ranking, ['fdr', 'source_recall'])
    counts = [int(((summary.kind == k) & (summary.metric == m) & (summary.endpoint == e)).sum())
              for k, m, e in strata]
    fig, axes = _panels(save.plt, counts, 6.4, unit=0.0, floor=2.6)
    for ax, (kind, metric, endpoint) in zip(axes, strata):
        cell = summary[(summary.kind == kind) & (summary.metric == metric) & (summary.endpoint == endpoint)]
        for _, row in cell.iterrows():
            _mark(ax, row.source_recall_median, row.fdr_median, colour[row.method],
                  row.n_eligible, row.n_cells)
        ax.axhline(alpha, color='black', linestyle='--', linewidth=1)
        ax.set_xlim(-0.02, _ceiling(cell.source_recall_median, 0.05))
        ax.set_ylim(-0.005, _ceiling(list(cell.fdr_median) + [alpha], 0.06))
        ax.set_xlabel('median recall of the design truth')
        ax.set_ylabel(f'median FDR (dashed: alpha = {alpha:g})')
        ax.set_title(f'{kind} / {metric} / {endpoint}   ({len(cell)} methods)', fontsize=9)
    return save(fig, 'f1_calibration_vs_sensitivity')


def figure_null_calibration(ranking, save, colour):
    """Null p-value rejection rate against the two preregistered null criteria."""
    if ranking.empty:
        return None
    strata = _strata(ranking)
    summary = _cell_summary(ranking, ['null_fpr'])
    panels = [summary[(summary.kind == k) & (summary.metric == m) & (summary.endpoint == e)]
              .sort_values('null_fpr_median') for k, m, e in strata]
    fig, axes = _panels(save.plt, [len(p) for p in panels], 6.8, unit=0.34)
    for ax, (kind, metric, endpoint), cell in zip(axes, strata, panels):
        positions = range(len(cell))
        for position, (_, row) in zip(positions, cell.iterrows()):
            shade = colour[row.method] if row.n_eligible >= row.n_cells else '0.80'
            hatch = None if row.n_eligible >= row.n_cells else '//'
            ax.barh(position, row.null_fpr_q75, color=shade, hatch=hatch, alpha=0.35, height=0.6)
            ax.barh(position, row.null_fpr_median, color=shade, hatch=hatch, height=0.6)
            ax.plot([row.null_fpr_q25, row.null_fpr_q75], [position, position],
                    color='black', linewidth=1, zorder=4)
            ax.text(row.null_fpr_q75, position, f'  {int(row.n_eligible)}/{int(row.n_cells)} cells',
                    va='center', fontsize=7, color='0.25')
        ax.axvline(NULL_CONSERVATIVE, color='0.45', linestyle=':', linewidth=1)
        ax.axvline(NULL_LIMIT, color='black', linestyle='--', linewidth=1)
        ax.set_yticks(list(positions), list(cell.method), fontsize=8)
        ax.set_xlabel('median null rejection rate; whisker is the cell IQR')
        ax.set_title(f'{kind} / {metric} / {endpoint}   (dotted {NULL_CONSERVATIVE:g} valid, '
                     f'dashed {NULL_LIMIT:g} eligibility limit)', fontsize=9)
        ax.set_xlim(0, _ceiling(list(cell.null_fpr_q75) + [NULL_LIMIT], NULL_LIMIT))
    return save(fig, 'f2_null_calibration')


def figure_spillover(ranking, save, colour):
    """How much of a method's design-truth error is a real closure shift instead."""
    if ranking.empty:
        return None
    strata = _strata(ranking)
    summary = _cell_summary(ranking, ['fdr', 'spillover_share'])
    fig, axes = _panels(save.plt, [1] * len(strata), 6.4, unit=0.0, floor=2.6)
    for ax, (kind, metric, endpoint) in zip(axes, strata):
        cell = summary[(summary.kind == kind) & (summary.metric == metric) & (summary.endpoint == endpoint)]
        for _, row in cell.iterrows():
            _mark(ax, row.fdr_median, row.spillover_share_median, colour[row.method],
                  row.n_eligible, row.n_cells)
        upper = _ceiling(cell.spillover_share_median, 0.02)
        reach = min(_ceiling(cell.fdr_median, 0.05), upper)
        ax.plot([0, reach], [0, reach], color='black', linestyle=':', linewidth=1)
        ax.set_xlim(-0.005, _ceiling(cell.fdr_median, 0.05))
        ax.set_ylim(-0.01, upper)
        ax.set_xlabel('median FDR against the design truth')
        ax.set_ylabel('median spillover share')
        ax.set_title(f'{kind} / {metric} / {endpoint}   '
                     '(dotted: every discovery is a closure shift)', fontsize=9)
    return save(fig, 'f5_closure_spillover')


def figure_variance(variance, save):
    """Sequential variance shares per factor, one bar per endpoint, metric and statistic."""
    if variance.empty:
        return None
    factors = [c for c in ['method', 'generator', 'depth', 'fraction_da', 'direction', 'residual']
               if c in variance.columns]
    order = variance.sort_values(STRATUM + ['statistic'], ascending=False).reset_index(drop=True)
    fig, ax = save.plt.subplots(figsize=(9.5, 0.34 * len(order) + 2.4), constrained_layout=True)
    positions = range(len(order))
    bottom = [0.0] * len(order)
    for factor in factors:
        values = pd.to_numeric(order[factor], errors='coerce').fillna(0.0).to_list()
        ax.barh(list(positions), values, left=bottom, height=0.72,
                label=factor, alpha=0.85, edgecolor='white', linewidth=0.4)
        bottom = [b + v for b, v in zip(bottom, values)]
    ax.set_yticks(list(positions),
                  [f"{r.endpoint} / {r.metric} / {r.statistic}" for _, r in order.iterrows()],
                  fontsize=7)
    ax.set_xlabel('share of between-cell variance, factors in the recorded order')
    ax.legend(ncol=len(factors), fontsize=8, frameon=False, loc='lower center',
              bbox_to_anchor=(0.5, 1.005))
    return save(fig, 'f3_variance_shares')


def figure_rank_uncertainty(ranks, save, colour):
    """Rank stability across conditions; the bootstrap interval is drawn behind each box."""
    if ranks.empty:
        return None
    frame = ranks.dropna(subset=['median_rank'])
    if frame.empty:
        return None
    strata = _strata(frame)
    panels = [sorted(set(frame.loc[(frame.kind == k) & (frame.metric == m) & (frame.endpoint == e), 'method']),
                     key=lambda name: float(frame.loc[frame.method.eq(name), 'median_rank'].median()))
              for k, m, e in strata]
    pairs = [(s, p) for s, p in zip(strata, panels) if len(p) > 1]
    if not pairs:
        return None
    fig, axes = _panels(save.plt, [len(p) for _, p in pairs], 6.8, unit=0.34)
    for ax, ((kind, metric, endpoint), methods) in zip(axes, pairs):
        cell = frame[(frame.kind == kind) & (frame.metric == metric) & (frame.endpoint == endpoint)]
        data = [cell.loc[cell.method.eq(name), 'median_rank'].to_numpy(dtype=float) for name in methods]
        positions = list(range(len(methods)))[::-1]
        ax.boxplot(data, positions=positions, widths=0.55, vert=True, showfliers=False,
                   medianprops=dict(color='black'), boxprops=dict(color='0.6'))
        for position, name in zip(positions, methods):
            rows = cell.loc[cell.method.eq(name)]
            low = float(pd.to_numeric(rows['low'], errors='coerce').median())
            high = float(pd.to_numeric(rows['high'], errors='coerce').median())
            if low == low and high == high:
                ax.plot([low, high], [position, position], color=colour[name],
                        linewidth=2.6, alpha=0.55, zorder=3)
        ax.set_yticks(positions, methods, fontsize=8)
        ax.set_ylim(-0.7, len(methods) - 0.3)
        ax.set_xlim(0.4, max(ax.get_xlim()[1], 2.6))
        ax.set_xlabel('median rank of recall within a condition, 1 is best')
        ax.set_title(f'{kind} / {metric} / {endpoint}   '
                     '(line: median 95% bootstrap interval)', fontsize=9)
    return save(fig, 'f4_rank_uncertainty')


def write_figures(ranking, variance, ranks, figures, alpha=0.05, formats=('png', 'svg'),
                  mpl_config=None):
    """Render the analysis figures into `figures`; empty list if matplotlib is absent."""
    figures = Path(figures)
    figures.mkdir(parents=True, exist_ok=True)
    try:
        plt = _pyplot(mpl_config or figures / '.mpl')
    except ImportError:
        print('matplotlib is not installed; skipping analysis figures', flush=True)
        return []
    save = _save_factory(plt, figures, formats)
    colour = _palette([f for frame in (ranking, ranks) if 'method' in frame.columns for f in frame.method])
    written = []
    for stem in (figure_calibration(ranking, save, colour, alpha),
                 figure_null_calibration(ranking, save, colour),
                 figure_spillover(ranking, save, colour),
                 figure_variance(variance, save),
                 figure_rank_uncertainty(ranks, save, colour)):
        if stem:
            written.append(stem)
    return sorted(written)


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
    frames, guide, rankings = [], [], []
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
            ranking['tier'] = tier
            rankings.append(ranking)
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
    variance = pd.DataFrame(variance_rows)
    save_table(variance, output / 'variance.tsv')
    save_table(pd.concat(contrasts, ignore_index=True) if contrasts else pd.DataFrame(), output / 'contrasts.tsv')
    rank_frame = pd.DataFrame(ranks)
    save_table(rank_frame, output / 'rank_uncertainty.tsv')
    write_json(output / 'consensus.json', dict(k=CONSENSUS_K, status='specified_only_not_evaluated'))
    figures = write_figures(pd.concat(rankings, ignore_index=True), variance, rank_frame,
                            output / 'figures')
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
        'Batch confounding and consensus performance have not been evaluated.\n\n'
        '## Figures\n\n'
        'Each panel pools every cell of one endpoint and input metric; a filled marker means the\n'
        'method is eligible in all of them, a hollow marker that eligibility is partial and grey\n'
        'that it is absent. Pooling over design conditions is descriptive.\n\n'
        + ''.join(f'- `figures/{name}`\n' for name in figures)
        + ('- no figures: matplotlib is not installed in this environment\n' if not figures else ''))
    print(f'Wrote conditional analysis to {output} ({len(figures)} figures)', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
