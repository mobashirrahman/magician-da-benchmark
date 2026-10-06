"""Compact standalone report, plots and actual resource measurements."""
from pathlib import Path
import html
import json
import os
import pandas as pd
from .io import read_json, write_json, save_table, tree_bytes
from .config import resolve, registry
from .storage import preflight

RESOURCE_TABLES = ("rule_resources.tsv", "method_settings.tsv")


def _read(path):
    return pd.read_csv(path, sep="\t", keep_default_na=False, na_values=["NA"])


def make_report(cfg, summary, benchmarks, output, usage_output):
    root = Path(summary)
    write_json(resolve(cfg["output_dir"]) / "provenance/config.json", cfg)
    ranking = _read(root / "ranking.tsv")
    recovery = _read(root / "recovery.tsv")
    scores = _read(root / "scores.tsv")
    settings = _read(root / "method_settings.tsv") if (root / "method_settings.tsv").exists() \
        else pd.DataFrame(columns=["method_id", "method_family", "variant", "endpoint"])
    measurements = []
    for file in benchmarks:
        df = pd.read_csv(file, sep="\t")
        if len(df):
            row = df.iloc[0].to_dict()
            row["job"] = str(Path(file).relative_to(resolve(cfg["output_dir"])))
            measurements.append(row)
    resources = pd.DataFrame(measurements)
    save_table(resources, root / "rule_resources.tsv")
    # Measured cost per method, from the job benchmark rather than the R timer.
    measured = []
    for _, row in ranking.drop_duplicates(["kind", "metric", "endpoint", "method"]).iterrows():
        pattern = f"da/{row.kind}/{row.metric}/{row.method}.tsv"
        jobs = resources.loc[resources.job.str.endswith(pattern)] if len(resources) else resources
        measured.append(dict(kind=row.kind, metric=row.metric, endpoint=row.endpoint,
                             method=row.method, method_family=row.method_family,
                             n_jobs=len(jobs),
                             total_seconds=float(pd.to_numeric(jobs.get("s"), errors="coerce").sum()) if len(jobs) else None,
                             longest_seconds=float(pd.to_numeric(jobs.get("s"), errors="coerce").max()) if len(jobs) else None,
                             largest_rss_kb=float(pd.to_numeric(jobs.get("max_rss"), errors="coerce").max()) if len(jobs) else None))
    cost = pd.DataFrame(measured)
    save_table(cost, root / "method_resources.tsv")
    unusable = _read(root / "scores.tsv").loc[lambda d: ~d.status.isin(["success"])]
    unusable = unusable[["case", "kind", "metric", "endpoint", "method", "status", "method_message"]].drop_duplicates()
    usages = dict(retained_output_bytes_before_report=tree_bytes(resolve(cfg["output_dir"])),
                  managed_cache_bytes=tree_bytes(resolve(cfg["cache_dir"])),
                  storage_budget_bytes=int(cfg["storage"]["budget_gb"] * 10**9),
                  total_measured_job_seconds=float(pd.to_numeric(resources.get("s", pd.Series(dtype=float)), errors="coerce").sum()),
                  largest_measured_job_rss_mb=float(pd.to_numeric(resources.get("max_rss", pd.Series(dtype=float)), errors="coerce").max()) if len(resources) else None,
                  note="Job RSS is not concurrent total RAM. Size snapshot excludes this report and its own unfinished benchmark.")
    write_json(usage_output, usages)
    # Standalone scientific plots, suitable for export; PNG and SVG are small.
    mpl_config = resolve(cfg["cache_dir"]) / "matplotlib"
    mpl_config.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(mpl_config)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    spike = scores.loc[scores.scenario.eq("spiked") & scores.status.eq("success")]
    strata = sorted(set(zip(scores.kind, scores.metric, scores.endpoint)))
    if strata:
        fig, axes = plt.subplots(len(strata), 2, figsize=(13, 3.2 * len(strata)), squeeze=False, constrained_layout=True)
        for row, (kind, metric, endpoint) in zip(axes, strata):
            subset = spike.loc[(spike.kind == kind) & (spike.metric == metric) & (spike.endpoint == endpoint)]
            methods = list(dict.fromkeys(subset.method))
            for ax, statistic, label in zip(row, ["fdr", "source_recall"],
                                            ["False discovery proportion against the design truth", "Recall of the design truth"]):
                if len(subset):
                    means = subset.groupby("method")[statistic].mean().reindex(methods)
                    ax.scatter(range(len(methods)), means.values, marker="o")
                    ax.set_xticks(range(len(methods)), methods, rotation=60, ha="right")
                ax.set_title(f"{kind} / {metric} / {endpoint}")
                ax.set_ylabel(label)
                ax.set_ylim(-0.02, 1.02)
                if statistic == "fdr":
                    ax.axhline(cfg["analysis"]["alpha"], color="grey", linestyle="--", linewidth=1)
        fig.savefig(root / "performance.svg")
        fig.savefig(root / "performance.png", dpi=160)
        plt.close(fig)
    recs = read_json(root / "recommendations.json")
    from .evaluation import CONDITION_COLUMNS
    condition_text = lambda r: "; ".join(f"{c}={r[c]}" for c in CONDITION_COLUMNS if c in r)
    conclusions = "".join(
        f"<li>{html.escape(r['kind'])}/{html.escape(r['metric'])}/"
        f"{html.escape(r['endpoint'])}: "
        f"{html.escape(r['best_method'] or ', '.join(r.get('candidate_methods', [])) or 'no eligible method')} "
        f"— {html.escape(r['conclusion'])} ({r['evidence']}) "
        f"{html.escape(condition_text(r))}</li>" for r in recs)
    corrected = "generator_version" in ranking.columns
    eligibility = (
        "Corrected Tier 1 candidates are assessed separately by condition. The upper FDR confidence bound "
        "must be ≤0.10, raw null p-value rejection rate must be ≤0.07 (below 0.03 is flagged conservative), "
        "and failure rate must be &lt;5%. Recall must be positive and significant unevaluable features absent. "
        "Recommendations also require at least 100 null and 50 spiked runs per condition; null runs are shared "
        "by spiked conditions that differ only in effect factors. Thin cells remain exploratory."
        if corrected else
        "Eligibility requires successful evaluable runs, positive source recall, mean spiked FDR ≤α, "
        "null probability of any false discovery ≤α, and no significant unassigned MAGs.")
    lane_note = "".join(
        f"<li><b>{html.escape(name)}</b>: truth lane <code>{html.escape(lane['lane'])}</code> — "
        f"{html.escape(lane['description'])}</li>"
        for name, lane in sorted(registry(cfg).endpoints.items()))
    settings_table = settings.to_html(index=False, escape=True) if len(settings) else "<p>none recorded</p>"
    body = f"""<!doctype html><html lang="en"><meta charset="utf-8"><title>MAGICIAN benchmark</title>
    <style>body{{font:15px system-ui;max-width:1500px;margin:30px auto;padding:0 20px}}table{{border-collapse:collapse;width:100%;font-size:12px}}th,td{{padding:6px;border:1px solid #ddd}}img{{max-width:100%}}section{{overflow-x:auto}}</style>
    <h1>MAGICIAN differential abundance benchmark</h1>
    <p>Contrast: Treatment / Control. Source genomes and MAGs are quantified separately. Discoveries are
    scored against the design truth, the features the design changed, which is the same set for every
    method. The truth lane for each method's declared endpoint is scored beside it: <code>endpoint_fdr</code>
    is unavailable where that scale has no null feature, and <code>spillover_share</code> is the share of
    discoveries that are real shifts on the method's own scale but lie outside the design truth.</p>
    <p>Rankings are specific to these simulations and recorded conditions. Small validation runs establish
    pipeline execution, not statistical reliability. Tier 1 does not test MAG reconstruction or real-data validity.</p>
    <h2>Endpoints and truth lanes</h2><ul>{lane_note}</ul>
    <h2>Candidates</h2><ul>{conclusions}</ul>
    <p>{eligibility} Ordering uses source recall,
    then average precision and FDR. Corrected Tier 1 uses raw p-value average precision;
    unavailable raw rankings are retained as unavailable. Variants of one method family are reported as one
    candidate. Bootstrap intervals resample whole independent seeds; the null rate also carries an exact
    binomial interval.</p>
    <h2>Performance by endpoint and input</h2><img src="performance.svg" alt="FDR and source recall by method">
    <p>The overview plots average across design conditions. Use the condition-specific table and candidate
    records for interpretation; pooled plot averages are descriptive.</p>
    <section>{ranking.to_html(index=False, escape=True)}</section>
    <h2>Declared method settings</h2><section>{settings_table}</section>
    <h2>Measured resources by method</h2><section>{cost.to_html(index=False, escape=True)}</section>
    <h2>Genome recovery</h2><section>{recovery.to_html(index=False, escape=True)}</section>
    <p>Feature FDR/recall are conditional on uniquely assigned MAGs. Source recall includes unrecovered DA
    genomes as missed changes. Multiple bins for a source are reported, and source discoveries are
    deduplicated. Ambiguous and unassigned features are excluded from truth scoring, with their significant
    calls reported explicitly. Fold-change error is only computed when a method's reported units match its
    truth lane's coordinates.</p>
    <h2>Runs that produced no usable result</h2><section>{unusable.to_html(index=False, escape=True)}</section>
    <h2>Resources</h2><pre>{html.escape(json.dumps(usages, indent=2, sort_keys=True))}</pre>
    <p>See scores.tsv, ranking.tsv, recovery.tsv, method_settings.tsv, method_resources.tsv,
    rule_resources.tsv, recommendations.json and per-experiment feature_evaluation.tsv for machine-readable
    results.</p></html>"""
    Path(output).write_text(body)
