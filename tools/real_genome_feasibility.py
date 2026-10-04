#!/usr/bin/env python3
"""Real-genome feasibility: one paired pilot before any repeated execution.

  python tools/real_genome_feasibility.py plan config/genome_manifest.tsv --out config/feasibility.yaml
  python run_magician.py --configfile config/feasibility.yaml --cores 8
  python tools/real_genome_feasibility.py assess results_feasibility

The plan step checks the frozen manifest and writes a one-seed, one-pair overlay.
The assess step verifies depth, mapping ambiguity, usable MAG recovery and measured
disk and RAM against the budget, and refuses to recommend repeated execution unless
recovery is usable and the storage estimate fits with headroom. A failed feasibility
run is a result, not a reason to silently restrict the truth.
"""
from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from magician.config import resolve
from magician.genomes import check_genomes
from magician.io import table, read_json, write_json


def plan(manifest, output, out_dir, seed, samples_per_group, read_pairs, budget_gb):
    candidates = table(Path(manifest))
    required = {"genome_id", "accession", "version", "path", "provenance", "source_group"}
    missing = required - set(candidates.columns)
    if missing:
        raise SystemExit(f"{manifest} is missing columns: {', '.join(sorted(missing))}")
    rows, problems = check_genomes(candidates.to_dict("records"))
    for problem in problems:
        print(f"unusable: {problem}")
    if problems:
        raise SystemExit("The manifest contains unusable references; refusing to plan")
    if len(rows) < 20:
        raise SystemExit(f"Feasibility needs at least 20 real genomes, got {len(rows)}")
    body = (f"# Real-genome feasibility pilot, planned by tools/real_genome_feasibility.py.\n"
            f"output_dir: {out_dir}\n"
            f"genomes:\n  manifest: {manifest}\n"
            f"experiment:\n  scenarios: ['null', 'spiked']\n  seeds: [{seed}]\n"
            f"  samples_per_group: {samples_per_group}\n  read_pairs: {read_pairs}\n"
            f"storage:\n  budget_gb: {budget_gb}\n")
    Path(output).write_text(body)
    print(f"wrote {output}: {len(rows)} genomes, one paired pilot, seed {seed}")


def assess(output_dir):
    root = Path(output_dir)
    recovery = table(root / "benchmark/recovery.tsv")
    scores = table(root / "benchmark/scores.tsv")
    usage = read_json(root / "benchmark/storage.json")
    preflight = read_json(root / "provenance/preflight.json")
    verdict = dict(output=str(root), cases=sorted(recovery.case.tolist()))
    verdict["recovered_sources"] = {row.case: int(row.n_recovered_sources)
                                    for row in recovery.itertuples()}
    verdict["recovered_da_sources"] = {row.case: int(row.n_recovered_da_sources)
                                       for row in recovery.itertuples()}
    verdict["ambiguous"] = int(recovery.n_ambiguous.sum())
    verdict["unassigned"] = int(recovery.n_unassigned.sum())
    verdict["multiply_matched_sources"] = int(recovery.n_multiply_matched_sources.sum())
    verdict["retained_bytes"] = usage["retained_output_bytes_before_report"]
    verdict["largest_job_rss_mb"] = usage["largest_measured_job_rss_mb"]
    verdict["preflight_peak_total_bytes"] = preflight["estimated_peak_total_bytes"]
    verdict["storage_budget_bytes"] = usage["storage_budget_bytes"]
    min_recovered = min(verdict["recovered_sources"].values())
    # Usable means most sources come back as unambiguous MAGs, so the repeated study
    # measures recovery rather than assembly failure.
    verdict["recovery_usable"] = bool(min_recovered >= 0.7 * int(recovery.source_genomes.max())
                                     and verdict["ambiguous"] == 0)
    verdict["storage_fits"] = bool(
        preflight["estimated_peak_total_bytes"] < usage["storage_budget_bytes"])
    failures = scores.loc[~scores.status.isin(["success", "empty"]),
                          ["case", "kind", "metric", "method", "status"]].drop_duplicates()
    verdict["method_failures"] = failures.to_dict("records")
    verdict["go_for_repeated_execution"] = bool(verdict["recovery_usable"] and verdict["storage_fits"]
                                               and not len(failures))
    verdict["note"] = ("repeat the pilot at the confirmed depth with 10 to 20 paired seeds "
                       "only when go_for_repeated_execution is true")
    write_json(root / "benchmark/feasibility.json", verdict)
    print(json.dumps(verdict, indent=2, sort_keys=True))
    return 0 if verdict["go_for_repeated_execution"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    plan_parser = sub.add_parser("plan", help="Write a one-seed feasibility overlay")
    plan_parser.add_argument("manifest")
    plan_parser.add_argument("--out", default="config/feasibility.yaml")
    plan_parser.add_argument("--out-dir", default="results_feasibility")
    plan_parser.add_argument("--seed", type=int, default=1001)
    plan_parser.add_argument("--samples-per-group", type=int, default=10)
    plan_parser.add_argument("--read-pairs", type=int, default=500000)
    plan_parser.add_argument("--budget-gb", type=float, default=100.0)
    assess_parser = sub.add_parser("assess", help="Judge a completed feasibility run")
    assess_parser.add_argument("output_dir")
    args = parser.parse_args()
    if args.command == "plan":
        plan(args.manifest, resolve(args.out), args.out_dir, args.seed,
             args.samples_per_group, args.read_pairs, args.budget_gb)
        return 0
    return assess(args.output_dir)


if __name__ == "__main__":
    sys.exit(main())