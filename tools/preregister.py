#!/usr/bin/env python3
"""Freeze thresholds/analysis code and tag the preregistration revision.

Rule (plan.md): nothing in Tier 3-4 tunes methods or thresholds. Thresholds
and analysis code are frozen and committed before Tier 2 starts.

  python tools/preregister.py --tag preregistration --methods config/methods.yaml

Records: active method list + parameters, metrics/eligibility rules, analysis
code hash (src/magician + workflow/scripts), git revision. Creates an
annotated git tag. Any post-hoc analysis after the tag is labelled exploratory.
Predefined consensus rule: discovery if >=3 eligible methods agree (frozen here).
"""
from pathlib import Path
import argparse
import hashlib
import subprocess
import sys
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

CONSENSUS_K = 3


def code_hash():
    h = hashlib.sha256()
    for path in sorted((ROOT / "src/magician").glob("*.py")):
        h.update(path.read_bytes())
    for path in sorted((ROOT / "workflow/scripts/da").glob("*.R")):
        h.update(path.read_bytes())
    h.update((ROOT / "workflow/scripts/run_da.R").read_bytes())
    h.update((ROOT / "workflow/scripts/da/contract.R").read_bytes())
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tag", default="preregistration")
    p.add_argument("--methods", default="config/methods.yaml")
    p.add_argument("--out", default="provenance/preregistration.json")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()
    reg = yaml.safe_load((ROOT / args.methods).read_text())
    active = [m["id"] for m in reg["methods"] if m.get("status") == "active"]
    from magician.io import write_json
    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(ROOT),
                         capture_output=True, text=True).stdout.strip()
    record = {
        "tag": args.tag,
        "git_revision": rev,
        "methods_active": sorted(active),
        "methods_file": args.methods,
        "alpha": 0.05,
        "primary_truth": "implanted design truth, common to all methods; endpoint lanes are secondary",
        "eligibility": {"fdr_ci_high_le": 0.10, "null_fpr_le": 0.07, "null_fpr_conservative_lt": 0.03,
                        "failure_rate_lt": 0.05, "min_null_runs": 100, "min_spiked_runs": 50},
        "consensus_rule": f"discovery if >={CONSENSUS_K} eligible methods agree",
        "analysis_code_sha256": code_hash(),
        "note": "Thresholds/analysis frozen before Tier 2; post-tag changes are exploratory.",
    }
    if args.dry_run:
        print(yaml.safe_dump(record, sort_keys=False))
        return 0
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    write_json(out, record)
    print(f"wrote {out} at {rev} with {len(active)} active methods")
    subprocess.run(["git", "add", args.methods, "src/magician", "workflow/scripts/da",
                    "config/zero_handling.yaml"], cwd=str(ROOT))
    rc = subprocess.run(["git", "commit", "-m", f"Freeze {args.tag}: {len(active)} methods, "
                         f"code {record['analysis_code_sha256'][:12]}"],
                        cwd=str(ROOT)).returncode
    if rc:
        print("commit failed (nothing to commit?); tagging current HEAD")
    subprocess.run(["git", "tag", "-a", args.tag, "-m",
                    f"Preregistration: {len(active)} methods, consensus k={CONSENSUS_K}"],
                   cwd=str(ROOT))
    print(f"tagged {args.tag}; LOCOM excluded (package failure) is recorded in methods.yaml")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
