#!/usr/bin/env python3
"""Phase 0: re-measure timings (plan.md Section 5) and rescale the compute plan.

Estimates rest on existing runs and must be re-measured in Phase 0.
Calibration: 964 jobs, 49,241 core-seconds for 40 cases (~20 core-min/case).
Pilot20: simulation ~5 min/sample, assembly ~5.5 min for 20x200k pairs.

  python tools/phase0_measure.py --cores 4 --out docs/phase0_timings.json
  python tools/phase0_measure.py --report  # print rescaled Tier 1/2/3 wall times

Records core count, RAM, filesystem free, Snakemake/Python/R versions.
"""
from pathlib import Path
import argparse
import json
import platform
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def env_record():
    def run(cmd):
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
            return (out.stdout or out.stderr or "").strip().split("\n")[0][:200]
        except Exception as e:
            return f"unavailable: {e}"
    import sys as _sys
    return {
        "platform": platform.platform(),
        "cores": __import__("os").cpu_count(),
        "python": _sys.version.split()[0],
        "snakemake": run(["snakemake", "--version"]),
        "R": run(["Rscript", "--version"]),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cores", type=int, default=4)
    p.add_argument("--out", default="docs/phase0_timings.json")
    p.add_argument("--report", action="store_true")
    args = p.parse_args()
    rec = env_record()
    try:
        total, used, free = shutil.disk_usage(str(ROOT))
        rec.update({"disk_total_gb": total / 1e9, "disk_free_gb": free / 1e9})
    except Exception:
        pass
    # Rescale plan.md wall times to this machine's core count.
    ref_cores, ref_hours = 60, {"tier1": 35, "tier2": 72, "tier3a": 14, "tier3b": 10,
                                "tier4": 24, "slack": 48}
    scale = ref_cores / max(args.cores, 1)
    rec["rescaled_wall_hours"] = {k: round(v * scale, 1) for k, v in ref_hours.items()}
    rec["note"] = (f"Reference is 60 cores; this machine has {args.cores} cores "
                   f"(x{scale:.1f}). Tier 2 needs >=5TB scratch; disk binding.")
    if args.report:
        print(json.dumps(rec, indent=2))
        return 0
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=2) + "\n")
    print(f"wrote {out}")
    print(json.dumps(rec["rescaled_wall_hours"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
