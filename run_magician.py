#!/usr/bin/env python3
"""Convenient v2 entry point; the standard Snakemake CLI also works."""
from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
from magician.config import load, resolve


def main():
    p = argparse.ArgumentParser(description="Run the MAGICIAN differential abundance benchmark")
    p.add_argument("--configfile", default="config/config.yaml")
    p.add_argument("--cores", type=int, default=4)
    p.add_argument("--dry-run", "-n", action="store_true")
    p.add_argument("--no-conda", action="store_true", help="Use tools on PATH (development only)")
    p.add_argument("--snakemake", default="snakemake")
    args, extra = p.parse_known_args()
    cfg = load(args.configfile)
    if args.cores < 1:
        p.error("cores must be positive")
    command = [args.snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"),
               "--directory", str(ROOT), "--workflow-profile", str(ROOT / ("profiles/testing" if args.no_conda else "profiles/default")),
               "--configfile", str(resolve(args.configfile)), "--cores", str(args.cores),
               "--conda-prefix", str(resolve(cfg["cache_dir"]) / "conda")]
    if args.dry_run:
        command += ["--dry-run"]
    command += extra
    env = os.environ.copy()
    env["CONDA_PKGS_DIRS"] = str(resolve(cfg["cache_dir"]) / "packages")
    env["XDG_CACHE_HOME"] = str(resolve(cfg["cache_dir"]) / "xdg")
    return subprocess.call(command, env=env)


if __name__ == "__main__":
    raise SystemExit(main())
