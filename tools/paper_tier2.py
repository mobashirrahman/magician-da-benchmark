#!/usr/bin/env python3
"""Tier 2 runner: 60 full-pipeline experiments (oracle/MAG/reference-70%).

Cells (plan.md): community 50/150/300 genomes x strain
(species-unique vs 20% species with 2-3 strains) x depth 2M/10M read pairs
(40 samples/experiment) x null/spiked (20% DA balanced) x 5 seeds.
That is 3x2x2x2x5=120; the 300-genome strain cells are dropped to fit budget,
leaving 60 experiments. Each cell is one co-assembly (MetaBAT2 default),
one oracle table (source genomes), one MAG table and one incomplete-reference
table (70% present).

Ablations on one design cell x5 seeds: second binner (SemiBin2/COMEBin),
completeness filter (none vs >=50/<10 vs >=90/<5), per-sample vs co-assembly.

  python tools/paper_tier2.py --plan
  python tools/paper_tier2.py --generate-only --manifest config/genome_manifest.tsv
  python tools/paper_tier2.py --run --cores 60 --manifest config/genome_manifest.tsv

Requires a frozen real manifest (tools/real_genomes.py freeze). At most 3
experiments run at a time (disk: ~100GB FASTQ per 40-sample 10M experiment);
reads are deleted after counting. Always run with --cores 60 on the server,
memory capped below 400GB.
"""
from pathlib import Path
import argparse
import subprocess
import sys
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def cells():
    out = []
    for community in (50, 150, 300):
        for strain in ("unique", "mixed"):
            if community == 300 and strain == "mixed":
                continue  # dropped to fit budget (plan.md)
            for depth in (2_000_000, 10_000_000):
                out.append(dict(community=community, strain=strain, depth=depth))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--plan", action="store_true")
    p.add_argument("--generate-only", action="store_true")
    p.add_argument("--run", action="store_true")
    p.add_argument("--manifest", default="config/genome_manifest.tsv")
    p.add_argument("--base", default="config/paper_tier2.yaml")
    p.add_argument("--outdir", default="cache/paper_tier2")
    p.add_argument("--cores", type=int, default=60)
    p.add_argument("--seeds", nargs="*", type=int, default=[101, 102, 103, 104, 105])
    args = p.parse_args()
    table = cells()
    if args.plan:
        print(f"{len(table)} cells x 5 seeds x 2 scenarios = {len(table) * 10} experiments")
        for i, c in enumerate(table):
            print(f"cell{i:02d}: community={c['community']} strain={c['strain']} depth={c['depth']}")
        print("Reduced-scope 60-expt option: --seeds 101 102 103 (10 cells x3 seeds x2 = 60).")
        print("Ablations (one cell x5 seeds): binner semibin2/comebin, "
              "completeness medium/strict, assembly single (if time allows).")
        return 0
    base = yaml.safe_load((ROOT / args.base).read_text())
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    configs = []
    for i, c in enumerate(table):
        cfg = yaml.safe_load((ROOT / args.base).read_text())
        cfg["genomes"]["manifest"] = args.manifest
        cfg["experiment"]["read_pairs"] = c["depth"]
        cfg["experiment"]["seeds"] = args.seeds
        cfg["output_dir"] = f"results_paper_tier2_cell{i:02d}"
        # Record the cell design for variance decomposition (method x factor).
        cfg.setdefault("tier2_cell", {}).update(c)
        cfg["tier2_cell"]["cell_id"] = f"cell{i:02d}"
        cfg["tier2_cell"]["strain_note"] = (
            "species-unique" if c["strain"] == "unique"
            else "20% species with 2-3 related strains (manifest source_group)")
        path = outdir / f"cell{i:02d}.yaml"
        path.write_text(yaml.safe_dump(cfg, sort_keys=False))
        configs.append(path)
    print(f"wrote {len(configs)} cell configs to {outdir} (manifest={args.manifest})")
    print("Strain axis uses manifest source_group: freeze with "
          "--group diverse=N related=M (see tools/real_genomes.py).")
    if args.generate_only or not args.run:
        return 0
    for path in configs:
        # At most 3 experiments in flight (disk binding); run sequentially here,
        # the scheduler parallelises within each experiment (assembly 16 threads).
        cmd = [sys.executable, str(ROOT / "run_magician.py"),
               "--configfile", str(path), "--cores", str(args.cores)]
        print("+", " ".join(cmd))
        rc = subprocess.call(cmd, cwd=str(ROOT))
        if rc:
            print(f"FAILED {path}, stopping tier.")
            return rc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
