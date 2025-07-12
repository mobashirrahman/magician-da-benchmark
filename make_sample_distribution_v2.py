#!/usr/bin/env python3
"""
Generate CAMISIM abundance tables from a YAML spec. (v2) (CAMISIM_complex)
Usage:  python make_sample_distribution_v2.py config.yaml
"""
import pathlib, csv, math, random, sys
from typing import Dict, List
import numpy as np
import yaml, typer

app = typer.Typer()

# ──────────────────────────────────────────────────────────────
def dirichlet_sample(baseline: np.ndarray, conc: float, rng):
    alpha = baseline / baseline.sum() * conc
    return rng.dirichlet(alpha)

def lognormal_vec(mu, sigma, size, rng):
    return rng.lognormal(mean=mu, sigma=sigma, size=size)

# ──────────────────────────────────────────────────────────────
@app.command()
def main(config_path: pathlib.Path):
    # 0 ▸ read config
    cfg = yaml.safe_load(config_path.read_text())
    rng = np.random.default_rng(cfg.get("seed", None))
    random.seed(cfg.get("seed", None))

    ref_dir   = pathlib.Path(cfg["reference_dir"])
    out_dir   = pathlib.Path(cfg["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    genomes   = sorted(ref_dir.glob("*.*"))
    n_genomes = len(genomes)
    if n_genomes == 0:
        typer.echo("❌ No genomes found in reference_dir"); sys.exit(1)

    # 1 ▸ baseline weights
    pri = cfg["abundance_prior"]
    if pri["distribution"] == "lognormal":
        baseline = lognormal_vec(pri["mu"], pri["sigma"], n_genomes, rng).round().astype(int)
    else:
        raise NotImplementedError("Only lognormal prior implemented")

    baseline[baseline < 1] = 1

    # 2 ▸ choose DA genomes & fold-changes
    da_cfg = cfg["differential_abundance"]
    da_idx = rng.choice(n_genomes, size=da_cfg["n_da"], replace=False)
    folds  = rng.choice(da_cfg["effect_sizes"], size=da_cfg["n_da"], replace=True)
    fc_vec = np.ones(n_genomes)
    fc_vec[da_idx] = folds

    # 3 ▸ build sample table
    samples, groups = [], []
    for g, n in cfg["groups"].items():
        for i in range(1, n + 1):
            sid = f"{g[:1].upper()}{str(i).zfill(2)}"
            samples.append(sid)
            groups.append(g)
    n_samples = len(samples)

    # 4 ▸ library depths
    depth_cfg = cfg["library_depth"]
    lib_sizes = lognormal_vec(depth_cfg["mu"], depth_cfg["sigma"], n_samples, rng).round().astype(int)

    # 5 ▸ per-sample abundances (Dirichlet–multinomial)
    conc_param = 1.0           # smaller → more dispersion
    weight_mat = np.zeros((n_samples, n_genomes), dtype=int)
    for s in range(n_samples):
        # multiply by fold-change if sample belongs to a treatment group
        is_treat = (groups[s] != "Control")
        baseline_this = baseline * (fc_vec if is_treat else 1)
        rel_abund = dirichlet_sample(baseline_this, conc_param, rng)
        weight_mat[s] = rng.multinomial(lib_sizes[s], rel_abund)

    # 6 ▸ write read_depth
    with (out_dir / "read_depth.txt").open("w") as fh:
        for sid, d in zip(samples, lib_sizes):
            print(sid, d, sep="\t", file=fh)

    # 7 ▸ sample_distributions.tsv
    header = ["genomes", "seq_type", *samples]
    with (out_dir / "sample_distributions.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for j, g in enumerate(genomes):
            w.writerow([g.resolve(), "chromosome", *weight_mat[:, j]])

    # 8 ▸ truth table
    with (out_dir / "truth_table.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["genome", "is_DA", "true_log2FC",
                    "baseline_weight"])
        for j, g in enumerate(genomes):
            is_da = j in da_idx
            log2fc = math.log2(fc_vec[j])
            w.writerow([g.name, is_da, log2fc, baseline[j]])

    typer.echo(f" Design written to {out_dir.resolve()}")

# ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    app()
