#!/usr/bin/env python3
"""
Generate CAMISIM abundance tables from a YAML spec. (v2) (CAMISIM_complex)
Usage:  python make_sample_distribution_v2.py config.yaml

Features:
- Configurable sample groups and replicates (e.g., C01_rep01, C01_rep02)
- Dirichlet-multinomial abundance modeling
- Automatic .gz genome decompression
- Differential abundance simulation
- Realistic noise sources:
  * Sample processing variation
  * Batch effects
  * PCR amplification bias
  * Contamination
  * Zero-inflation (dropouts)
"""
import pathlib, csv, math, random, sys, gzip, shutil
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

    # 1 ▸ decompress any *.gz genome files first so CAMISIM can read them
    for gz_path in ref_dir.glob("*.gz"):
        dest_path = gz_path.with_suffix("")  # strip the .gz extension
        if not dest_path.exists():
            typer.echo(f"⚙️  Decompressing {gz_path.name} → {dest_path.name}")
            with gzip.open(gz_path, "rb") as f_in, open(dest_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        # Remove the original compressed file to avoid duplicate inputs
        try:
            gz_path.unlink()
        except Exception as e:
            typer.echo(f"⚠️  Could not remove {gz_path.name}: {e}")

    genomes   = sorted(p for p in ref_dir.iterdir() if p.is_file() and not p.suffix.endswith('gz'))
    n_genomes = len(genomes)
    if n_genomes == 0:
        typer.echo("❌ No genomes found in reference_dir"); sys.exit(1)

    # 2 ▸ baseline weights
    pri = cfg["abundance_prior"]
    if pri["distribution"] == "lognormal":
        baseline = lognormal_vec(pri["mu"], pri["sigma"], n_genomes, rng).round().astype(int)
    else:
        raise NotImplementedError("Only lognormal prior implemented")

    baseline[baseline < 1] = 1

    # 3 ▸ choose DA genomes & fold-changes
    da_cfg = cfg["differential_abundance"]
    da_idx = rng.choice(n_genomes, size=da_cfg["n_da"], replace=False)
    folds  = rng.choice(da_cfg["effect_sizes"], size=da_cfg["n_da"], replace=True)
    fc_vec = np.ones(n_genomes)
    fc_vec[da_idx] = folds

    # 4 ▸ build sample table
    n_replicates = cfg.get("replicates", 1)  # default to 1 if not specified
    samples, groups = [], []
    for g, n in cfg["groups"].items():
        for i in range(1, n + 1):
            base_sid = f"{g[:1].upper()}{str(i).zfill(2)}"
            if n_replicates > 1:
                for r in range(1, n_replicates + 1):
                    sid = f"{base_sid}_rep{str(r).zfill(2)}"
                    samples.append(sid)
                    groups.append(g)
            else:
                samples.append(base_sid)
                groups.append(g)
    n_samples = len(samples)

    # 5 ▸ library depths
    depth_cfg = cfg["library_depth"]
    lib_sizes = lognormal_vec(depth_cfg["mu"], depth_cfg["sigma"], n_samples, rng).round().astype(int)

    # 6 ▸ per-sample abundances (Dirichlet–multinomial)
    conc_param = 1.0           # smaller → more dispersion
    weight_mat = np.zeros((n_samples, n_genomes), dtype=int)
    
    # Apply noise sources if configured
    noise_cfg = cfg.get("noise_sources", {})
    
    # PCR amplification bias (per-genome)
    pcr_bias = np.ones(n_genomes)
    if noise_cfg.get("pcr_amplification", {}).get("enabled", False):
        pcr_cv = noise_cfg["pcr_amplification"]["cv"]
        pcr_sigma = np.sqrt(np.log(1 + pcr_cv**2))  # convert CV to lognormal sigma
        pcr_bias = rng.lognormal(0, pcr_sigma, size=n_genomes)
        typer.echo(f"📊 Applied PCR amplification bias (CV={pcr_cv:.1%})")
    
    # Batch effects (per-sample)
    sample_batch_effects = np.ones(n_samples)
    if noise_cfg.get("batch_effects", {}).get("enabled", False):
        n_batches = noise_cfg["batch_effects"]["n_batches"]
        batch_var = noise_cfg["batch_effects"]["batch_variance"]
        # Assign samples to batches
        batch_assignments = rng.choice(n_batches, size=n_samples)
        # Generate batch effects
        batch_multipliers = rng.lognormal(0, batch_var, size=n_batches)
        sample_batch_effects = batch_multipliers[batch_assignments]
        typer.echo(f"📊 Applied batch effects ({n_batches} batches, var={batch_var:.1%})")
    
    for s in range(n_samples):
        # multiply by fold-change if sample belongs to a treatment group
        is_treat = (groups[s] != "Control")
        baseline_this = baseline * (fc_vec if is_treat else 1)
        
        # Apply PCR bias
        baseline_this = baseline_this * pcr_bias
        
        # Generate base abundances
        rel_abund = dirichlet_sample(baseline_this, conc_param, rng)
        counts = rng.multinomial(lib_sizes[s], rel_abund)
        
        # Apply sample processing noise
        if noise_cfg.get("sample_processing", {}).get("enabled", False):
            proc_cv = noise_cfg["sample_processing"]["cv"]
            proc_sigma = np.sqrt(np.log(1 + proc_cv**2))
            proc_noise = rng.lognormal(0, proc_sigma)
            counts = (counts * proc_noise).astype(int)
        
        # Apply batch effects
        counts = (counts * sample_batch_effects[s]).astype(int)
        
        weight_mat[s] = counts

    # Post-processing: apply contamination and zero-inflation
    
    # Contamination (add low-level contaminant reads)
    if noise_cfg.get("contamination", {}).get("enabled", False):
        n_contam = noise_cfg["contamination"]["n_contaminants"]
        contam_rate = noise_cfg["contamination"]["contamination_rate"]
        # Select random genomes as contaminants
        contam_idx = rng.choice(n_genomes, size=n_contam, replace=False)
        for s in range(n_samples):
            total_reads = weight_mat[s].sum()
            contam_reads = int(total_reads * contam_rate)
            # Distribute contamination reads among contaminant genomes
            contam_counts = rng.multinomial(contam_reads, 
                                          np.ones(n_contam) / n_contam)
            weight_mat[s, contam_idx] += contam_counts
        typer.echo(f"📊 Added contamination ({n_contam} species, {contam_rate:.1%} rate)")
    
    # Zero-inflation (random dropouts)
    if noise_cfg.get("zero_inflation", {}).get("enabled", False):
        dropout_prob = noise_cfg["zero_inflation"]["dropout_prob"]
        # Generate dropout mask (1 = keep, 0 = dropout)
        dropout_mask = rng.binomial(1, 1-dropout_prob, size=(n_samples, n_genomes))
        weight_mat = weight_mat * dropout_mask
        n_dropouts = np.sum(dropout_mask == 0)
        typer.echo(f"📊 Applied zero-inflation ({n_dropouts} dropouts, {dropout_prob:.1%} prob)")

    # 7 ▸ write read_depth
    with (out_dir / "read_depth.txt").open("w") as fh:
        for sid, d in zip(samples, lib_sizes):
            print(sid, d, sep="\t", file=fh)
    
    # Write batch assignments if batch effects are enabled
    if noise_cfg.get("batch_effects", {}).get("enabled", False):
        with (out_dir / "batch_assignments.txt").open("w") as fh:
            print("sample_id", "batch", sep="\t", file=fh)
            for sid, batch in zip(samples, batch_assignments):
                print(sid, f"batch_{batch+1}", sep="\t", file=fh)

    # 8 ▸ sample_distributions.tsv
    header = ["genomes", "seq_type", *samples]
    with (out_dir / "sample_distributions.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for j, g in enumerate(genomes):
            w.writerow([g.resolve(), "chromosome", *weight_mat[:, j]])

    # 9 ▸ truth table
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
