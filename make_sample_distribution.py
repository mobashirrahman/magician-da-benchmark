#!/usr/bin/env python3
"""
make_design.py  –  realistic CAMISIM abundance table + truth map
─────────────────────────────────────────────────────────────────
• 24 genomes               (all FASTA/GBK files in reference_dir)
• 10 Control + 10 Treatment samples
• Baseline weight per genome ~ Γ(shape=2, scale=5)  → median ~9
• 4 genomes are differential: multiplied by 2×, 4× or 8×
• Each sample pulls a random library size 15–60 M reads
  (written to a separate metadata file – CAMISIM’s read_depth.txt)
Outputs
  - sample_distributions.tsv          (for CAMISIM --metadata)
  - truth_table.csv                   (is_DA, true_log2FC, baseline, fc)
  - read_depth.txt                    (sample depth in reads)
"""

import pathlib, random, csv, math, json
import numpy as np

# ─── user-tunable knobs ──────────────────────────────────────────────
reference_dir   = pathlib.Path("data/camisim_toy")       # 24 genome FASTAs
n_controls      = 10
n_treats        = 10
gamma_shape     = 2
gamma_scale     = 5
fold_choices    = [2, 4, 8]      # sampled w/out replacement if possible
n_differential  = 4
depth_min, depth_max = 15e6, 60e6   # per-sample read depth range
seed            = 42
out_prefix      = pathlib.Path("data/camisim_realistic")
# ─────────────────────────────────────────────────────────────────────

rng = np.random.default_rng(seed)
random.seed(seed)
out_prefix.mkdir(parents=True, exist_ok=True)

# 1 ▸ collect genomes
genomes = sorted(reference_dir.glob("*.*"))
if len(genomes) != 24:
    raise SystemExit(f"❌ expected 24 genomes but found {len(genomes)} in {reference_dir}")

n_genomes = len(genomes)

# 2 ▸ baseline weights ~ Gamma
baseline_wt = rng.gamma(gamma_shape, gamma_scale, size=n_genomes).round().astype(int)
baseline_wt[baseline_wt < 5] = 5    # never below 5

# 3 ▸ choose DA genomes and assign fold-changes
da_idx = rng.choice(n_genomes, size=n_differential, replace=False)
folds  = rng.choice(fold_choices, size=n_differential, replace=True)
enriched_wt = baseline_wt.copy()
enriched_wt[da_idx] = baseline_wt[da_idx] * folds

# 4 ▸ library sizes per sample
samples_all = [f"C{str(i).zfill(2)}" for i in range(1, n_controls+1)] + \
              [f"T{str(i).zfill(2)}" for i in range(1, n_treats+1)]
depths = rng.uniform(depth_min, depth_max, size=len(samples_all)).round().astype(int)

# 5 ▸ write read_depth.txt (CAMISIM takes this directly)
with (out_prefix / "read_depth.txt").open("w") as fh:
    for sid, d in zip(samples_all, depths):
        print(sid, d, sep="\t", file=fh)

# 6 ▸ sample_distributions.tsv
header = ["genomes", "seq_type"] + samples_all
with (out_prefix / "sample_distributions.tsv").open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(header)
    for i, g in enumerate(genomes):
        weights = (
            [baseline_wt[i]] * n_controls +
            [enriched_wt[i] if i in da_idx else baseline_wt[i]] * n_treats
        )
        w.writerow([g.resolve(), "chromosome", *weights])

# 7 ▸ ground truth table
with (out_prefix / "truth_table.csv").open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["genome", "is_DA", "true_log2FC",
                "baseline_weight", "enriched_weight"])
    for i, g in enumerate(genomes):
        is_da = i in da_idx
        fc    = enriched_wt[i] / baseline_wt[i]
        log2fc = math.log2(fc)
        w.writerow([g.name, is_da, log2fc,
                    baseline_wt[i], enriched_wt[i]])

print("✅ design written to", out_prefix.resolve())
print("   differential genomes:",
      ", ".join(genomes[i].name for i in da_idx))
