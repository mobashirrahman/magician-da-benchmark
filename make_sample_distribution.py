#!/usr/bin/env python3
"""
make_design.py  –  realistic CAMISIM abundance table + truth map
─────────────────────────────────────────────────────────────────
• genomes: all FASTA/GBK files found in `reference_dir` (no fixed count)
• 10 Control + 10 Treatment samples
• Baseline weight per genome ~ Γ(shape=2, scale=5)  → median ~9
• 4 genomes are differential: multiplied by 2×, 4× or 8×
• Each sample pulls a random library size 15–60 M reads
  (written to a separate metadata file – CAMISIM's read_depth.txt)
Outputs
  - sample_distributions.tsv          (for CAMISIM --metadata)
  - truth_table.csv                   (is_DA, true_log2FC, baseline, fc)
  - read_depth.txt                    (sample depth in reads)
"""

import pathlib, random, csv, math, json, gzip, shutil
import numpy as np

# ─── user-tunable knobs ──────────────────────────────────────────────
# Path to the folder containing downloaded genomes (adjust if needed)
reference_dir   = pathlib.Path("data/HOMD/genomes")
n_controls      = 10
n_treats        = 10
n_replicates    = 2  # number of technical/biological replicates per sample
gamma_shape     = 2
gamma_scale     = 5
fold_choices    = [2, 4, 8]      # sampled w/out replacement if possible
n_differential  = 80
depth_min, depth_max = 15e6, 60e6   # per-sample read depth range
seed            = 42
out_prefix      = pathlib.Path("data/HOMD/design")
# ─────────────────────────────────────────────────────────────────────

rng = np.random.default_rng(seed)
random.seed(seed)
out_prefix.mkdir(parents=True, exist_ok=True)

# 1 ▸ decompress any *.gz genome files first so CAMISIM can read them
for gz_path in reference_dir.glob("*.gz"):
    dest_path = gz_path.with_suffix("")  # strip the .gz extension
    if not dest_path.exists():
        print(f"⚙️  Decompressing {gz_path.name} → {dest_path.name}")
        with gzip.open(gz_path, "rb") as f_in, open(dest_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    # Remove the original compressed file to avoid duplicate inputs
    try:
        gz_path.unlink()
    except Exception as e:
        print(f"⚠️  Could not remove {gz_path.name}: {e}")

# 2 ▸ collect genomes (exclude any remaining .gz)
genomes = sorted(p for p in reference_dir.iterdir() if p.is_file() and not p.suffix.endswith('gz'))
if not genomes:
    raise SystemExit(f"❌ No genome FASTA files found in {reference_dir} – aborting.")

n_genomes = len(genomes)

# ensure we don't request more differential genomes than available
n_differential = min(n_differential, n_genomes)

# 3 ▸ baseline weights ~ Gamma
baseline_wt = rng.gamma(gamma_shape, gamma_scale, size=n_genomes).round().astype(int)
baseline_wt[baseline_wt < 5] = 5    # never below 5

# 4 ▸ choose DA genomes and assign fold-changes
da_idx = rng.choice(n_genomes, size=n_differential, replace=False)
folds  = rng.choice(fold_choices, size=n_differential, replace=True)
enriched_wt = baseline_wt.copy()
enriched_wt[da_idx] = baseline_wt[da_idx] * folds

# 5 ▸ library sizes per sample
samples_base = [f"C{str(i).zfill(2)}" for i in range(1, n_controls+1)] + \
               [f"T{str(i).zfill(2)}" for i in range(1, n_treats+1)]

# explode base IDs into replicate-specific IDs, e.g. C01_rep01 … C01_rep05
samples_all = [f"{sid}_rep{str(r).zfill(2)}"
               for sid in samples_base
               for r in range(1, n_replicates + 1)]

# draw one library size per replicate
depths = rng.uniform(depth_min, depth_max, size=len(samples_all)).round().astype(int)

# 6 ▸ write read_depth.txt (CAMISIM takes this directly)
with (out_prefix / "read_depth.txt").open("w") as fh:
    for sid, d in zip(samples_all, depths):
        print(sid, d, sep="\t", file=fh)

# 7 ▸ sample_distributions.tsv
header = ["genomes", "seq_type"] + samples_all
with (out_prefix / "sample_distributions.tsv").open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(header)
    for i, g in enumerate(genomes):
        weights = (
            [baseline_wt[i]] * n_controls * n_replicates +
            [enriched_wt[i] if i in da_idx else baseline_wt[i]] * n_treats * n_replicates
        )
        w.writerow([g.resolve(), "chromosome", *weights])

# 8 ▸ ground truth table
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
