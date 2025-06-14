#!/usr/bin/env python3
"""
Create sample_distributions.tsv with 24 genomes × 20 samples
(10 controls C01-C10, 10 treatments T01-T10)
"""

import pathlib, random, csv, os

# -------- user-editable settings ----------
reference_dir = pathlib.Path("data/camisim_toy/")   # folder with *.fa or *.gb
n_controls    = 5
n_treats      = 10  
FC_DA         = 4          # fold-change for differential genomes
seed          = 42         # reproducibility
# -----------------------------------------

random.seed(seed)

# 1) collect genome files
genomes = sorted(reference_dir.glob("*.*"))         # all files
assert len(genomes) == 24, f"Expected 24 genomes, found {len(genomes)}"

# 2) choose which genomes are differential
da_indices = set(random.sample(range(24), k=4))      # 4 / 24 ≈ 17 %

# 3) build header
header = ["genomes", "seq_type"] \
       + [f"C{str(i).zfill(2)}" for i in range(1, n_controls+1)] \
       + [f"T{str(i).zfill(2)}" for i in range(1, n_treats+1)]

out_tsv = pathlib.Path("data/sample_distributions_camisim_toy.tsv")
out_tsv.parent.mkdir(exist_ok=True, parents=True)

with out_tsv.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(header)

    for idx, g in enumerate(genomes):
        baseline = 1
        treat    = FC_DA if idx in da_indices else 1
        row = [g.resolve(), "chromosome"] \
            + [baseline]*n_controls \
            + [treat]*n_treats
        writer.writerow(row)

print("✅  sample_distributions.tsv written to", out_tsv.resolve())
print("   Differential genomes (4):", [genomes[i].name for i in da_indices])
