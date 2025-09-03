#!/usr/bin/env python3
"""
make_sample_distribution_from_real.py – Build CAMISIM design from real taxa mapping
──────────────────────────────────────────────────────────────────────────────────
Inputs
  - data/metaphlan4_genome_mapping.csv  (columns: original_metaphlan_taxon, species_search_term,
                                         abundance, relative_abundance, genome_status,
                                         genome_accession, genome_filename)
  - --genome_dir pointing to FASTA files (compressed .fna.gz allowed; will be decompressed if needed)

Outputs (under --out_dir)
  - sample_distributions.tsv   (CAMISIM abundance weights per sample)
  - truth_table.csv            (DA status and intended log2FC per genome)
  - read_depth.txt             (per-sample read depth for CAMISIM)

Design choices
  - Base composition p0 from real relative_abundance for rows with usable genomes
  - 40 samples total (default): 20 Control, 20 Treatment
  - One replicate per sample (names without _rep suffix)
  - Per-sample compositions sampled via Dirichlet(alpha = c * p0), c configurable
  - Differential abundance: balanced up/down across abundance strata; modest effect sizes
  - Convert compositions to CAMISIM weights with inverse genome-size scaling; integer weights
  - Tight library size band to avoid confounding
"""

import argparse
import csv
import gzip
import math
import pathlib
import random
import shutil
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate CAMISIM design from real mapping")
    p.add_argument("--mapping_csv", type=pathlib.Path, default=pathlib.Path("data/metaphlan4_genome_mapping.csv"),
                   help="CSV with taxon to genome mapping and relative_abundance")
    p.add_argument("--genome_dir", type=pathlib.Path, default=pathlib.Path("/local/mrahman/magician/data/mooh_representative"),
                   help="Directory containing genome FASTA (.fna or .fna.gz)")
    p.add_argument("--out_dir", type=pathlib.Path, default=pathlib.Path("data/from_real_v1"),
                   help="Output directory for design files")

    # Sampling & groups
    p.add_argument("--n_controls", type=int, default=20)
    p.add_argument("--n_treats", type=int, default=20)
    p.add_argument("--seed", type=int, default=1337)
    p.add_argument("--dirichlet_concentration", type=float, default=400.0,
                   help="Dirichlet total concentration c; alpha = c * p0")

    # DA configuration
    p.add_argument("--n_da", type=int, default=30, help="Total number of DA genomes (balanced up/down)")
    p.add_argument("--effect_sizes_up", nargs="*", type=float, default=[1.2, 1.3, 1.5, 1.6],
                   help="Fold-changes for up-regulated features")
    p.add_argument("--effect_sizes_down", nargs="*", type=float, default=[1/1.2, 1/1.3, 1/1.5],
                   help="Fold-changes for down-regulated features (must be < 1)")

    # Library sizes & weights
    p.add_argument("--depth_min", type=float, default=28e6)
    p.add_argument("--depth_max", type=float, default=32e6)
    p.add_argument("--min_weight", type=int, default=5, help="Minimum integer weight per genome per sample")

    return p.parse_args()


def _decompress_if_needed(fpath: pathlib.Path) -> pathlib.Path:
    """If fpath endswith .gz, produce an uncompressed .fna next to it and return the .fna path."""
    if fpath.suffix == ".gz":
        dest = fpath.with_suffix("")  # strip .gz
        if not dest.exists():
            try:
                with gzip.open(fpath, "rb") as fin, open(dest, "wb") as fout:
                    shutil.copyfileobj(fin, fout)
            except Exception as e:
                print(f"⚠️  Failed to decompress {fpath.name}: {e}")
                return fpath
        return dest
    return fpath


def _fasta_total_size(fpath: pathlib.Path) -> int:
    total = 0
    try:
        for rec in SeqIO.parse(fpath, "fasta"):
            total += len(rec.seq)
    except Exception as e:
        print(f"⚠️  Could not read {fpath.name}: {e}")
        return 1_000_000
    return total if total > 0 else 1_000_000


def load_mapping(mapping_csv: pathlib.Path, genome_dir: pathlib.Path) -> pd.DataFrame:
    df = pd.read_csv(mapping_csv)
    # Keep rows with usable genome_filename and non-failed status
    usable = (
        df["genome_filename"].notna()
        & (df["genome_filename"].astype(str) != "Not Available")
        & (~df["genome_status"].astype(str).str.contains("Download Failed|Species Not Extractable", case=False))
    )
    sub = df.loc[usable, ["original_metaphlan_taxon", "relative_abundance", "genome_filename"]].copy()
    # Deduplicate by genome_filename (sum relative_abundance)
    sub = (
        sub.groupby(["genome_filename"], as_index=False)
           .agg({
               "relative_abundance": "sum",
               # Keep one representative taxon (first) for provenance
               "original_metaphlan_taxon": "first",
           })
    )
    sub["genome_path"] = sub["genome_filename"].apply(lambda n: (genome_dir / str(n)).resolve())
    # Filter to existing files
    exists_mask = sub["genome_path"].apply(lambda p: p.exists())
    if int((~exists_mask).sum()) > 0:
        missing = sub.loc[~exists_mask, "genome_filename"].head(5).tolist()
        print(f"⚠️  Missing {int((~exists_mask).sum())} genomes. Examples: {missing}")
    sub = sub.loc[exists_mask].reset_index(drop=True)
    # Normalize relative abundance to sum to 1 (across available genomes)
    total = sub["relative_abundance"].sum()
    if total <= 0:
        raise SystemExit("❌ relative_abundance sums to zero – cannot proceed")
    sub["p0"] = sub["relative_abundance"] / total
    return sub


def pick_da_indices(p0: np.ndarray, n_da: int, rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray]:
    """Abundance-stratified selection: choose ~n_da/2 up and ~n_da/2 down across deciles."""
    n = len(p0)
    n_da = min(n_da, n)
    ranks = np.argsort(p0)  # ascending
    # Split into 10 strata (deciles)
    strata = np.array_split(ranks, 10)
    up, down = [], []
    target_up = n_da // 2
    target_down = n_da - target_up
    # Alternate drawing from each stratum to spread selection
    s_idx = 0
    while (len(up) < target_up or len(down) < target_down) and s_idx < 10000:
        for s in strata:
            if len(up) < target_up and len(s) > 0:
                up.append(int(rng.choice(s)))
            if len(down) < target_down and len(s) > 0:
                down.append(int(rng.choice(s)))
            if len(up) >= target_up and len(down) >= target_down:
                break
        s_idx += 1
    return np.array(up[:target_up]), np.array(down[:target_down])


def main() -> None:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    random.seed(args.seed)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # 1) Load mapping and build base composition p0 over available genomes
    map_df = load_mapping(args.mapping_csv, args.genome_dir)
    genomes: List[pathlib.Path] = []
    p0_list: List[float] = []
    for _, row in map_df.iterrows():
        fpath = pathlib.Path(row["genome_path"])  # may be .gz
        fpath_unz = _decompress_if_needed(fpath)
        genomes.append(fpath_unz)
        p0_list.append(float(row["p0"]))
    p0 = np.array(p0_list, dtype=float)
    p0 = p0 / p0.sum()

    n_genomes = len(genomes)
    print(f"📊 Using {n_genomes} genomes from mapping (after filtering & existence check)")

    # 2) Compute genome sizes for inverse-size scaling
    sizes = np.array([_fasta_total_size(g) for g in genomes], dtype=float)
    med_size = float(np.median(sizes))
    size_factors = med_size / sizes

    # 3) Build sample names
    ctrl_ids = [f"C{str(i).zfill(2)}" for i in range(1, args.n_controls + 1)]
    trt_ids  = [f"T{str(i).zfill(2)}" for i in range(1, args.n_treats + 1)]
    samples  = ctrl_ids + trt_ids

    # 4) Select DA indices (balanced up/down)
    up_idx, down_idx = pick_da_indices(p0, args.n_da, rng)
    up_idx_set = set(up_idx.tolist())
    down_idx_set = set(down_idx.tolist())

    # 5) Sample per-sample compositions and apply DA to treatment
    alpha = p0 * float(args.dirichlet_concentration)
    comp_by_sample: Dict[str, np.ndarray] = {}
    for sid in samples:
        comp = rng.dirichlet(alpha)
        # Apply DA effects to treatment samples only
        if sid.startswith("T"):
            # Up-regulated
            if len(up_idx) > 0:
                ups = rng.choice(np.array(args.effect_sizes_up, dtype=float), size=len(up_idx), replace=True)
                comp[up_idx] = comp[up_idx] * ups
            # Down-regulated (apply multiplicative factor < 1)
            if len(down_idx) > 0:
                dns = rng.choice(np.array(args.effect_sizes_down, dtype=float), size=len(down_idx), replace=True)
                comp[down_idx] = comp[down_idx] * dns
            # Renormalize
            s = comp.sum()
            if s <= 0:
                comp = p0.copy()
            else:
                comp = comp / s
        comp_by_sample[sid] = comp

    # 6) Convert compositions to CAMISIM integer weights with inverse-size scaling
    # Choose a per-sample scale so typical weight magnitudes are in tens
    def to_weights(comp: np.ndarray) -> np.ndarray:
        scaled = comp * size_factors
        # Normalize to sum 1, then scale by a constant factor
        s = scaled.sum()
        if s <= 0:
            scaled = comp
            s = scaled.sum()
        scaled = scaled / s
        weights = np.rint(np.maximum(args.min_weight, scaled * 2000)).astype(int)
        return weights

    # 7) Write read_depth.txt
    depths = rng.uniform(args.depth_min, args.depth_max, size=len(samples)).round().astype(int)
    with (args.out_dir / "read_depth.txt").open("w") as fh:
        for sid, d in zip(samples, depths):
            print(sid, d, sep="\t", file=fh)

    # 8) Write sample_distributions.tsv
    header = ["genomes", "seq_type", *samples]
    with (args.out_dir / "sample_distributions.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for gi, g in enumerate(genomes):
            row = [str(g.resolve()), "chromosome"]
            for sid in samples:
                comp = comp_by_sample[sid]
                row.append(int(to_weights(comp)[gi]))
            w.writerow(row)

    # 9) Write truth_table.csv – DA status and intended log2FC (group-level intent)
    # For up_idx, assign a representative fold-change as median of effect_sizes_up; for down likewise
    up_fc_repr = float(np.median(np.array(args.effect_sizes_up, dtype=float))) if len(args.effect_sizes_up) else 1.0
    down_fc_repr = float(np.median(np.array(args.effect_sizes_down, dtype=float))) if len(args.effect_sizes_down) else 1.0

    with (args.out_dir / "truth_table.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["genome", "is_DA", "true_log2FC", "genome_size_bp", "size_scaling_factor", "direction"])
        for i, g in enumerate(genomes):
            if i in up_idx_set:
                is_da = True
                log2fc = math.log2(up_fc_repr)
                direction = "up"
            elif i in down_idx_set:
                is_da = True
                log2fc = math.log2(down_fc_repr)
                direction = "down"
            else:
                is_da = False
                log2fc = 0.0
                direction = "none"
            w.writerow([g.name, is_da, log2fc, int(sizes[i]), float(size_factors[i]), direction])

    # Summary
    print("✅ From-real design written to", args.out_dir.resolve())
    print(f"   Samples: {len(samples)}  (Controls={args.n_controls}, Treatment={args.n_treats})")
    print(f"   DA genomes: {args.n_da} (balanced up/down)")
    print(f"   Dirichlet concentration: {args.dirichlet_concentration}")
    print(f"   Genome size range: {int(sizes.min()):,}–{int(sizes.max()):,} bp; median={int(med_size):,} bp")


if __name__ == "__main__":
    main()




