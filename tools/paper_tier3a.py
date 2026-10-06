#!/usr/bin/env python3
"""Tier 3a real nulls: 200 random label splits on one real assembly (any discovery false).

Real steps (Phase 0 picks the cohort; see tools/pick_datasets.py):
  1. Take healthy-cohort metagenomes (>=100 samples, one study/site/similar depth).
  2. Subsample reads to a common depth, assemble + bin once (existing workflow).
  3. This tool builds a matrix bundle with 200 cases sharing one count matrix
     but different group labels (balanced; confounded with batch where metadata
     allows). Any discovery is false: type I error with real structure.

  python tools/paper_tier3a.py --counts counts.tsv --features features.tsv \\
      --samples samples.tsv --out cache/paper_tier3a --splits 200 --seed 0
  python run_magician.py --configfile config/paper_tier3a.yaml --cores 60

Smoke test (no real reads needed):
  python tools/paper_tier3a.py --smoke --cores 4
"""
from pathlib import Path
import argparse
import subprocess
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from magician.replication import random_splits
from magician.io import save_table, write_json


def build_null_bundle(counts_file, features_file, samples_file, output, n_splits=200, seed=0):
    counts = pd.read_csv(counts_file, sep="\t")
    features = pd.read_csv(features_file, sep="\t")
    samples = pd.read_csv(samples_file, sep="\t")
    if "sample_id" not in samples.columns:
        raise ValueError("samples.tsv needs sample_id")
    batch = None
    if "batch" in samples.columns:
        batch = dict(zip(samples.sample_id, samples.batch))
    ids = samples.sample_id.tolist()
    # Subsample reads to a common depth is done upstream; here we only relabel.
    splits = random_splits(ids, n_splits, seed, batch)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    listing = []
    feature_ids = counts["feature_id"].tolist() if "feature_id" in counts.columns else features["feature_id"].tolist()
    for i, labels in enumerate(splits):
        case = f"null_s{seed + i}"
        dest = out / case
        dest.mkdir(parents=True, exist_ok=True)
        counts.to_csv(dest / "counts.tsv", sep="\t", index=False)
        features.to_csv(dest / "features.tsv", sep="\t", index=False)
        rows = [dict(sample_id=s, group=labels[s],
                     read_pairs=int(samples.loc[samples.sample_id == s, "read_pairs"].iloc[0])
                     if "read_pairs" in samples.columns else 1000000)
                for s in ids]
        save_table(pd.DataFrame(rows), dest / "samples.tsv")
        truth = pd.DataFrame({"feature_id": feature_ids, "is_da": False,
                              "true_log2fc": 0.0,
                              "baseline_proportion": 1.0 / len(feature_ids),
                              "treatment_proportion": 1.0 / len(feature_ids),
                              "length_bp": 12000})
        truth.to_csv(dest / "truth.tsv", sep="\t", index=False)
        # Latent expectations: uniform (null truth needs no effect sizes).
        exp_rows = []
        for s in ids:
            for f in feature_ids:
                exp_rows.append(dict(sample_id=s, feature_id=f, genomic_proportion=1.0 / len(feature_ids),
                                     target_proportion=1.0 / len(feature_ids),
                                     read_proportion=1.0 / len(feature_ids),
                                     library_read_pairs=rows[0]["read_pairs"],
                                     length_bp=12000, expected_count=100,
                                     absolute_copies=float("nan")))
        save_table(pd.DataFrame(exp_rows), dest / "expectations.tsv")
        pd.DataFrame({"feature_id": feature_ids, "source_id": feature_ids,
                      "match_status": "matched"}).to_csv(dest / "matching.tsv", sep="\t", index=False)
        listing.append(dict(case=case, scenario="null", seed=seed + i, design="tier3a_real_null"))
    save_table(pd.DataFrame(listing), out / "cases.tsv")
    write_json(out / "designs.json", {"tier": "3a", "n_cases": len(listing), "splits": n_splits, "seed": seed,
                                      "note": "Any discovery is false; real structure at no assumption cost."})
    return listing


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts", default=None)
    p.add_argument("--features", default=None)
    p.add_argument("--samples", default=None)
    p.add_argument("--out", default="cache/paper_tier3a")
    p.add_argument("--splits", type=int, default=200)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--cores", type=int, default=4)
    p.add_argument("--generate-only", action="store_true")
    args = p.parse_args()
    if args.smoke:
        # Stand-in: reuse the Tier 1 smoke bundle's first case counts.
        import tempfile
        from magician.generators import g1_case
        tmp = Path(args.out + "_smoke_src")
        tmp.mkdir(parents=True, exist_ok=True)
        spec = dict(samples_per_group=10, n_features=50, depth=20000, fraction_da=0.0,
                    log2fc=1.0, direction="balanced", total_load="none",
                    zero_structure="sampling", confounder="none")
        g1_case(spec, "null", 999, tmp / "donor")
        listing = build_null_bundle(tmp / "donor" / "counts.tsv", tmp / "donor" / "features.tsv",
                                    tmp / "donor" / "samples.tsv", args.out, 8, args.seed)
        print(f"smoke Tier 3a bundle: {len(listing)} null splits in {args.out}")
        if args.generate_only:
            return 0
        cmd = [sys.executable, str(ROOT / "run_magician.py"), "--configfile",
               "config/paper_tier3a.yaml", "--cores", str(args.cores)]
        print("+", " ".join(cmd), "(matrix bundle must be cache/paper_tier3a; "
              "smoke wrote 8 splits there)")
        return subprocess.call(cmd, cwd=str(ROOT))
    if not (args.counts and args.features and args.samples):
        p.error("Provide --counts/--features/--samples (from one real assembly) or --smoke")
    listing = build_null_bundle(args.counts, args.features, args.samples, args.out,
                                args.splits, args.seed)
    print(f"Tier 3a bundle: {len(listing)} null splits in {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
