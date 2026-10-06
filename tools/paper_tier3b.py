#!/usr/bin/env python3
"""Tier 3b read-level implants (Nearing-style implants extended to MAG space).

In group-B samples, up/down-sample reads aligned to chosen MAGs by a known
factor (4 effects x 3 DA fractions x 10 seeds = 120 DA sets on one assembly),
conserving library size. Truth exact by construction; assembly/bins/noise/batch
real. Deferred in the reduced scope (Tier 1+2+3a first); tooling is frozen here
so thresholds cannot be tuned on it later.

  python tools/paper_tier3b.py --counts counts.tsv --mags MAG_00001,MAG_00002 \\
      --factor 2.0 --fraction 0.1 --seeds 10 --out cache/paper_tier3b
"""
from pathlib import Path
import argparse
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from magician.replication import implant_reads
from magician.io import save_table, write_json


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts", required=True)
    p.add_argument("--features", default=None)
    p.add_argument("--mags", default="", help="Comma-separated MAG ids to implant")
    p.add_argument("--factor", type=float, default=2.0)
    p.add_argument("--fraction", type=float, default=0.1)
    p.add_argument("--seeds", type=int, default=10)
    p.add_argument("--out", default="cache/paper_tier3b")
    args = p.parse_args()
    counts = pd.read_csv(args.counts, sep="\t", index_col=0)
    mags = [m for m in args.mags.split(",") if m] or list(counts.index[: max(2, int(len(counts) * args.fraction))])
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    listing = []
    for seed in range(args.seeds):
        implanted, truth = implant_reads(counts, mags, args.factor, seed)
        case = out / f"implanted_s{seed}"
        case.mkdir(parents=True, exist_ok=True)
        implanted.reset_index().to_csv(case / "counts.tsv", sep="\t", index=False)
        truth.rename("is_da").reset_index().rename(
            columns={"index": "feature_id"}).to_csv(case / "truth.tsv", sep="\t", index=False)
        listing.append(dict(case=case.name, factor=args.factor, seed=seed, n_da=int(truth.sum())))
    save_table(pd.DataFrame(listing), out / "cases.tsv")
    write_json(out / "designs.json", {"tier": "3b", "effects": [0.5, 1, 2, 4],
                                      "fractions": [0.05, 0.10, 0.20], "note": "FROZEN - do not tune on this"})
    print(f"wrote {len(listing)} implant sets to {out} (reduced scope: tooling frozen, run deferred)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
