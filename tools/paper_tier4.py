#!/usr/bin/env python3
"""Tier 4 cohort replication (IBD + CRC): reliability, not accuracy.

Without ground truth, report only cross-cohort replication of each method's
discoveries (rank agreement, replicated-discovery rate between independent
cohorts, per method) and agreement between methods. Also report the raw-data
null from Tier 3a (permuted labels in the same cohort). Deferred in the
reduced scope; tooling frozen here.

  python tools/paper_tier4.py --a cohortA_discoveries.tsv --b cohortB_discoveries.tsv
"""
from pathlib import Path
import argparse
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from magician.replication import replication_rate, cross_cohort_table


def ranked_list(path):
    df = pd.read_csv(path, sep="\t")
    if "qvalue" in df.columns:
        df = df.sort_values("qvalue")
    return df["feature_id"].tolist() if "feature_id" in df.columns else []


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--a", default=None)
    p.add_argument("--b", default=None)
    p.add_argument("--method", default="method")
    p.add_argument("--out", default="analysis/tier4_replication.tsv")
    args = p.parse_args()
    if not (args.a and args.b):
        print("Tier 4 deferred in reduced scope (need 2 independent cohorts). "
              "Provide --a/--b discovery tables to compute replication.")
        print("Metric: Jaccard replicated-discovery rate + Spearman rank agreement.")
        return 0
    res = replication_rate(ranked_list(args.a), ranked_list(args.b))
    res["method"] = args.method
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([res]).to_csv(out, sep="\t", index=False)
    print(f"replication: jaccard={res['jaccard']:.3f} spearman={res['spearman']:.3f} -> {out}")
    print("NOTE: replication is reliability, not accuracy (must state in paper).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
