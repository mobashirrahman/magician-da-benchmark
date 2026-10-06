#!/usr/bin/env python3
"""Phase 0 dataset picker for Tiers 3-4 (default: auto-pick IBD + CRC).

Dataset accessions are NOT pinned in plan.md. Phase 0 selects and records them
with URL, licence, checksum and read counts after checking they are public and
redistributable. Nothing is downloaded by this tool; it records the choice.

  python tools/pick_datasets.py --list
  python tools/pick_datasets.py --choose ibd_default crc_default --out config/datasets.yaml

Defaults (suggestion, unverified): one IBD set (IBDMDB/HMP2 + independent IBD)
and one CRC set (>=3 public cohorts); Tier 3a needs >=100 healthy samples from
one study/site/similar depth. Verify before running: this tool only scaffolds
the record that Phase 0 must fill in.
"""
from pathlib import Path
import argparse
import sys
import yaml

ROOT = Path(__file__).resolve().parents[1]

CANDIDATES = {
    "ibd_default": {"disease": "IBD", "cohorts": ["IBDMDB/HMP2", "independent IBD cohort (TBD)"],
                    "url": "TBD in Phase 0", "licence": "TBD", "checksum": "TBD",
                    "note": "Cross-cohort replication; no ground truth (Tier 4)."},
    "crc_default": {"disease": "CRC", "cohorts": ["3+ public CRC cohorts (TBD)"],
                    "url": "TBD in Phase 0", "licence": "TBD", "checksum": "TBD",
                    "note": "Cross-cohort replication; no ground truth (Tier 4)."},
    "healthy_default": {"disease": "healthy", "cohorts": ["one study, one body site, >=100 samples (TBD)"],
                        "url": "TBD in Phase 0", "licence": "TBD", "checksum": "TBD",
                        "note": "Tier 3a real nulls + 3b implants; subsampled to common depth."},
    "mock_default": {"disease": "mock", "cohorts": ["ZymoBIOMICS-type log ratios, public CAMI data (TBD)"],
                     "url": "TBD in Phase 0", "licence": "TBD", "checksum": "TBD",
                     "note": "Tier 3c sanity anchor; never used for ranking."},
}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--list", action="store_true")
    p.add_argument("--choose", nargs="*", default=["healthy_default", "ibd_default", "crc_default"])
    p.add_argument("--out", default="config/datasets.yaml")
    args = p.parse_args()
    if args.list:
        for name, rec in CANDIDATES.items():
            print(f"{name}: {rec['disease']} {rec['cohorts']} -- {rec['note']}")
        return 0
    chosen = {}
    for name in args.choose:
        if name not in CANDIDATES:
            p.error(f"Unknown candidate: {name} (see --list)")
        chosen[name] = CANDIDATES[name]
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump({"datasets": chosen,
                                   "warning": "Accessions unverified; Phase 0 must fill URL/licence/checksum/read counts."},
                                  sort_keys=False))
    print(f"wrote {out} with {len(chosen)} dataset slots (UNVERIFIED - Phase 0 must verify)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
