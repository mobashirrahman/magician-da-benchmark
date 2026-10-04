#!/usr/bin/env python3
"""Check, select and freeze a real-genome manifest.

  python tools/real_genomes.py check config/my_genomes.tsv
  python tools/real_genomes.py freeze config/my_genomes.tsv --out config/genome_manifest.tsv \\
      --rule "one genome per source group, longest first" --group diverse=20

Nothing is downloaded. Candidate paths may point at an existing local reference
collection; only the chosen entries are referenced, by path.
"""
from pathlib import Path
import argparse
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from magician.genomes import check_genomes, select_manifest, write_manifest
from magician.io import table


def parse_groups(values):
    groups = {}
    for item in values or []:
        if "=" not in item:
            raise SystemExit(f"--group expects NAME=COUNT, got {item!r}")
        name, count = item.split("=", 1)
        groups[name] = int(count)
    return groups


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["check", "freeze"])
    parser.add_argument("candidates", help="TSV with genome_id, accession, version, path, provenance, source_group")
    parser.add_argument("--out", default="config/genome_manifest.tsv")
    parser.add_argument("--rule", default="all candidates in file order")
    parser.add_argument("--group", action="append", default=[],
                        help="NAME=COUNT selection target, repeatable")
    parser.add_argument("--no-require-groups", action="store_true")
    args = parser.parse_args()
    candidates = table(Path(args.candidates))
    required = {"genome_id", "accession", "version", "path", "provenance", "source_group"}
    missing = required - set(candidates.columns)
    if missing:
        raise SystemExit(f"{args.candidates} is missing columns: {', '.join(sorted(missing))}")
    rows, problems = check_genomes(candidates.to_dict("records"))
    for problem in problems:
        print(f"unusable: {problem}", file=sys.stderr)
    print(f"{len(rows) - len(problems)} of {len(candidates)} candidates are usable")
    if args.command == "check":
        return 1 if problems else 0
    groups = parse_groups(args.group)
    if not groups:
        chosen, problems, rule = rows, [], args.rule
    else:
        chosen, problems, rule = select_manifest(rows, groups, args.rule,
                                                 require_groups=not args.no_require_groups)
    for problem in problems:
        print(f"selection: {problem}", file=sys.stderr)
    if not chosen:
        raise SystemExit("Selection produced no genomes")
    frame = write_manifest(chosen, Path(args.out), rule)
    print(f"wrote {args.out} with {len(frame)} genomes "
          f"({int(frame.length_bp.min())} to {int(frame.length_bp.max())} bp)")
    return 0


if __name__ == "__main__":
    sys.exit(main())