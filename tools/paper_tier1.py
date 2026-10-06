#!/usr/bin/env python3
"""Build and run a versioned, corrected Tier 1 benchmark.

  python tools/paper_tier1.py --smoke --cores 16
  python tools/paper_tier1.py --full --cores 16
  python tools/paper_tier1.py --full --generate-only

Smoke covers all three generators and balanced/all-up/all-down implants.
Full uses 24 spiked core cells (100 replicates each) on 6 null conditions (200
replicates each) and a one-factor-at-a-time survey (20 spiked per cell, 20 null
per new null condition). Null data do not depend on effect factors, so null
cases are generated once per null condition. Batch-confounding scenarios are
excluded. Nonempty destinations are rejected; earlier experiment data is preserved.
"""
from pathlib import Path
import argparse
import hashlib
import subprocess
import sys
import shutil
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from magician.generators import (GENERATOR_VERSION, PROFILE_DIR, design_plan, tier1_plan,
                                 build_bundle, realism_check)
from magician.config import resolve
from magician.io import table, write_json


def main():
    p = argparse.ArgumentParser(description=__doc__)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--smoke', action='store_true')
    mode.add_argument('--full', action='store_true')
    p.add_argument('--bundle')
    p.add_argument('--output-dir')
    p.add_argument('--config-out')
    p.add_argument('--null-replicates', type=int)
    p.add_argument('--spiked-replicates', type=int)
    p.add_argument('--fractional-replicates', type=int, default=20)
    p.add_argument('--first-seed', type=int, default=1000)
    p.add_argument('--cores', type=int, default=16)
    p.add_argument('--scheduler', choices=['greedy', 'ilp'], default='greedy')
    p.add_argument('--generate-only', action='store_true')
    p.add_argument('--no-realism', action='store_true')
    args = p.parse_args()
    if not 1 <= args.cores <= 60:
        p.error('cores must be within [1, 60]')
    suffix = '_validation' if args.smoke else ''
    bundle = resolve(args.bundle or f'cache/paper_tier1_v3{suffix}')
    output = args.output_dir or f'results_paper_tier1_v3{suffix}'
    config_path = resolve(args.config_out or f'config/paper_tier1_v3{suffix}.yaml')
    if bundle.exists() and any(bundle.iterdir()):
        p.error(f'Bundle {bundle} is not empty; choose a new versioned directory')
    if config_path.exists() or resolve(output).exists():
        p.error('Config/results destination already exists; choose new destinations')
    for replicas in (args.null_replicates, args.spiked_replicates, args.fractional_replicates):
        if replicas is not None and not 1 <= replicas < 100000:
            p.error('replicate counts must be within [1, 99999]')
    config = yaml.safe_load((ROOT / 'config/paper_tier1.yaml').read_text())
    config['output_dir'] = output
    config['matrix_input']['bundle'] = str(bundle)
    if args.full:
        # Retain both input generations and all durable analysis outputs.
        config['storage']['budget_gb'] = 200
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(f'# Corrected Tier 1 ({GENERATOR_VERSION}); written by tools/paper_tier1.py.\n' + yaml.safe_dump(config, sort_keys=False))
    if args.smoke:
        cells = [dict(generator=g, samples_per_group=10, n_features=100, depth=500000,
                      fraction_da=0.10, log2fc=1.0, direction=d, total_load='none',
                      zero_structure='sampling', confounder='none')
                 for g in ('g1', 'g2', 'g3') for d in ('balanced', 'all-up', 'all-down')]
        plan = design_plan(cells, args.null_replicates or 2, args.spiked_replicates or 2, 'validation')
    else:
        plan = tier1_plan(args.null_replicates or 200, args.spiked_replicates or 100,
                          args.fractional_replicates)
    build_bundle(plan, bundle, args.first_seed, progress=True)
    sources = sorted((ROOT / 'src/magician').glob('*.py'))
    sources += sorted((ROOT / 'workflow/scripts').rglob('*.R'))
    sources += sorted((ROOT / 'workflow/scripts').rglob('*.py'))
    sources += [ROOT / 'tools/paper_tier1.py', ROOT / 'tools/paper_analysis.py', ROOT / 'workflow/Snakefile',
                ROOT / 'tools/fetch_real_profiles.py', *sorted(PROFILE_DIR.iterdir())]
    sources += sorted((ROOT / 'workflow/rules').glob('*.smk'))
    sources += sorted((ROOT / 'config').glob('*.yaml'))
    sources += sorted((ROOT / 'workflow/schemas').glob('*.yaml'))
    sources += sorted((ROOT / 'workflow/envs').glob('*.yaml'))
    sources += sorted((ROOT / 'workflow/envs/locks').glob('*.json'))
    sources += [ROOT / 'run_magician.py', ROOT / 'pyproject.toml', ROOT / 'environment.yaml']
    hashes = {}
    for source in sources:
        relative = source.relative_to(ROOT)
        target = bundle / 'code_snapshot' / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = hashlib.sha256(source.read_bytes()).hexdigest()
    write_json(bundle / 'code_snapshot.json', hashes)
    git = lambda *a: subprocess.run(['git', *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    write_json(bundle / 'revision.json', dict(git_revision=git('rev-parse', 'HEAD'),
                                              git_describe=git('describe', '--tags', '--always'),
                                              uncommitted_changes=bool(git('status', '--porcelain', '--untracked-files=no'))))
    if not args.no_realism:
        realism_check(bundle, bundle / 'diagnostics')
    print(f'{GENERATOR_VERSION}: {len(table(bundle / "cases.tsv"))} cases; config {config_path}', flush=True)
    command = [sys.executable, str(ROOT / 'run_magician.py'), '--configfile', str(config_path),
               '--cores', str(args.cores), '--scheduler', args.scheduler]
    print('Run:', ' '.join(command), flush=True)
    if args.generate_only:
        return 0
    status = subprocess.call(command, cwd=ROOT)
    if status:
        return status
    return subprocess.call([sys.executable, str(ROOT / 'tools/paper_analysis.py'),
                            '--tier1', str(resolve(output) / 'benchmark'),
                            '--out', str(resolve(output) / 'analysis')], cwd=ROOT)


if __name__ == '__main__':
    raise SystemExit(main())
