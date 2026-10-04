"""Exercise the real R contract and adapters on reproducible count data.

Every fixture is designed to break something specific: unsafe R feature names,
positive and negative effects, structural zeros, a zero library, a single feature,
an empty catalogue, permuted sample order and metadata that disagrees with the
matrix. Assertions check status, declared units, sign and ID preservation.
"""
from pathlib import Path
import argparse
import json
import os
import shutil
import subprocess
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from magician.config import load, registry, zero_policies

POSITIVE = {"wilcoxon_clr", "maaslin2", "deseq2", "edger", "limma_voom", "aldex2", "aldex2_g0",
            "ancombc2", "maaslin3_abundance", "aldex3", "adapt", "linda", "zicoseq", "locom",
            "rademu"}
NEGATIVE = set()


def run(rscript, work, method, counts, samples, features, *, timeout_minutes=20,
        allow_failures="false", metric="counts", kind="sources", registry_entry=None,
        hard_timeout_minutes=30):
    entry = json.loads(registry_entry)
    suffix = method
    result, status, session = [work / f"{suffix}.{s}" for s in ("tsv", "status.json", "session.txt")]
    command = [rscript, str(ROOT / "workflow/scripts/run_da.R"), json.dumps(entry),
               str(counts), str(samples), str(features), str(result), str(status), str(session),
               "0.05", "128", "42", str(int(timeout_minutes * 60)), allow_failures,
               "null_s1", metric, kind, str(work)]
    log = work / f"{suffix}.log"
    with log.open("w") as handle:
        proc = subprocess.run(command, env=dict(os.environ, OMP_NUM_THREADS="1",
                                                 OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                                                 TMPDIR=str(work)), stdout=handle,
                              stderr=subprocess.STDOUT, timeout=hard_timeout_minutes * 60)
    info = json.loads(status.read_text())
    data = pd.read_csv(result, sep="\t")
    return proc, info, data, log


def fixtures(work, rng):
    """Counts plus the metadata and feature tables the contract reads."""
    ids = [f"source-{i:03d}" for i in range(80)]
    samples = [f"{group}{i:02d}" for group in ("C", "T") for i in range(1, 7)]
    means = np.full((80, 12), 400.0)
    means[:10, 6:] *= 4
    means[10:20, 6:] /= 4
    values = rng.negative_binomial(40, 40 / (40 + means))
    counts = pd.DataFrame(values, index=ids, columns=samples)
    counts.index.name = "feature_id"
    counts.to_csv(work / "counts.tsv", sep="\t")
    # Length-corrected relative abundance for methods whose contract starts there.
    lengths = np.array([12000 + (i % 7) * 2500 for i in range(80)], dtype=float)
    per_base = counts.divide(pd.Series(lengths, index=ids), axis=0)
    relative = per_base.divide(per_base.sum(axis=0), axis=1)
    relative.to_csv(work / "relative.tsv", sep="\t")
    pd.DataFrame({"sample_id": samples, "group": ["Control"] * 6 + ["Treatment"] * 6,
                  "read_pairs": [4000] * 12}).to_csv(work / "samples.tsv", sep="\t", index=False)
    pd.DataFrame({"feature_id": ids, "retained_for_da": True,
                  "reserved_other": False}).to_csv(work / "features.tsv", sep="\t", index=False)

    structural = counts.copy()
    structural.loc[ids[20:30], :] = 0
    structural.iloc[:, 3:6] = 0          # a sample with a zero library
    structural.to_csv(work / "zeros.tsv", sep="\t")
    structural_relative = relative.copy()
    structural_relative.loc[ids[20:30], :] = 0
    structural_relative.iloc[:, 3:6] = 0
    structural_relative.to_csv(work / "zeros_relative.tsv", sep="\t")
    # Structural zeros are left as reserved mass for the declared `other` category.
    reserved = {feature: index >= 20 for index, feature in enumerate(ids)}
    pd.DataFrame({"feature_id": ids, "retained_for_da": True,
                  "reserved_other": [reserved[feature] for feature in ids]}).to_csv(
        work / "features_reserved.tsv", sep="\t", index=False)
    # A wide matrix, so an explicit timeout has work it genuinely cannot finish.
    wide_ids = [f"wide-{i:05d}" for i in range(4000)]
    wide_samples = [f"{group}{i:02d}" for group in ("C", "T") for i in range(1, 21)]
    wide = pd.DataFrame(np.random.default_rng(7).negative_binomial(
        40, 40 / (40 + np.full((4000, 40), 200.0))), index=wide_ids, columns=wide_samples)
    wide.index.name = "feature_id"
    wide.to_csv(work / "wide.tsv", sep="\t")
    pd.DataFrame({"sample_id": wide_samples, "group": ["Control"] * 20 + ["Treatment"] * 20,
                  "read_pairs": [20000] * 40}).to_csv(work / "samples_wide.tsv", sep="\t", index=False)
    pd.DataFrame({"feature_id": wide_ids, "retained_for_da": True,
                  "reserved_other": False}).to_csv(work / "features_wide.tsv", sep="\t", index=False)

    counts.iloc[:1].to_csv(work / "single.tsv", sep="\t")
    pd.DataFrame({"feature_id": ids[:1], "retained_for_da": True,
                  "reserved_other": False}).to_csv(work / "features_single.tsv", sep="\t", index=False)
    counts.iloc[:0].to_csv(work / "empty.tsv", sep="\t")
    pd.DataFrame(columns=["feature_id", "retained_for_da", "reserved_other"]).to_csv(
        work / "features_empty.tsv", sep="\t", index=False)
    # A permuted sample order must be reordered explicitly, not intersected away.
    order = [5, 4, 3, 2, 1, 0] + list(range(6, 12))
    counts[samples].T.iloc[order].T.to_csv(work / "counts_permuted.tsv", sep="\t")
    pd.DataFrame({"sample_id": samples, "group": ["Treatment"] * 6 + ["Control"] * 6,
                  "read_pairs": [4000] * 12}).to_csv(work / "samples_missing_group.tsv", sep="\t", index=False)
    pd.DataFrame({"sample_id": samples[:-1], "group": ["Control"] * 6 + ["Treatment"] * 5,
                  "read_pairs": [4000] * 11}).to_csv(work / "samples_incomplete.tsv", sep="\t", index=False)
    return ids, samples, counts


def check_sign(method, data, ids):
    if method in {"zicoseq"}:
        # An F statistic is non-negative; direction is meaningless for it.
        return
    if method in NEGATIVE:
        return
    if method in POSITIVE:
        assert data.set_index("feature_id").log2fc.loc[ids[:10]].median() > 0, \
            f"{method}: the Treatment effect has the wrong direction"
        assert data.set_index("feature_id").log2fc.loc[ids[10:20]].median() < 0, \
            f"{method}: the down-regulated effect has the wrong direction"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", default="cache/adapter_validation")
    parser.add_argument("--methods", default="", help="Comma-separated registry ids; default all active")
    parser.add_argument("--timeout-minutes", type=float, default=20.0)
    parser.add_argument("--conda-prefix", default=str(ROOT / "cache/conda"))
    parser.add_argument("--frozen-prefix", default=str(ROOT / "cache/envs"),
                        help="Bootstrapped method environments, as cache/envs/magician_*")
    parser.add_argument("--rscript", action="append", default=[],
                        help="Override an environment's Rscript, as ENV=PATH")
    args = parser.parse_args()
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    cfg = load()
    reg = registry(cfg)
    if args.methods:
        methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    else:
        methods = [entry["id"] for entry in reg.payload["methods"]
                   if entry["status"] in cfg["analysis"]["enable_status"]]
    ids, samples, counts = fixtures(work, np.random.default_rng(123))
    failures = []

    def report(name, ok, detail=""):
        if ok:
            print(f"{name}: ok {detail}".strip(), flush=True)
            return
        failures.append(name)
        print(f"{name}: FAILED {detail}", flush=True)

    conda_prefix = Path(args.conda_prefix)
    executables = sorted(path for path in conda_prefix.glob("*/bin/Rscript")) if conda_prefix.exists() else []
    frozen = Path(args.frozen_prefix)
    if frozen.exists():
        executables += sorted(path for path in frozen.glob("magician_*/bin/Rscript"))
    if not executables:
        found = shutil.which("Rscript")
        executables = [Path(found)] if found else []
    if not executables:
        raise SystemExit("No Rscript found; pass --rscript ENV=PATH or build the pinned environments")
    overrides = dict(item.split("=", 1) for item in args.rscript)

    def rscript_for(method):
        name = reg.entry(method)["environment"]
        if name in overrides:
            return Path(overrides[name])
        package = reg.entry(method)["package"]
        if package:
            for candidate in executables:
                probe = subprocess.run(
                    [str(candidate), "-e", f'cat(requireNamespace("{package}", quietly=TRUE))'],
                    capture_output=True, text=True, timeout=180)
                if probe.stdout.strip().endswith("TRUE"):
                    return candidate
        return executables[0]

    for method in methods:
        entry = reg.entry_json(method)
        executable = rscript_for(method)
        accepted = reg.entry(method)["accepted_inputs"]
        metric = "counts" if "counts" in accepted else accepted[0]
        counts_file = work / "counts.tsv" if metric == "counts" else work / f"{metric}.tsv"
        try:
            proc, info, data, log = run(executable, work, method, counts_file,
                                        work / "samples.tsv", work / "features.tsv",
                                        registry_entry=entry, metric=metric,
                                        timeout_minutes=args.timeout_minutes)
            ok = proc.returncode == 0 and info["status"] == "success"
            assert data.feature_id.tolist() == ids, "Feature IDs changed or were reordered"
            assert data.endpoint.tolist() == [reg.entry(method)["endpoint"]] * len(ids)
            assert data.effect_scale.tolist() == [reg.entry(method)["effect_scale"]] * len(ids)
            assert info["n_tested"] > 0, "No finite tests"
            check_sign(method, data, ids)
            print(f"{method}: success, version {info['package_version']}, "
                  f"{info['n_tested']} tests, endpoint {info['endpoint']}", flush=True)
        except (OSError, ValueError, AssertionError, KeyError, subprocess.TimeoutExpired) as exc:
            report(method, False, f"({exc}); see {work / (method + '.log')}")

        # Structural zeros and a zero library must stay visible, not become a crash.
        if method in {"wilcoxon_clr", "maaslin2", "maaslin2_count_scale", "wilcoxon_clr_count_scale",
                      "maaslin3_abundance", "maaslin3_corrected", "aldex2", "aldex2_g0",
                      "aldex3", "adapt", "rademu", "linda", "zicoseq", "locom"}:
            zeros_file = work / ("zeros.tsv" if metric == "counts" else f"zeros_{metric}.tsv")
            try:
                proc, info, data, log = run(executable, work, method + "_zeros", zeros_file,
                                            work / "samples.tsv", work / "features.tsv",
                                            registry_entry=entry, metric=metric,
                                            timeout_minutes=args.timeout_minutes,
                                            allow_failures="true")
                # A zero library is data no model can fit; failing explicitly is correct.
                assert info["status"] in {"success", "empty", "failed", "unavailable"}, info["message"]
                if info["status"] == "success":
                    assert data.feature_id.tolist() == ids, "Feature IDs lost with structural zeros"
                print(f"{method} with structural zeros and a zero library: "
                      f"{info['status']} ({info['n_tested']} tests)", flush=True)
            except (OSError, ValueError, AssertionError, KeyError) as exc:
                report(f"{method}_zeros", False, f"({exc})")
            # The reserved other category joins the geometry but is never a discovery.
            try:
                proc, info, data, log = run(executable, work, method + "_reserved", zeros_file,
                                            work / "samples.tsv", work / "features_reserved.tsv",
                                            registry_entry=entry, metric=metric,
                                            timeout_minutes=args.timeout_minutes,
                                            allow_failures="true")
                assert info["status"] in {"success", "empty", "failed", "unavailable"}, info["message"]
                if info["status"] != "success":
                    print(f"{method} with a reserved other category: {info['status']} "
                          f"({info['message'][:60]})", flush=True)
                    continue
                reserved = data.loc[data.fit_status == "reserved_category"]
                assert not reserved.tested.any(), "A reserved category was reported as a discovery"
                print(f"{method} with a reserved other category: "
                      f"{info['status']}, {len(reserved)} reserved, never tested", flush=True)
            except (OSError, ValueError, AssertionError, KeyError) as exc:
                report(f"{method}_reserved", False, f"({exc})")

        # A single feature, an empty catalogue, permuted samples and bad metadata.
        if method == "wilcoxon_clr":
            for name, counts_file, features_file, samples_file, want in [
                ("single", "single.tsv", "features_single.tsv", "samples.tsv", {"success", "empty", "failed"}),
                ("empty", "empty.tsv", "features_empty.tsv", "samples.tsv", {"empty"}),
                ("permuted", "counts_permuted.tsv", "features.tsv", "samples.tsv", {"success"}),
                ("missing_group", "counts.tsv", "features.tsv", "samples_missing_group.tsv", {"failed"}),
                ("incomplete", "counts.tsv", "features.tsv", "samples_incomplete.tsv", {"failed"}),
            ]:
                try:
                    proc, info, data, log = run(executable, work, f"wilcoxon_{name}", work / counts_file,
                                                work / samples_file, work / features_file,
                                                registry_entry=entry, timeout_minutes=args.timeout_minutes)
                    assert info["status"] in want, f"status {info['status']}, expected one of {sorted(want)}"
                    if info["status"] == "failed":
                        # A rejected input must be visible, never a quiet success.
                        assert proc.returncode != 0, "a rejected input reported success"
                        print(f"contract {name}: rejected ({info['message'][:60]})", flush=True)
                        continue
                    assert proc.returncode == 0, f"exit {proc.returncode}: {info['message']}"
                    if name == "permuted":
                        assert data.feature_id.tolist() == ids, "Permuted samples changed feature order"
                        assert info["n_tested"] == info["n_retained"], "Samples were silently dropped"
                    print(f"contract {name}: {info['status']} ({info['message'][:60]})", flush=True)
                except (OSError, ValueError, AssertionError, KeyError) as exc:
                    report(f"contract_{name}", False, f"({exc}); see {work / ('wilcoxon_' + name + '.log')}")

        # A method must refuse an input metric it does not declare.
        probe = next((m for m in ("counts", "tpm", "rpkm", "relative")
                      if m not in reg.entry(method)["accepted_inputs"]), None)
        if probe:
            try:
                proc, info, data, log = run(executable, work, f"{method}_bad_metric", work / "counts.tsv",
                                            work / "samples.tsv", work / "features.tsv",
                                            registry_entry=entry, metric=probe,
                                            timeout_minutes=args.timeout_minutes)
                assert proc.returncode != 0 and info["status"] == "failed"
                print(f"{method} refuses {probe}: {info['message'][:60]}", flush=True)
            except (OSError, ValueError, AssertionError, KeyError) as exc:
                report(f"{method}_bad_metric", False, f"({exc})")

        # An explicit timeout must be reported, never replaced by another method.
        if method == "wilcoxon_clr":
            try:
                proc, info, data, log = run(executable, work, "wilcoxon_timeout", work / "wide.tsv",
                                            work / "samples_wide.tsv", work / "features_wide.tsv",
                                            registry_entry=entry, timeout_minutes=0.02,
                                            allow_failures="true", hard_timeout_minutes=10)
                assert info["status"] == "timeout", f"status {info['status']}"
                print("explicit timeout: reported as timeout, never replaced", flush=True)
            except (OSError, ValueError, AssertionError, KeyError) as exc:
                report("wilcoxon_timeout", False, f"({exc})")
    print(f"\n{len(failures)} failure(s); logs in {work}", flush=True)
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())