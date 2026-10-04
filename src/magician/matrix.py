"""Matrix-only execution mode: quantify, filter and test matrices without reads.

A bundle is a versioned directory of cases. Each case supplies raw counts, feature
lengths, sample metadata, the latent expectations the truth lanes are computed
from, and an identity MAG/source matching table. Nothing here needs FASTQ, an
assembly or a binning catalogue, so new methods and matrix simulations never pay
for the genome pipeline.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from .config import case_info, resolve
from .io import table, save_table, sha256, write_json
from .truth import read_expectations

REQUIRED = ("samples.tsv", "features.tsv", "counts.tsv", "truth.tsv",
            "expectations.tsv", "matching.tsv")


def bundle_root(cfg):
    path = resolve(cfg["matrix_input"]["bundle"])
    if not path.is_dir():
        raise ValueError(f"Matrix bundle directory not found: {path}")
    return path


def cases(cfg):
    """Case identifiers declared by the bundle, in file order."""
    listing = table(bundle_root(cfg) / "cases.tsv")
    if not {"case", "scenario", "seed"} <= set(listing.columns):
        raise ValueError("cases.tsv must contain case, scenario and seed columns")
    names = listing.case.astype(str).tolist()
    if len(set(names)) != len(names):
        raise ValueError("cases.tsv lists a case more than once")
    for name in names:
        case_info(name)
    return names


def digest(cfg):
    root = bundle_root(cfg)
    files = sorted(p for p in root.rglob("*") if p.is_file())
    return {str(p.relative_to(root)): sha256(p) for p in files}


def validate_case(cfg, case):
    """Reject a bundle case whose parts disagree before any job is built."""
    source = bundle_root(cfg) / case
    missing = [name for name in REQUIRED if not (source / name).exists()]
    if missing:
        raise ValueError(f"Matrix bundle case {case} is missing: {', '.join(missing)}")
    samples = table(source / "samples.tsv")
    features = table(source / "features.tsv")
    counts = table(source / "counts.tsv")
    truth = table(source / "truth.tsv")
    expectations = read_expectations(source / "expectations.tsv")
    matching = table(source / "matching.tsv")
    if not {"sample_id", "group", "read_pairs"} <= set(samples.columns):
        raise ValueError(f"{case}: samples.tsv needs sample_id, group and read_pairs")
    if not samples.sample_id.is_unique:
        raise ValueError(f"{case}: duplicate sample IDs")
    if set(samples.group.astype(str)) != {"Control", "Treatment"}:
        raise ValueError(f"{case}: exactly Control and Treatment samples are required")
    if not features.feature_id.is_unique:
        raise ValueError(f"{case}: duplicate feature IDs")
    counts = counts.set_index("feature_id")
    if set(counts.columns) != set(samples.sample_id):
        raise ValueError(f"{case}: count columns and sample IDs disagree")
    if set(counts.index) - set(features.feature_id):
        raise ValueError(f"{case}: counts contain unknown feature IDs")
    if set(truth.feature_id) - set(features.feature_id):
        raise ValueError(f"{case}: truth contains unknown feature IDs")
    if expectations.duplicated(["sample_id", "feature_id"]).any():
        raise ValueError(f"{case}: duplicated sample and genome rows in expectations.tsv")
    if set(expectations.sample_id) - set(samples.sample_id):
        raise ValueError(f"{case}: expectations contain unknown sample IDs")
    if set(expectations.feature_id) - set(features.feature_id):
        raise ValueError(f"{case}: expectations contain unknown feature IDs")
    if set(matching.feature_id) - set(features.feature_id):
        raise ValueError(f"{case}: matching contains unknown feature IDs")
    if (counts < 0).any().any():
        raise ValueError(f"{case}: counts must be nonnegative")
    return dict(case=case, n_features=len(features), n_samples=len(samples))


def import_case(cfg, case, output):
    """Validate the bundle and place its durable raw counts and latent truth."""
    source = bundle_root(cfg) / case
    validate_case(cfg, case)
    out = Path(output)
    (out / "design" / case).mkdir(parents=True, exist_ok=True)
    (out / "catalogues" / case).mkdir(parents=True, exist_ok=True)
    raw = out / "matrices" / case / "raw_sources"
    raw.mkdir(parents=True, exist_ok=True)
    samples = table(source / "samples.tsv")
    features = table(source / "features.tsv")
    expectations = read_expectations(source / "expectations.tsv")
    counts = table(source / "counts.tsv").set_index("feature_id").reindex(features.feature_id)
    counts.index.name = "feature_id"
    counts = counts.apply(pd.to_numeric, errors="raise")
    if (counts % 1 != 0).any().any():
        raise ValueError(f"{case}: a matrix bundle must supply integer counts")
    counts = counts.astype(int)
    counts.reset_index().to_csv(raw / "counts.tsv", sep="\t", index=False, na_rep="NA")
    features.to_csv(raw / "features.tsv", sep="\t", index=False, na_rep="NA")
    totals = counts.sum(axis=0)
    library = pd.DataFrame({"sample_id": samples.sample_id,
                            "read_pairs": samples.read_pairs.astype(int),
                            "mapped_pairs": totals.reindex(samples.sample_id).to_numpy(),
                            "catalogue_pairs": totals.reindex(samples.sample_id).to_numpy(),
                            "unassigned_pairs": 0,
                            "unmapped_pairs": (samples.read_pairs.astype(int)
                                               - totals.reindex(samples.sample_id).to_numpy())})
    library.to_csv(raw / "library.tsv", sep="\t", index=False, na_rep="NA")
    write_json(raw / "raw.json", {"kind": "sources", "n_features": int(len(features)),
                                  "n_samples": int(len(samples)), "source": "matrix bundle",
                                  "count_unit": "integer counts as supplied by the bundle",
                                  "raw_retained": True})
    for name in ("samples.tsv", "truth.tsv"):
        table(source / name).to_csv(out / "design" / case / name, sep="\t", index=False, na_rep="NA")
    expectations.rename(columns={"feature_id": "genome_id"}).to_csv(
        out / "design" / case / "expectations.tsv", sep="\t", index=False, na_rep="NA")
    table(source / "matching.tsv").to_csv(out / "catalogues" / case / "matching.tsv",
                                          sep="\t", index=False, na_rep="NA")
    write_json(out / "design" / case / "design.json",
               dict(case=case, scenario=case_info(case)[0], seed=case_info(case)[1],
                    input_mode="matrix", bundle=str(source), config=cfg))


def import_manifest(cfg, output):
    """Immutable record of exactly which bundle files this run consumed."""
    checksums = digest(cfg)
    write_json(Path(output), {"input_mode": "matrix",
                              "bundle": str(bundle_root(cfg)),
                              "n_cases": len(cases(cfg)),
                              "files": checksums})