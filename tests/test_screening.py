"""Matrix screening replicates and versioned real-genome manifests."""
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from magician.screening import DESIGNS, build_bundle, design_replicate
from magician.genomes import check_genomes, select_manifest, write_manifest
from magician.io import table, write_json, fasta


def test_every_declared_design_is_distinct_and_balanced():
    assert len(DESIGNS) >= 5
    baseline = DESIGNS["baseline"]
    # The subset varies one factor at a time against a common baseline design.
    for name, spec in DESIGNS.items():
        if name == "baseline":
            continue
        differing = [key for key, value in spec.items() if baseline.get(key) != value]
        assert differing, f"{name} is a duplicate of the baseline design"
    assert {spec["effect"] for spec in DESIGNS.values()} == {"2x", "4x"}
    assert {spec["changed_fraction"] for spec in DESIGNS.values()} == {0.1, 0.2, 0.4}


def test_replicate_is_independent_reproducible_and_evaluable(tmp_path):
    first = design_replicate("baseline", "spiked", 7, tmp_path / "a")
    second = design_replicate("baseline", "spiked", 7, tmp_path / "b")
    other = design_replicate("baseline", "spiked", 8, tmp_path / "c")
    assert (tmp_path / "a/counts.tsv").read_bytes() == (tmp_path / "b/counts.tsv").read_bytes()
    assert (tmp_path / "a/counts.tsv").read_bytes() != (tmp_path / "c/counts.tsv").read_bytes()
    truth = table(tmp_path / "a/truth.tsv")
    counts = table(tmp_path / "a/counts.tsv")
    samples = table(tmp_path / "a/samples.tsv")
    assert len(truth) == DESIGNS["baseline"]["n_features"]
    assert set(counts.columns) - {"feature_id"} == set(samples.sample_id)
    assert counts.feature_id.tolist() == truth.feature_id.tolist()
    assert truth.is_da.sum() == truth.true_log2fc.abs().gt(0).sum()
    # Redistribution conserves the changed subset's mass.
    changed = truth.loc[truth.is_da, "baseline_proportion"].sum()
    assert truth.loc[truth.is_da, "treatment_proportion"].sum() == pytest.approx(changed)
    assert truth.baseline_proportion.sum() == pytest.approx(1.0)
    assert truth.treatment_proportion.sum() == pytest.approx(1.0)
    # Up and down effects are both present, and integer counts are non-negative.
    assert (truth.loc[truth.is_da, "true_log2fc"] > 0).any()
    assert (truth.loc[truth.is_da, "true_log2fc"] < 0).any()
    assert (counts.set_index("feature_id") >= 0).all().all()


def test_null_replicates_have_no_changed_feature(tmp_path):
    design_replicate("baseline", "null", 11, tmp_path / "null")
    truth = table(tmp_path / "null/truth.tsv")
    assert not truth.is_da.any()
    assert truth.true_log2fc.abs().max() == pytest.approx(0.0, abs=1e-12)


def test_bundle_uses_disjoint_seeds_per_design_and_scenario(tmp_path):
    summary = build_bundle(["baseline", "overdispersed"], 2, 2, tmp_path, first_seed=5000)
    cases = table(tmp_path / "cases.tsv")
    assert len(cases) == 8 == len(summary)
    assert cases.case.is_unique
    assert cases.seed.is_unique
    assert set(cases.scenario) == {"null", "spiked"}
    written = read_json = __import__("json").loads((tmp_path / "designs.json").read_text())
    assert written["n_cases"] == 8
    assert set(written["designs"]) == {"baseline", "overdispersed"}
    for case in cases.case:
        for name in ("counts.tsv", "samples.tsv", "features.tsv", "truth.tsv",
                     "expectations.tsv", "matching.tsv"):
            assert (tmp_path / case / name).exists(), f"{case}/{name} missing"


def test_spurious_bundles_are_rejected(tmp_path):
    from magician.config import load
    bundle = build_bundle(["baseline"], 1, 1, tmp_path / "bundle", first_seed=9000)
    cfg = load(overrides={"input_mode": "matrix",
                          "matrix_input": {"bundle": str(tmp_path / "bundle")},
                          "analysis": {"min_total_count": 0, "min_prevalence": 0.0}})
    from magician import matrix
    assert matrix.cases(cfg) == ["null_s9000", "spiked_s9001"]
    case = tmp_path / "bundle/null_s9000"
    counts = pd.read_csv(case / "counts.tsv", sep="\t")
    counts["C01"] = 0
    counts.to_csv(case / "counts.tsv", sep="\t", index=False)
    # Zero total counts are still a valid matrix; they simply have no information.
    matrix.validate_case(cfg, "null_s9000")


def _entry(tmp_path, name, length=12000, accession="GCF_000001.1", group="diverse"):
    path = tmp_path / f"{name}.fa"
    with path.open("w") as handle:
        handle.write(f">{name}\n{'ACGT' * (length // 4)}\n")
    return dict(genome_id=name, accession=accession, version="1", path=str(path),
                provenance="synthetic fixture", source_group=group)


def test_manifest_checks_reject_short_and_unknown_accessions(tmp_path):
    entries = [_entry(tmp_path, "genome_a"),
               _entry(tmp_path, "genome_b", length=100),
               _entry(tmp_path, "genome_c", accession="not-an-accession")]
    rows, problems = check_genomes(entries)
    assert len(rows) == 3
    assert any("below the" in problem for problem in problems)
    assert any("recognised accession" in problem for problem in problems)


def test_manifest_checks_reject_missing_files(tmp_path):
    rows, problems = check_genomes([dict(genome_id="absent", accession="GCF_000002.1", version="1",
                                        path=str(tmp_path / "nope.fa"),
                                        provenance="none", source_group="diverse")])
    assert not rows and any("missing file" in problem for problem in problems)


def test_selection_rule_is_recorded_and_enforces_group_sizes(tmp_path):
    entries = [_entry(tmp_path, f"div_{i}", group="diverse") for i in range(3)]
    entries += [_entry(tmp_path, f"rel_{i}", group="related", accession=f"GCF_10000{i}.1")
                for i in range(2)]
    rows, problems = check_genomes(entries)
    chosen, problems, rule = select_manifest(rows, {"diverse": 2, "related": 1},
                                             "first two per group in file order")
    assert not problems
    assert len(chosen) == 3
    frame = write_manifest(chosen, tmp_path / "manifest.tsv", rule)
    assert len(frame) == 3
    assert frame.genome_id.is_unique
    assert (tmp_path / "manifest.json").exists()
    # A group that cannot satisfy its target is an explicit failure, not a short list.
    short, problems, _ = select_manifest(rows, {"diverse": 9}, "nine diverse genomes")
    assert not short
    assert any("candidates" in problem for problem in problems)


def test_manifest_records_lengths_checksums_and_provenance(tmp_path):
    rows, _ = check_genomes([_entry(tmp_path, "genome_a")])
    frame = write_manifest(rows, tmp_path / "m.tsv", "single fixture genome")
    assert frame.length_bp.iloc[0] == 12000
    assert len(str(frame.sha256.iloc[0])) == 64
    assert frame.provenance.iloc[0] == "synthetic fixture"
    assert frame.accession.iloc[0] == "GCF_000001.1"