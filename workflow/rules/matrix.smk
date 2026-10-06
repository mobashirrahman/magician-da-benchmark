"""Matrix-only entry point: import a versioned bundle, then use the shared analysis path."""
MATRIX_FILES = ["samples.tsv", "features.tsv", "counts.tsv", "truth.tsv",
                "expectations.tsv", "matching.tsv"]

rule import_matrix_case:
    input:
        listing=lambda wc: [str(REPO / CFG["matrix_input"]["bundle"] / "cases.tsv")],
        bundle=lambda wc: ([str(REPO / CFG["matrix_input"]["bundle"] / wc.case / name)
                           for name in MATRIX_FILES] +
                          ([str(REPO / CFG["matrix_input"]["bundle"] / wc.case / "generator.json")]
                           if (REPO / CFG["matrix_input"]["bundle"] / wc.case / "generator.json").exists() else []))
    output:
        raw=directory(f"{OUT}/matrices/{{case}}/raw_sources"),
        samples=f"{OUT}/design/{{case}}/samples.tsv",
        expectations=f"{OUT}/design/{{case}}/expectations.tsv",
        truth=f"{OUT}/design/{{case}}/truth.tsv",
        design=f"{OUT}/design/{{case}}/design.json",
        matching=f"{OUT}/catalogues/{{case}}/matching.tsv"
    params:
        implementation=implementation("matrix"),
        cfg=CFG, step="import_matrix", scratch=lambda wc: f"{OUT}/work/{wc.case}/scratch"
    threads: 1
    resources: mem_mb=1024, disk_mb=200
    log: f"{OUT}/logs/{{case}}/import_matrix.log"
    benchmark: f"{OUT}/benchmarks/design/{{case}}.tsv"
    conda: PYENV
    script: STEP

rule import_matrix_manifest:
    input:
        listing=lambda wc: [str(REPO / CFG["matrix_input"]["bundle"] / "cases.tsv")],
        bundles=lambda wc: [str(REPO / CFG["matrix_input"]["bundle"] / wc.case / name)
                            for wc.case in CASES for name in MATRIX_FILES]
    output:
        manifest=f"{OUT}/provenance/matrix_bundle.json"
    params:
        implementation=implementation("matrix"),
        cfg=CFG, step="import_manifest", scratch=f"{OUT}/work/scratch/report"
    threads: 1
    resources: mem_mb=256, disk_mb=50
    log: f"{OUT}/logs/import_matrix_manifest.log"
    benchmark: f"{OUT}/benchmarks/import_matrix.tsv"
    conda: PYENV
    script: STEP
