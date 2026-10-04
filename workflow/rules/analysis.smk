rule quantify_features:
    input:
        counts=lambda wc: expand(f"{OUT}/work/{{case}}/counts/{{sample}}.{'sources' if wc.kind == 'sources' else 'contigs'}.tsv", case=wc.case, sample=SAMPLES),
        samples=f"{OUT}/design/{{case}}/samples.tsv",
        features=lambda wc: f"{OUT}/provenance/references.tsv" if wc.kind == "sources" else f"{OUT}/catalogues/{wc.case}/mags.tsv",
        assignments=lambda wc: [f"{OUT}/catalogues/{wc.case}/assignments.tsv"] if wc.kind == "mags" else []
    output:
        raw=directory(f"{OUT}/matrices/{{case}}/raw_{{kind}}")
    params:
        implementation=implementation("matrices"),
        cfg=MATRIX_CFG, step="quantify", scratch=lambda wc: f"{OUT}/work/{wc.case}/scratch"
    threads: 1
    resources: mem_mb=1024, disk_mb=200
    log: f"{OUT}/logs/{{case}}/quantify/{{kind}}.log"
    benchmark: f"{OUT}/benchmarks/{{case}}/quantify/{{kind}}.tsv"
    conda: PYENV
    script: STEP

# Filtering and normalization are derived from the durable raw counts, so changing a
# threshold re-runs only this rule and the analysis, never simulation or assembly.
rule derive_matrices:
    input:
        raw=f"{OUT}/matrices/{{case}}/raw_{{kind}}",
        settings=REGISTRY,
        policies=POLICY_FILE
    output:
        matrices=directory(f"{OUT}/matrices/{{case}}/{{kind}}")
    params:
        implementation=implementation("matrices", "truth"),
        cfg=MATRIX_CFG, step="derive", scratch=lambda wc: f"{OUT}/work/{wc.case}/scratch"
    threads: 1
    resources: mem_mb=1024, disk_mb=200
    log: f"{OUT}/logs/{{case}}/derive/{{kind}}.log"
    benchmark: f"{OUT}/benchmarks/{{case}}/derive/{{kind}}.tsv"
    conda: PYENV
    script: STEP

rule truth_lanes:
    input:
        expectations=f"{OUT}/design/{{case}}/expectations.tsv",
        samples=f"{OUT}/design/{{case}}/samples.tsv",
        matrices=f"{OUT}/matrices/{{case}}/{{kind}}",
        matching=lambda wc: [f"{OUT}/catalogues/{wc.case}/matching.tsv"] if wc.kind == "mags" and not MATRIX_MODE else [],
        policies=POLICY_FILE
    output:
        lanes=f"{OUT}/truth/{{case}}/{{kind}}.tsv"
    params:
        implementation=implementation("truth"),
        cfg=TRUTH_CFG, step="truth", scratch=lambda wc: f"{OUT}/work/{wc.case}/scratch"
    threads: 1
    resources: mem_mb=1024, disk_mb=100
    log: f"{OUT}/logs/{{case}}/truth/{{kind}}.log"
    benchmark: f"{OUT}/benchmarks/{{case}}/truth/{{kind}}.tsv"
    conda: PYENV
    script: STEP

rule merge_truth_lanes:
    input:
        lanes=lambda wc: [f"{OUT}/truth/{wc.case}/{kind}.tsv" for kind in KINDS]
    output:
        lanes=f"{OUT}/truth/{{case}}/lanes.tsv"
    params:
        implementation=implementation("truth"),
        cfg=TRUTH_CFG, step="merge_truth", scratch=lambda wc: f"{OUT}/work/{wc.case}/scratch"
    threads: 1
    resources: mem_mb=512, disk_mb=50
    log: f"{OUT}/logs/{{case}}/truth/merge.log"
    benchmark: f"{OUT}/benchmarks/{{case}}/truth_merge.tsv"
    conda: PYENV
    script: STEP

rule differential_abundance:
    input:
        script=str(REPO / "workflow/scripts/run_da.R"),
        contract=str(ADAPTER_DIR / "contract.R"),
        adapter=lambda wc: str(ADAPTER_DIR / f"{REG.lookup(wc.method)['adapter']}.R"),
        registry=REGISTRY,
        policies=POLICY_FILE,
        matrices=f"{OUT}/matrices/{{case}}/{{kind}}",
        samples=f"{OUT}/design/{{case}}/samples.tsv"
    output:
        results=f"{OUT}/da/{{case}}/{{kind}}/{{metric}}/{{method}}.tsv",
        status=f"{OUT}/da/{{case}}/{{kind}}/{{metric}}/{{method}}.status.json",
        session=f"{OUT}/da/{{case}}/{{kind}}/{{metric}}/{{method}}.session.txt"
    params:
        matrix=lambda wc: f"{OUT}/matrices/{wc.case}/{wc.kind}/{wc.metric}.tsv",
        features=lambda wc: f"{OUT}/matrices/{wc.case}/{wc.kind}/features.tsv",
        settings=lambda wc: REG.entry_json(wc.method),
        alpha=CFG["analysis"]["alpha"],
        mc=CFG["analysis"]["mc_samples"],
        seed=lambda wc: int(wc.case.split("_s")[1]),
        timeout=int(CFG["analysis"]["da_timeout_minutes"]) * 60,
        allow_failures=str(CFG["analysis"]["allow_method_failures"]).lower(),
        tmp=lambda wc: f"{OUT}/work/{wc.case}/scratch/da_{wc.kind}_{wc.metric}_{wc.method}"
    threads: 1
    resources:
        mem_mb=lambda wc: int(REG.lookup(wc.method).get("mem_mb") or CFG["resources"]["da_mem_mb"]),
        disk_mb=1000,
        runtime=lambda wc: int(REG.lookup(wc.method).get("timeout_minutes")
                               or CFG["analysis"]["da_timeout_minutes"]) * 60
    log: f"{OUT}/logs/{{case}}/da/{{kind}}/{{metric}}/{{method}}.log"
    benchmark: f"{OUT}/benchmarks/{{case}}/da/{{kind}}/{{metric}}/{{method}}.tsv"
    conda: lambda wc: str(REG.execution_env(wc.method, ENV_PREFIX))
    shell:
        """
        mkdir -p {params.tmp:q}
        env TMPDIR={params.tmp:q} OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
            Rscript {input.script:q} {params.settings:q} {params.matrix:q} {input.samples:q} \
            {params.features:q} {output.results:q} {output.status:q} {output.session:q} \
            {params.alpha} {params.mc} {params.seed} {params.timeout} {params.allow_failures} \
            {wildcards.case:q} {wildcards.metric:q} {wildcards.kind:q} {params.tmp:q} > {log:q} 2>&1
        """

rule evaluate_case:
    input:
        truth=f"{OUT}/design/{{case}}/truth.tsv",
        lanes=f"{OUT}/truth/{{case}}/lanes.tsv",
        matching=f"{OUT}/catalogues/{{case}}/matching.tsv",
        results=lambda wc: [f"{OUT}/da/{wc.case}/{kind}/{metric}/{method}.tsv" for kind in KINDS for metric, method in COMBINATIONS],
        statuses=lambda wc: [f"{OUT}/da/{wc.case}/{kind}/{metric}/{method}.status.json" for kind in KINDS for metric, method in COMBINATIONS]
    output:
        scores=f"{OUT}/evaluation/{{case}}/scores.tsv",
        features=f"{OUT}/evaluation/{{case}}/feature_evaluation.tsv",
        recovery=f"{OUT}/evaluation/{{case}}/recovery.tsv",
        statuses=f"{OUT}/evaluation/{{case}}/method_statuses.json"
    params:
        implementation=implementation("evaluation"),
        cfg=CFG, step="evaluate", scratch=lambda wc: f"{OUT}/work/{wc.case}/scratch"
    threads: 1
    resources: mem_mb=1024, disk_mb=200
    log: f"{OUT}/logs/{{case}}/evaluate.log"
    benchmark: f"{OUT}/benchmarks/{{case}}/evaluate.tsv"
    conda: PYENV
    script: STEP