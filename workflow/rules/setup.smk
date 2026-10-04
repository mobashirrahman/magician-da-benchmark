rule prepare_references:
    input:
        manifest=lambda wc: [str(REPO / CFG["genomes"]["manifest"])] if CFG["genomes"]["manifest"] else [],
        sources=lambda wc: __import__("magician.io", fromlist=["table"]).table(REPO / CFG["genomes"]["manifest"]).path.map(lambda p: str(REPO / p)).tolist() if CFG["genomes"]["manifest"] else []
    output:
        refs=temp(directory(REF)),
        provenance=f"{OUT}/provenance/references.tsv"
    params:
        implementation=implementation("design"),
        cfg=cfg_for("genomes", "simulation", "assembly"), step="prepare_references", scratch=f"{OUT}/work/scratch/setup"
    threads: 1
    resources: mem_mb=1024, disk_mb=1024
    log: f"{OUT}/logs/prepare_references.log"
    benchmark: f"{OUT}/benchmarks/prepare_references.tsv"
    conda: PYENV
    script: STEP

rule design:
    input:
        refs=REF,
        provenance=f"{OUT}/provenance/references.tsv"
    output:
        truth=f"{OUT}/design/{{case}}/truth.tsv",
        samples=f"{OUT}/design/{{case}}/samples.tsv",
        allocations=f"{OUT}/design/{{case}}/allocations.tsv",
        expectations=f"{OUT}/design/{{case}}/expectations.tsv",
        design=f"{OUT}/design/{{case}}/design.json"
    params:
        implementation=implementation("design"),
        cfg=cfg_for("experiment"), step="design", scratch=f"{OUT}/work/scratch/setup"
    threads: 1
    resources: mem_mb=512, disk_mb=100
    log: f"{OUT}/logs/design/{{case}}.log"
    benchmark: f"{OUT}/benchmarks/design/{{case}}.tsv"
    conda: PYENV
    script: STEP

rule preflight:
    input:
        refs=f"{OUT}/provenance/references.tsv",
        designs=expand(f"{OUT}/design/{{case}}/design.json", case=CASES)
    output:
        preflight=f"{OUT}/provenance/preflight.json",
        config=f"{OUT}/provenance/preflight_config.json"
    params:
        implementation=implementation("storage"),
        cfg=cfg_for("genomes", "experiment", "simulation"), step="preflight", scratch=f"{OUT}/work/scratch/setup"
    threads: 1
    resources: mem_mb=512, disk_mb=100
    log: f"{OUT}/logs/preflight.log"
    benchmark: f"{OUT}/benchmarks/preflight.tsv"
    conda: PYENV
    script: STEP