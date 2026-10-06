"""Tier 2 third table: deliberately incomplete reference (70% present).

Only included when reference_table.enabled is true. The present set is fixed
per run; per-sample reference counts filter source counts to present genomes.
"""

rule select_reference:
    input:
        refs=f"{OUT}/provenance/references.tsv"
    output:
        genomes=f"{OUT}/provenance/reference_genomes.tsv"
    params:
        implementation=implementation("reference"),
        cfg=cfg_for("reference_table"), step="select_reference",
        scratch=f"{OUT}/work/scratch/setup"
    threads: 1
    resources: mem_mb=512, disk_mb=100
    log: f"{OUT}/logs/select_reference.log"
    benchmark: f"{OUT}/benchmarks/select_reference.tsv"
    conda: PYENV
    script: STEP


rule make_reference_counts:
    input:
        source=f"{OUT}/work/{{case}}/counts/{{sample}}.sources.tsv",
        genomes=f"{OUT}/provenance/reference_genomes.tsv"
    output:
        counts=temp(f"{OUT}/work/{{case}}/counts/{{sample}}.reference.tsv")
    params:
        implementation=implementation("reference"),
        cfg=cfg_for("reference_table"), step="filter_reference_counts",
        scratch=lambda wc: f"{OUT}/work/{wc.case}/scratch"
    threads: 1
    resources: mem_mb=512, disk_mb=500
    log: f"{OUT}/logs/{{case}}/reference_counts/{{sample}}.log"
    benchmark: f"{OUT}/benchmarks/{{case}}/reference_counts/{{sample}}.tsv"
    conda: PYENV
    script: STEP


rule make_reference_catalogue:
    input:
        refs=f"{OUT}/provenance/references.tsv",
        genomes=f"{OUT}/provenance/reference_genomes.tsv"
    output:
        features=f"{OUT}/catalogues/{{case}}/reference.tsv",
        matching=f"{OUT}/catalogues/{{case}}/reference_matching.tsv"
    params:
        implementation=implementation("reference"),
        cfg=cfg_for("reference_table"), step="reference_catalogue",
        scratch=lambda wc: f"{OUT}/work/{wc.case}/scratch"
    threads: 1
    resources: mem_mb=512, disk_mb=100
    log: f"{OUT}/logs/{{case}}/reference_catalogue.log"
    benchmark: f"{OUT}/benchmarks/{{case}}/reference_catalogue.tsv"
    conda: PYENV
    script: STEP
