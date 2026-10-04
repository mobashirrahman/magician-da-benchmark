rule simulate_reads:
    input:
        refs=REF,
        allocations=f"{OUT}/design/{{case}}/allocations.tsv",
        samples=f"{OUT}/design/{{case}}/samples.tsv",
        preflight=f"{OUT}/provenance/preflight.json",
        previous=lambda wc: [f"{OUT}/matrices/{PREVIOUS[wc.case]}/{kind}" for kind in KINDS] if PREVIOUS[wc.case] else []
    output:
        r1=temp(f"{OUT}/work/{{case}}/reads/{{sample}}.R1.fq.gz"),
        r2=temp(f"{OUT}/work/{{case}}/reads/{{sample}}.R2.fq.gz")
    params:
        implementation=implementation("reads"),
        cfg=cfg_for("simulation"), step="simulate", scratch=lambda wc: f"{OUT}/work/{wc.case}/scratch"
    threads: 1
    resources: mem_mb=512, disk_mb=lambda wc: work_disk_mb(wc, 2), io_jobs=1
    log: f"{OUT}/logs/{{case}}/simulate/{{sample}}.log"
    benchmark: f"{OUT}/benchmarks/{{case}}/simulate/{{sample}}.tsv"
    conda: BIOENV
    script: STEP

rule assemble:
    input:
        r1=lambda wc: expand(f"{OUT}/work/{{case}}/reads/{{sample}}.R1.fq.gz", case=wc.case, sample=SAMPLES),
        r2=lambda wc: expand(f"{OUT}/work/{{case}}/reads/{{sample}}.R2.fq.gz", case=wc.case, sample=SAMPLES)
    output: contigs=temp(f"{OUT}/work/{{case}}/assembly/contigs.fa")
    params:
        implementation=implementation("reads"),
        cfg=cfg_for("assembly"), step="assemble", scratch=lambda wc: f"{OUT}/work/{wc.case}/scratch"
    threads: CFG["resources"]["threads"]
    resources: mem_mb=CFG["resources"]["assembly_mem_mb"], disk_mb=lambda wc: work_disk_mb(wc, 6), io_jobs=1
    log: f"{OUT}/logs/{{case}}/assemble.log"
    benchmark: f"{OUT}/benchmarks/{{case}}/assemble.tsv"
    conda: BIOENV
    script: STEP

rule map_contigs:
    input:
        ref=f"{OUT}/work/{{case}}/assembly/contigs.fa",
        r1=f"{OUT}/work/{{case}}/reads/{{sample}}.R1.fq.gz",
        r2=f"{OUT}/work/{{case}}/reads/{{sample}}.R2.fq.gz"
    output:
        bam=temp(f"{OUT}/work/{{case}}/alignments/{{sample}}.bam"),
        counts=temp(f"{OUT}/work/{{case}}/counts/{{sample}}.contigs.tsv")
    params:
        implementation=implementation("reads"),
        cfg=cfg_for("mapping"), step="map_contigs", scratch=lambda wc: f"{OUT}/work/{wc.case}/scratch"
    threads: CFG["resources"]["threads"]
    resources: mem_mb=CFG["resources"]["mapping_mem_mb"], disk_mb=lambda wc: work_disk_mb(wc, 3), io_jobs=1
    log: f"{OUT}/logs/{{case}}/map_contigs/{{sample}}.log"
    benchmark: f"{OUT}/benchmarks/{{case}}/map_contigs/{{sample}}.tsv"
    conda: BIOENV
    script: STEP

rule map_sources:
    input:
        ref=REF,
        r1=f"{OUT}/work/{{case}}/reads/{{sample}}.R1.fq.gz",
        r2=f"{OUT}/work/{{case}}/reads/{{sample}}.R2.fq.gz"
    output: counts=temp(f"{OUT}/work/{{case}}/counts/{{sample}}.sources.tsv")
    params:
        implementation=implementation("reads"),
        cfg=cfg_for("mapping"), step="map_sources", scratch=lambda wc: f"{OUT}/work/{wc.case}/scratch"
    threads: CFG["resources"]["threads"]
    resources: mem_mb=CFG["resources"]["mapping_mem_mb"], disk_mb=1000, io_jobs=1
    log: f"{OUT}/logs/{{case}}/map_sources/{{sample}}.log"
    benchmark: f"{OUT}/benchmarks/{{case}}/map_sources/{{sample}}.tsv"
    conda: BIOENV
    script: STEP

rule contig_depth:
    input:
        contigs=f"{OUT}/work/{{case}}/assembly/contigs.fa",
        bams=lambda wc: expand(f"{OUT}/work/{{case}}/alignments/{{sample}}.bam", case=wc.case, sample=SAMPLES)
    output: depth=temp(f"{OUT}/work/{{case}}/depth.tsv")
    params:
        implementation=implementation("reads"),
        cfg=cfg_for(), step="depth", scratch=lambda wc: f"{OUT}/work/{wc.case}/scratch"
    threads: 1
    resources: mem_mb=1024, disk_mb=1000, io_jobs=1
    log: f"{OUT}/logs/{{case}}/contig_depth.log"
    benchmark: f"{OUT}/benchmarks/{{case}}/contig_depth.tsv"
    conda: BIOENV
    script: STEP

rule bin_catalogue:
    input:
        contigs=f"{OUT}/work/{{case}}/assembly/contigs.fa",
        depth=f"{OUT}/work/{{case}}/depth.tsv"
    output: catalogue=temp(directory(f"{OUT}/work/{{case}}/catalogue"))
    params:
        implementation=implementation("catalog"),
        cfg=cfg_for("binning"), step="bin", scratch=lambda wc: f"{OUT}/work/{wc.case}/scratch"
    threads: CFG["resources"]["threads"]
    resources: mem_mb=CFG["resources"]["binning_mem_mb"], disk_mb=1000
    log: f"{OUT}/logs/{{case}}/bin_catalogue.log"
    benchmark: f"{OUT}/benchmarks/{{case}}/bin_catalogue.tsv"
    conda: BIOENV
    script: STEP

rule align_mags:
    input:
        refs=REF,
        catalogue=f"{OUT}/work/{{case}}/catalogue"
    output: paf=temp(f"{OUT}/work/{{case}}/matches.paf")
    params:
        implementation=implementation("catalog"),
        cfg=cfg_for(), step="align_mags", scratch=lambda wc: f"{OUT}/work/{wc.case}/scratch"
    threads: CFG["resources"]["threads"]
    resources: mem_mb=CFG["resources"]["mapping_mem_mb"], disk_mb=1000
    log: f"{OUT}/logs/{{case}}/align_mags.log"
    benchmark: f"{OUT}/benchmarks/{{case}}/align_mags.tsv"
    conda: BIOENV
    script: STEP

rule match_mags:
    input:
        paf=f"{OUT}/work/{{case}}/matches.paf",
        catalogue=f"{OUT}/work/{{case}}/catalogue",
        references=f"{OUT}/provenance/references.tsv"
    output:
        matching=f"{OUT}/catalogues/{{case}}/matching.tsv",
        evidence=f"{OUT}/catalogues/{{case}}/match_evidence.tsv",
        features=f"{OUT}/catalogues/{{case}}/mags.tsv",
        assignments=f"{OUT}/catalogues/{{case}}/assignments.tsv",
        status=f"{OUT}/catalogues/{{case}}/catalogue.json"
    params:
        implementation=implementation("catalog"),
        cfg=cfg_for("matching"), step="match_mags", scratch=lambda wc: f"{OUT}/work/{wc.case}/scratch"
    threads: 1
    resources: mem_mb=1024, disk_mb=100
    log: f"{OUT}/logs/{{case}}/match_mags.log"
    benchmark: f"{OUT}/benchmarks/{{case}}/match_mags.tsv"
    conda: PYENV
    script: STEP
