rule aggregate_results:
    input:
        scores=expand(f"{OUT}/evaluation/{{case}}/scores.tsv", case=CASES),
        recovery=expand(f"{OUT}/evaluation/{{case}}/recovery.tsv", case=CASES)
    output:
        scores=f"{OUT}/benchmark/scores.tsv",
        recovery=f"{OUT}/benchmark/recovery.tsv",
        ranking=f"{OUT}/benchmark/ranking.tsv",
        settings=f"{OUT}/benchmark/method_settings.tsv",
        recommendations=f"{OUT}/benchmark/recommendations.json"
    params:
        implementation=implementation("evaluation"),
        cfg=CFG, step="aggregate", scratch=f"{OUT}/work/scratch/report"
    threads: 1
    resources: mem_mb=1024, disk_mb=100
    log: f"{OUT}/logs/aggregate.log"
    benchmark: f"{OUT}/benchmarks/aggregate.tsv"
    conda: PYENV
    script: STEP

rule benchmark_report:
    input:
        ranking=f"{OUT}/benchmark/ranking.tsv",
        recommendations=f"{OUT}/benchmark/recommendations.json",
        settings=f"{OUT}/benchmark/method_settings.tsv",
        benchmarks=all_benchmarks()
    output:
        report=f"{OUT}/benchmark/report.html",
        storage=f"{OUT}/benchmark/storage.json",
        resources=f"{OUT}/benchmark/rule_resources.tsv",
        method_resources=f"{OUT}/benchmark/method_resources.tsv",
        config=f"{OUT}/provenance/config.json",
        svg=f"{OUT}/benchmark/performance.svg",
        png=f"{OUT}/benchmark/performance.png"
    params:
        implementation=implementation("report"),
        cfg=CFG, step="report", scratch=f"{OUT}/work/scratch/report"
    threads: 1
    resources: mem_mb=1024, disk_mb=100
    log: f"{OUT}/logs/report.log"
    benchmark: f"{OUT}/benchmarks/report.tsv"
    conda: PYENV
    script: STEP
