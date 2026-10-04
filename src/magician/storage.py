"""Storage preflight independent of plotting and report formatting."""
from pathlib import Path
import platform
import pandas as pd
from .io import read_json, write_json


def preflight(cfg, references, design_files, output):
    # Retained lengths suffice; requesting temporary FASTAs here would rebuild
    # reads/assemblies when only the analysis configuration changes.
    genomes = pd.read_csv(references, sep="\t")
    designs = [read_json(path) for path in design_files]
    pairs = max(d["actual_read_pairs"] for d in designs)
    # Deliberately conservative: raw FASTQ plus BAM, sorting and assembler scratch.
    read_bytes = pairs * (4 * cfg["simulation"]["read_length"] + 300)
    peak = read_bytes * 8 + int(genomes.length_bp.sum()) * 20
    retained = len(designs) * max(len(genomes), 1) * cfg["experiment"]["samples_per_group"] * 100000
    reserve = (cfg["storage"]["environment_reserve_gb"] + cfg["storage"]["headroom_gb"]) * 10**9
    estimate = dict(storage_budget_bytes=int(cfg["storage"]["budget_gb"] * 10**9),
                    estimated_peak_work_bytes=int(peak), estimated_retained_bytes=int(retained),
                    environment_and_headroom_bytes=int(reserve), cases_run_sequentially=True,
                    estimated_peak_total_bytes=int(peak + retained + reserve),
                    synthetic_fixture=cfg["genomes"]["manifest"] is None, config=cfg,
                    python=platform.python_version())
    write_json(output, estimate)
    if estimate["estimated_peak_total_bytes"] > estimate["storage_budget_bytes"]:
        raise ValueError("Estimated storage exceeds budget. Reduce read pairs, samples or genomes; see preflight.json")
    from .execution import check_storage
    check_storage(cfg, required=int(peak))
