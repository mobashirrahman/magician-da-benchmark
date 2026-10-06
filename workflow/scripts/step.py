"""Snakemake adapter: all scientific logic is importable and independently testable."""
from pathlib import Path
import os
import sys
import time
import shutil
import traceback

REPO = Path(__file__).resolve().parents[2]
# Snakemake relocates script copies; params contains original config and source path.
if not (REPO / "src/magician").exists():
    REPO = Path(snakemake.scriptdir).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
from magician import design, reads, catalog, matrices, evaluation, report, truth, matrix
from magician.config import case_info, zero_policies
from magician.io import write_json, read_json
from magician.execution import check_storage

cfg = snakemake.params.cfg
step = snakemake.params.step
scratch = Path(snakemake.params.scratch)
scratch.mkdir(parents=True, exist_ok=True)
os.environ["TMPDIR"] = str(scratch)
Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
with open(snakemake.log[0], "w", buffering=1) as log:
    os.dup2(log.fileno(), 1)
    os.dup2(log.fileno(), 2)
    print(f"START {step} {time.strftime('%Y-%m-%dT%H:%M:%S%z')}", flush=True)
    try:
        check_storage(cfg)
        i, o = snakemake.input, snakemake.output
        if step == "prepare_references":
            design.prepare_references(cfg, o.refs)
            Path(o.provenance).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(Path(o.refs) / "genomes.tsv", o.provenance)
        elif step == "design":
            design.make_design(cfg, snakemake.wildcards.case, i.refs, Path(o.design).parent)
        elif step == "preflight":
            report.preflight(cfg, i.refs, i.designs, o.preflight)
            write_json(o.config, cfg)
        elif step == "simulate":
            reads.simulate(cfg, snakemake.wildcards.sample, Path(i.allocations).parent, i.refs, o.r1, o.r2, scratch)
        elif step == "assemble":
            reads.assemble(cfg, i.r1, i.r2, o.contigs, scratch, snakemake.threads, snakemake.resources.mem_mb)
        elif step == "map_contigs":
            reads.map_contigs(cfg, i.ref, i.r1, i.r2, o.bam, o.counts, scratch, snakemake.threads, snakemake.resources.mem_mb)
        elif step == "map_sources":
            reads.map_sources(cfg, Path(i.ref) / "combined.fa", i.r1, i.r2, o.counts, snakemake.threads)
        elif step == "depth":
            reads.depth(cfg, i.bams, i.contigs, o.depth)
        elif step == "bin":
            catalog.bin_catalogue(cfg, i.contigs, i.depth, o.catalogue, scratch, snakemake.threads, case_info(snakemake.wildcards.case)[1])
        elif step == "align_mags":
            catalog.align_mags(cfg, Path(i.refs) / "combined.fa", Path(i.catalogue) / "combined.fa", o.paf, snakemake.threads)
        elif step == "match_mags":
            catalog.match_paf(cfg, i.paf, Path(i.catalogue) / "mags.tsv", i.references, o.matching)
            for name in ("mags.tsv", "assignments.tsv", "catalogue.json"):
                shutil.copyfile(Path(i.catalogue) / name, Path(o.matching).parent / name)
            # Tier 2 ablation: completeness/contamination filter (none/medium/strict).
            kept, filt_report = catalog.apply_completeness_filter(
                cfg, Path(o.matching).parent / "mags.tsv", o.matching)
            if filt_report.get("n_removed", 0):
                import pandas as pd
                from magician.io import table as _table, save_table as _save
                mags = _table(Path(o.matching).parent / "mags.tsv")
                mags = mags.loc[mags.feature_id.isin(kept)]
                _save(mags, Path(o.matching).parent / "mags.tsv")
                assign = _table(Path(o.matching).parent / "assignments.tsv")
                assign = assign.loc[assign.feature_id.isin(kept)]
                _save(assign, Path(o.matching).parent / "assignments.tsv")
                matching = _table(o.matching)
                matching = matching.loc[matching.feature_id.isin(kept)]
                _save(matching, o.matching)
            write_json(Path(o.matching).parent / "completeness_filter.json", filt_report)
        elif step in {"matrices", "derive"}:
            matrices.derive(cfg, snakemake.wildcards.kind, i.raw, o.matrices)
        elif step == "quantify":
            matrices.quantify(cfg, snakemake.wildcards.kind, i.counts, i.samples, i.features,
                              i.assignments[0] if i.assignments else None, o.raw)
        elif step == "truth":
            analysis = cfg["analysis"]
            matching = Path(i.matching[0]) if hasattr(i, "matching") and len(i.matching) else None
            truth.build(cfg, snakemake.wildcards.case, snakemake.wildcards.kind, i.expectations,
                        i.samples, Path(i.matrices) / "features.tsv", o.lanes,
                        zero_policies=zero_policies(cfg), reference=analysis["truth_reference"],
                        reference_feature=analysis["truth_reference_feature"],
                        matching_file=matching)
        elif step == "merge_truth":
            import pandas as pd
            from magician.io import save_table
            save_table(pd.concat([pd.read_csv(path, sep="\t", keep_default_na=False,
                                              na_values=["NA"]) for path in i.lanes], ignore_index=True), o.lanes)
        elif step == "import_matrix":
            matrix.import_case(cfg, snakemake.wildcards.case, Path(o.raw).parents[2])
        elif step == "import_manifest":
            matrix.import_manifest(cfg, o.manifest)
        elif step == "select_reference":
            from magician import reference
            reference.select_present(i.refs, cfg["reference_table"]["fraction_present"],
                                     cfg["reference_table"]["seed"], o.genomes)
        elif step == "filter_reference_counts":
            from magician import reference
            reference.filter_counts(i.source, i.genomes, o.counts)
        elif step == "reference_catalogue":
            from magician import reference
            reference.build_catalogue(i.refs, i.genomes, o.features, o.matching)
        elif step == "evaluate":
            ref_match = None
            try:
                cand = Path(i.matching).parent / "reference_matching.tsv"
                if cand.exists():
                    ref_match = str(cand)
            except Exception:
                ref_match = None
            evaluation.evaluate_case(cfg, snakemake.wildcards.case, i.truth, i.lanes, i.matching,
                                     i.results, i.statuses, Path(o.scores).parent,
                                     reference_matching_file=ref_match)
        elif step == "aggregate":
            from magician.config import registry
            reg = registry(cfg)
            evaluation.aggregate(cfg, i.scores, i.recovery, Path(o.scores).parent,
                                 settings_rows=reg.settings_table(cfg["analysis"]["methods"],
                                                                  cfg["analysis"]["metrics"]))
        elif step == "report":
            report.make_report(cfg, Path(i.ranking).parent, i.benchmarks, o.report, o.storage)
        else:
            raise ValueError(f"Unknown step: {step}")
        print(f"SUCCESS {step} {time.strftime('%Y-%m-%dT%H:%M:%S%z')}", flush=True)
    except BaseException:
        traceback.print_exc()
        raise
