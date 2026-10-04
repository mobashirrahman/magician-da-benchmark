"""Generate independent compositions and truth for relative genomic abundance."""
from pathlib import Path
import math
import re
import numpy as np
import pandas as pd
from .config import resolve, case_info, samples
from .io import fasta, write_fasta, table, save_table, sha256, write_json


def prepare_references(cfg, output):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    manifest = cfg["genomes"]["manifest"]
    contigs, sources = [], []
    if manifest:
        df = table(resolve(manifest))
        if not {"genome_id", "path"} <= set(df.columns) or df.empty:
            raise ValueError("Genome manifest must contain genome_id and path columns")
        entries = [(str(r.genome_id), resolve(r.path)) for r in df.itertuples()]
    else:
        rng = np.random.default_rng(73021)
        entries = []
        n = cfg["genomes"]["synthetic_count"]
        for i in range(n):
            gid = f"synthetic_{i + 1:03d}"
            gc = 0.25 + 0.5 * i / max(n - 1, 1)
            seq = "".join(rng.choice(list("ACGT"), cfg["genomes"]["synthetic_length"],
                                     p=[(1-gc)/2, gc/2, gc/2, (1-gc)/2]))
            path = out / f"{gid}.fa"
            minimum = max(cfg["simulation"]["insert_size"] + 4 * cfg["simulation"]["insert_sd"],
                          cfg["assembly"]["min_contig_length"])
            chunks = min(cfg["genomes"]["synthetic_contigs"], max(1, len(seq) // minimum))
            boundaries = np.linspace(0, len(seq), chunks + 1, dtype=int)
            write_fasta([(f"{gid}_{j+1}", seq[boundaries[j]:boundaries[j+1]]) for j in range(chunks)], path)
            entries.append((gid, path))
    ids = [gid for gid, _ in entries]
    if len(ids) != len(set(ids)) or any(not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*", gid) for gid in ids):
        raise ValueError("Genome IDs must be unique and contain only letters, digits, _ or -")
    if len(ids) < 4:
        raise ValueError("At least four source genomes are required")
    threshold = cfg["simulation"]["insert_size"] + 4 * cfg["simulation"]["insert_sd"]
    for gid, path in entries:
        length = 0
        normalized = []
        for j, (_, seq) in enumerate(fasta(path)):
            if len(seq) < threshold or set(seq) - set("ACGTN") or seq.count("N") / len(seq) > 0.05:
                raise ValueError(f"{gid}: reference contigs must be >= {threshold} bp with <=5% N")
            cid = f"{gid}__c{j + 1:05d}"
            normalized.append((cid, seq))
            contigs.append(dict(contig_id=cid, genome_id=gid, length_bp=len(seq)))
            length += len(seq)
        if not normalized:
            raise ValueError(f"Empty reference genome: {path}")
        from .execution import check_storage
        check_storage(cfg, required=length * 2)
        target = out / f"{gid}.fa"
        write_fasta(normalized, target)
        with target.open() as source, (out / "combined.fa").open("a" if sources else "w") as combined:
            import shutil
            shutil.copyfileobj(source, combined, length=1024 * 1024)
        sources.append(dict(genome_id=gid, path=str(path.resolve()), length_bp=length,
                            sha256=sha256(path), synthetic=not bool(manifest)))
    save_table(pd.DataFrame(sources), out / "genomes.tsv")
    save_table(pd.DataFrame(contigs), out / "contigs.tsv")


def allocate(total, probabilities):
    expected = np.asarray(probabilities) * total
    result = np.floor(expected).astype(int)
    order = np.argsort(-(expected - result), kind="stable")
    result[order[:total - int(result.sum())]] += 1
    assert result.sum() == total
    return result


def make_design(cfg, case, references, output):
    scenario, seed = case_info(case)
    exp = cfg["experiment"]
    # Same baseline, noise, and library sizes in null/spiked runs for paired comparison.
    rng = np.random.default_rng(seed)
    ref = table(Path(references) / "genomes.tsv")
    ids, lengths = ref.genome_id.tolist(), ref.length_bp.to_numpy(dtype=float)
    baseline = rng.lognormal(0, 0.5, len(ids))
    baseline /= baseline.sum()
    treated = baseline.copy()
    n_da = max(2, int(round(len(ids) * exp["differential_fraction"])))
    selected = rng.choice(len(ids), min(len(ids), n_da), replace=False)
    if scenario == "spiked":
        treated[selected[:len(selected)//2]] *= exp["effect_size"]
        treated[selected[len(selected)//2:]] /= exp["effect_size"]
        treated[selected] *= baseline[selected].sum() / treated[selected].sum()
    effect = np.log2(treated / baseline)
    truth = pd.DataFrame({"feature_id": ids, "is_da": np.abs(effect) > 1e-10,
                          "true_log2fc": effect, "baseline_proportion": baseline,
                          "treatment_proportion": treated, "length_bp": lengths.astype(int)})
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    metadata, allocations, expectations = [], [], []
    bsigma = math.sqrt(math.log1p(exp["biological_cv"]**2))
    lsigma = math.sqrt(math.log1p(exp["library_cv"]**2))
    # Total genome copies per sample. Only simulated when load information is
    # requested; without it the absolute abundance lane stays explicitly unavailable.
    lcv = math.sqrt(math.log1p(exp["load_cv"]**2)) if exp["load_information"] else 0.0
    for i, sample in enumerate(samples(cfg)):
        group = "Control" if sample.startswith("C") else "Treatment"
        target = baseline if group == "Control" else treated
        abundance = target * rng.lognormal(-bsigma**2 / 2, bsigma, len(ids))
        abundance /= abundance.sum()
        pair_count = max(1, round(exp["read_pairs"] * rng.lognormal(-lsigma**2 / 2, lsigma)))
        dna = abundance * lengths
        dna /= dna.sum()
        counts = allocate(pair_count, dna)
        load = float(rng.lognormal(-lcv**2 / 2, lcv)) if lcv else float("nan")
        metadata.append(dict(sample_id=sample, group=group, read_pairs=pair_count,
                             simulation_seed=(seed + 1009 * (i + 1)) % (2**31 - 1) or 1))
        for gid, intended, p, q, count in zip(ids, target, abundance, dna, counts):
            allocations.append(dict(sample_id=sample, genome_id=gid,
                                    genomic_proportion=p, read_proportion=q, read_pairs=int(count)))
            # target_proportion is the group composition the design intends; p is the
            # noisy draw around it. Truth labels come from the target, truth effect sizes
            # from the expectation of the measured statistic.
            expectations.append(dict(sample_id=sample, genome_id=gid,
                                     genomic_proportion=p, target_proportion=float(intended),
                                     read_proportion=q, library_read_pairs=pair_count,
                                     length_bp=float(lengths[ids.index(gid)]),
                                     expected_count=count, absolute_copies=p * load))
    save_table(truth, out / "truth.tsv")
    save_table(pd.DataFrame(metadata), out / "samples.tsv")
    save_table(pd.DataFrame(allocations), out / "allocations.tsv")
    # Latent expectations for every truth lane, computed from the generator rather
    # than from the observed counts, so truth never depends on method filtering.
    save_table(pd.DataFrame(expectations), out / "expectations.tsv")
    write_json(out / "design.json", {"case": case, "scenario": scenario, "seed": seed,
                                    "estimand": "relative genomic abundance",
                                    "truth_lanes": ["relative_genomic_abundance", "expected_read_fraction",
                                                    "clr_log_ratio", "reference_relative_change",
                                                    "absolute_abundance"],
                                    "load_information": bool(exp["load_information"]),
                                    "independent_samples": True,
                                    "n_source_genomes": len(ids),
                                    "actual_read_pairs": sum(x["read_pairs"] for x in metadata),
                                    "config": cfg})
