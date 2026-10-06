def implementation(*modules):
    import hashlib
    paths = [REPO / "src/magician" / f"{name}.py" for name in ("config", "io", "execution", *modules)]
    return hashlib.sha256(b"".join(path.read_bytes() for path in paths)).hexdigest()

def cfg_for(*sections):
    # DA method changes must not invalidate simulation/assembly/mapping jobs.
    return {key: CFG[key] for key in ("output_dir", "cache_dir", "storage", *sections)}

MATRIX_CFG = cfg_for()
MATRIX_CFG["analysis"] = {key: CFG["analysis"][key] for key in
                          ("metrics", "min_total_count", "min_prevalence", "reserve_other")}
TRUTH_CFG = cfg_for("analysis", "experiment", "genomes")
TRUTH_CFG["analysis"]["metrics"] = CFG["analysis"]["metrics"]

def estimated_pairs(wc):
    import math
    from magician.io import table
    path = Path(OUT) / "design" / wc.case / "samples.tsv"
    if path.exists():
        df = table(path)
        if hasattr(wc, "sample"):
            return int(df.set_index("sample_id").loc[wc.sample, "read_pairs"])
        return int(df.read_pairs.sum())
    # Dry-run fallback; real reservations use the generated allocation.
    exp = CFG["experiment"]
    upper = exp["read_pairs"] * math.exp(5 * math.sqrt(math.log1p(exp["library_cv"]**2)))
    return int(upper * (1 if hasattr(wc, "sample") else max(len(SAMPLES), 1)))

def work_disk_mb(wc, factor):
    return max(1000, int(estimated_pairs(wc) * (4 * CFG["simulation"]["read_length"] + 300) * factor / 1e6))

def all_benchmarks():
    files = [f"{OUT}/benchmarks/aggregate.tsv"]
    if MATRIX_MODE:
        files.append(f"{OUT}/benchmarks/import_matrix.tsv")
    else:
        files += [f"{OUT}/benchmarks/prepare_references.tsv", f"{OUT}/benchmarks/preflight.tsv"]
        if REF_TABLE:
            files.append(f"{OUT}/benchmarks/select_reference.tsv")
    for case in CASES:
        files += [f"{OUT}/benchmarks/design/{case}.tsv", f"{OUT}/benchmarks/{case}/evaluate.tsv"]
        if not MATRIX_MODE:
            files += [f"{OUT}/benchmarks/{case}/{rule}.tsv"
                      for rule in ["assemble", "contig_depth", "bin_catalogue", "align_mags", "match_mags"]]
            files += [f"{OUT}/benchmarks/{case}/{rule}/{sample}.tsv"
                      for rule in ["simulate", "map_contigs", "map_sources"] for sample in SAMPLES]
            if REF_TABLE:
                files += [f"{OUT}/benchmarks/{case}/reference_counts/{sample}.tsv" for sample in SAMPLES]
                files.append(f"{OUT}/benchmarks/{case}/reference_catalogue.tsv")
            files += [f"{OUT}/benchmarks/{case}/quantify/{kind}.tsv" for kind in KINDS]
        files += [f"{OUT}/benchmarks/{case}/derive/{kind}.tsv" for kind in KINDS]
        files += [f"{OUT}/benchmarks/{case}/truth/{kind}.tsv" for kind in KINDS]
        files.append(f"{OUT}/benchmarks/{case}/truth_merge.tsv")
        files += [f"{OUT}/benchmarks/{case}/da/{kind}/{metric}/{method}.tsv"
                  for kind in KINDS for metric, method in COMBINATIONS]
    return files