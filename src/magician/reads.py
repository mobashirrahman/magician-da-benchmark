"""Simulate only requested pairs and stream mapping outputs; never retain SAM."""
from pathlib import Path
import gzip
import shutil
import subprocess
import tempfile
import pandas as pd
from .execution import run, version
from .io import table, save_table, fasta, write_fasta
from .design import allocate


def simulate(cfg, sample, design, references, r1, r2, scratch):
    print(version(["wgsim"]))
    allocations = table(Path(design) / "allocations.tsv")
    alloc = allocations.loc[allocations.sample_id == sample]
    metadata = table(Path(design) / "samples.tsv").set_index("sample_id").loc[sample]
    sim = cfg["simulation"]
    Path(scratch).mkdir(parents=True, exist_ok=True)
    for p in (r1, r2):
        Path(p).parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=scratch, prefix=f"sim_{sample}_") as tmp:
        paths = [Path(tmp) / "R1.fq", Path(tmp) / "R2.fq"]
        handles = [gzip.open(r1, "wb", compresslevel=1), gzip.open(r2, "wb", compresslevel=1)]
        try:
            for i, row in enumerate(alloc.itertuples()):
                if not row.read_pairs:
                    continue
                sequences = list(fasta(Path(references) / f"{row.genome_id}.fa"))
                lengths = [len(seq) for _, seq in sequences]
                per_contig = allocate(int(row.read_pairs), [n / sum(lengths) for n in lengths])
                for j, (record, n_pairs) in enumerate(zip(sequences, per_contig)):
                    if not n_pairs:
                        continue
                    # wgsim rounds independently per input contig. One-contig calls
                    # guarantee exact allocations even for fragmented real genomes.
                    single = Path(tmp) / "source.fa"
                    write_fasta([record], single)
                    seed = (int(metadata.simulation_seed) + 7919 * i + 104729 * j) % (2**31 - 1) or 1
                    command = ["wgsim", "-N", str(n_pairs), "-1", str(sim["read_length"]),
                               "-2", str(sim["read_length"]), "-d", str(sim["insert_size"]),
                               "-s", str(sim["insert_sd"]), "-e", str(sim["error_rate"]),
                               "-r", "0", "-S", str(seed), str(single), *map(str, paths)]
                    run(command, cfg)
                    for path, handle in zip(paths, handles):
                        lines = 0
                        with path.open("rb") as source:
                            for chunk in iter(lambda: source.read(1024 * 1024), b""):
                                lines += chunk.count(b"\n")
                                handle.write(chunk)
                        if lines != int(n_pairs) * 4:
                            raise ValueError("Simulator produced an unexpected number of FASTQ records")
                        path.unlink()
        finally:
            for handle in handles:
                handle.close()


def assemble(cfg, r1s, r2s, output, scratch, threads, mem_mb):
    print(version(["megahit", "--version"]))
    Path(scratch).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=scratch, prefix="assembly_") as tmp:
        dest = Path(tmp) / "megahit"
        run(["megahit", "-1", ",".join(map(str, r1s)), "-2", ",".join(map(str, r2s)),
             "-o", str(dest), "-t", str(threads), "-m", str(int(mem_mb * 1024**2 * 0.85)),
             "--k-list", ",".join(map(str, cfg["assembly"]["k_list"])),
             "--min-contig-len", str(cfg["assembly"]["min_contig_length"])], cfg)
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(dest / "final.contigs.fa", output)


def count_sam(lines, feature_by_contig, mapq):
    counts = {feature: 0 for feature in dict.fromkeys(feature_by_contig.values())}
    for line in lines:
        if line.startswith("@"):
            continue
        fields = line.rstrip().split("\t")
        flag = int(fields[1])
        # Count one properly paired primary first mate, rather than both reads.
        if flag & 67 != 67 or flag & (4 | 256 | 2048) or int(fields[4]) < mapq:
            continue
        feature = feature_by_contig.get(fields[2])
        if feature is None:
            raise ValueError(f"Alignment references unknown contig {fields[2]}")
        counts[feature] += 1
    return counts


def map_contigs(cfg, ref, r1, r2, bam, counts, scratch, threads, mem_mb):
    print(version(["minimap2", "--version"]))
    print(version(["samtools", "--version"]))
    Path(bam).parent.mkdir(parents=True, exist_ok=True)
    contigs = {cid: len(seq) for cid, seq in fasta(ref)}
    if contigs:
        with Path(bam).open("wb") as handle:
            run([["minimap2", "-ax", "sr", "--secondary=no", "-t", str(max(1, threads - 1)),
                  str(ref), str(r1), str(r2)],
                 ["samtools", "sort", "-@", "0", "-m", f"{max(128, mem_mb//4)}M",
                  "-T", str(Path(scratch) / f"sort_{Path(bam).stem}"), "-O", "BAM", "-"]], cfg, handle)
        proc = subprocess.Popen(["samtools", "view", str(bam)], stdout=subprocess.PIPE, text=True)
        try:
            values = count_sam(proc.stdout, {c: c for c in contigs}, cfg["mapping"]["min_mapq"])
        finally:
            proc.stdout.close()
            code = proc.wait()
        if code:
            raise subprocess.CalledProcessError(code, "samtools view")
    else:
        with tempfile.TemporaryFile(mode="w+") as header:
            header.write("@HD\tVN:1.6\tSO:coordinate\n")
            header.seek(0)
            with Path(bam).open("wb") as handle:
                subprocess.run(["samtools", "view", "-b", "-"], stdin=header, stdout=handle, check=True)
        values = {}
    save_table(pd.DataFrame([dict(feature_id=k, count=v, length_bp=contigs[k]) for k, v in values.items()],
                           columns=["feature_id", "count", "length_bp"]), counts)


def map_sources(cfg, ref, r1, r2, counts, threads):
    refs = table(Path(ref).parent / "contigs.tsv")
    lengths = refs.groupby("genome_id").length_bp.sum().to_dict()
    mapping = dict(zip(refs.contig_id, refs.genome_id))
    print(version(["minimap2", "--version"]))
    # Stream SAM to the counter. It is never written to disk.
    proc = subprocess.Popen(["minimap2", "-ax", "sr", "--secondary=no", "-t", str(threads),
                             str(ref), str(r1), str(r2)], stdout=subprocess.PIPE, text=True)
    try:
        values = count_sam(proc.stdout, mapping, cfg["mapping"]["min_mapq"])
    except BaseException:
        proc.terminate()
        raise
    finally:
        proc.stdout.close()
        code = proc.wait()
    if code:
        raise subprocess.CalledProcessError(code, "minimap2 source mapping")
    save_table(pd.DataFrame([dict(feature_id=k, count=v, length_bp=lengths[k]) for k, v in values.items()]), counts)


def depth(cfg, bams, contigs, output):
    if not any(True for _ in fasta(contigs)):
        Path(output).write_text("contigName\tcontigLen\ttotalAvgDepth\n")
    else:
        print(version(["jgi_summarize_bam_contig_depths", "--help"]))
        run(["jgi_summarize_bam_contig_depths", "--outputDepth", str(output), *map(str, bams)], cfg)
