"""Shared MAG catalogue and explicit alignment-based source assignments."""
from pathlib import Path
from collections import defaultdict
import tempfile
import shutil
import pandas as pd
from .io import fasta, write_fasta, table, save_table, write_json
from .execution import run, version


def bin_catalogue(cfg, contigs, depth, output, scratch, threads, seed):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    rows, assignments = [], []
    records = []
    if any(True for _ in fasta(contigs)):
        print(version(["metabat2", "--help"]))
        with tempfile.TemporaryDirectory(dir=scratch, prefix="binning_") as tmp:
            prefix = str(Path(tmp) / "bin")
            run(["metabat2", "-i", str(contigs), "-a", str(depth), "-o", prefix,
                 "-t", str(threads), "-m", str(cfg["binning"]["min_contig_length"]),
                 "-s", str(cfg["binning"]["min_bin_size"]), "--seed", str(seed)], cfg)
            for i, path in enumerate(sorted(Path(tmp).glob("bin.*.fa"))):
                mag = f"MAG_{i + 1:05d}"
                total = 0
                normalized = []
                for j, (cid, seq) in enumerate(fasta(path)):
                    qid = f"{mag}__c{j + 1:05d}"
                    assignments.append(dict(contig_id=cid, feature_id=mag, length_bp=len(seq)))
                    normalized.append((qid, seq))
                    total += len(seq)
                if normalized:
                    write_fasta(normalized, out / f"{mag}.fa")
                    records.extend(normalized)
                    rows.append(dict(feature_id=mag, length_bp=total, n_contigs=len(normalized)))
    write_fasta(records, out / "combined.fa")
    save_table(pd.DataFrame(rows, columns=["feature_id", "length_bp", "n_contigs"]), out / "mags.tsv")
    save_table(pd.DataFrame(assignments, columns=["contig_id", "feature_id", "length_bp"]), out / "assignments.tsv")
    write_json(out / "catalogue.json", {"n_mags": len(rows), "status": "recovered" if rows else "no_mags"})


def union_length(intervals):
    total, end = 0, -1
    for start, stop in sorted(intervals):
        total += max(0, stop - max(start, end))
        end = max(end, stop)
    return total


def match_paf(cfg, paf, mags, references, output):
    metadata = table(mags).set_index("feature_id")
    refs = table(references).set_index("genome_id")
    query_intervals, target_intervals, identities = defaultdict(list), defaultdict(list), defaultdict(list)
    with Path(paf).open() as handle:
        for line in handle:
            f = line.rstrip().split("\t")
            mag = f[0].rsplit("__c", 1)[0]
            source = f[5].rsplit("__c", 1)[0]
            identity = int(f[9]) / max(int(f[10]), 1)
            if identity < cfg["matching"]["min_identity"]:
                continue
            if mag not in metadata.index or source not in refs.index:
                raise ValueError("PAF contains unknown MAG or reference ID")
            query_intervals[(mag, source, f[0])].append((int(f[2]), int(f[3])))
            target_intervals[(mag, source, f[5])].append((int(f[7]), int(f[8])))
            identities[(mag, source)].append((int(f[9]), int(f[10])))
    evidence = []
    for mag, source in identities:
        aligned = sum(union_length(v) for (m, s, _), v in query_intervals.items() if m == mag and s == source)
        recovered = sum(union_length(v) for (m, s, _), v in target_intervals.items() if m == mag and s == source)
        pairs = identities[(mag, source)]
        evidence.append(dict(feature_id=mag, source_id=source,
                             aligned_fraction=min(1, aligned / metadata.loc[mag, "length_bp"]),
                             source_completeness=min(1, recovered / refs.loc[source, "length_bp"]),
                             identity=sum(x[0] for x in pairs) / sum(x[1] for x in pairs)))
    rows = []
    for mag in metadata.index:
        candidates = sorted((e for e in evidence if e["feature_id"] == mag),
                            key=lambda x: (-x["aligned_fraction"], -x["identity"], x["source_id"]))
        status, source = "unassigned", ""
        best = candidates[0] if candidates else dict(aligned_fraction=0, source_completeness=0, identity=0)
        if candidates:
            if best["aligned_fraction"] < cfg["matching"]["min_aligned_fraction"]:
                status = "low_coverage"
            elif len(candidates) > 1 and best["aligned_fraction"] - candidates[1]["aligned_fraction"] <= cfg["matching"]["ambiguity_margin"]:
                status = "ambiguous"
            else:
                status, source = "matched", best["source_id"]
        rows.append(dict(feature_id=mag, source_id=source, match_status=status,
                         aligned_fraction=best["aligned_fraction"], identity=best["identity"],
                         source_completeness=best["source_completeness"]))
    df = pd.DataFrame(rows, columns=["feature_id", "source_id", "match_status", "aligned_fraction", "identity", "source_completeness"])
    frequencies = df.loc[df.match_status == "matched", "source_id"].value_counts()
    df["multiple_bins_per_source"] = df.source_id.map(frequencies).fillna(0).gt(1)
    save_table(df, output)
    save_table(pd.DataFrame(evidence, columns=["feature_id", "source_id", "aligned_fraction", "source_completeness", "identity"]),
               Path(output).with_name("match_evidence.tsv"))


def align_mags(cfg, sources, mags, output, threads):
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with Path(output).open("w") as handle:
        if any(True for _ in fasta(mags)):
            print(version(["minimap2", "--version"]))
            run(["minimap2", "-x", "asm5", "--secondary=yes", "-N", "50", "-t", str(threads),
                 str(sources), str(mags)], cfg, handle)
