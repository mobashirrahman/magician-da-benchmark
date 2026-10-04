"""Versioned real-genome manifests with an explicit, recorded selection rule.

A manifest is immutable evidence: accession and version, local FASTA path, checksum,
length and provenance. Nothing is downloaded automatically. Existing local references
are reused by path after their provenance and completeness are checked, so a large
historical collection is never copied into this project.
"""
from pathlib import Path
import re
import shutil
import pandas as pd
from .config import resolve
from .io import table, save_table, sha256, fasta, write_json

REQUIRED = ("genome_id", "accession", "version", "path", "provenance", "source_group")
ACCESSION = re.compile(r"^(GC[AF]_\d+|NC_\d+\.\d+|[A-Z]{2}_\d+|U_\d+|GCA_\d+\.\d+|RS_\d+)$")
# Reads shorter than this cannot be assembled by this pipeline.
MINIMUM_LENGTH = 10000


def check_genomes(entries, maximum_n_content=0.05, minimum_length=MINIMUM_LENGTH):
    """Reject unusable references explicitly rather than failing inside assembly."""
    problems = []
    rows = []
    for entry in entries:
        path = Path(entry["path"])
        if not path.exists():
            problems.append(f"{entry['genome_id']}: missing file {path}")
            continue
        if not ACCESSION.match(str(entry["accession"])):
            problems.append(f"{entry['genome_id']}: accession {entry['accession']} is not a recognised accession")
        length, contigs, n_content = 0, 0, 0
        for _, sequence in fasta(path):
            contigs += 1
            length += len(sequence)
            n_content += sequence.count("N")
        if length < minimum_length:
            problems.append(f"{entry['genome_id']}: {length} bp is below the {minimum_length} bp minimum")
        if length and n_content / length > maximum_n_content:
            problems.append(f"{entry['genome_id']}: more than {maximum_n_content:.0%} N bases")
        rows.append(dict(genome_id=entry["genome_id"], accession=entry["accession"],
                         version=entry["version"], path=str(path.resolve()),
                         length_bp=length, n_contigs=contigs,
                         n_fraction=n_content / length if length else None,
                         provenance=entry["provenance"], source_group=entry["source_group"],
                         sha256=sha256(path)))
    return rows, problems


def select_manifest(entries, group_sizes, selection_rule, require_groups=True):
    """Apply a declared selection rule to a candidate pool.

    `group_sizes` names how many genomes to take from each source group. The rule is
    recorded verbatim so the selection is reproducible and auditable.
    """
    chosen, problems = [], []
    for group, wanted in group_sizes.items():
        candidates = sorted((entry for entry in entries if entry["source_group"] == group),
                            key=lambda entry: str(entry["genome_id"]))
        if require_groups and len(candidates) < wanted:
            problems.append(f"{group}: {len(candidates)} candidates, {wanted} required")
            continue
        chosen.extend(candidates[:wanted])
    ids = [entry["genome_id"] for entry in chosen]
    if len(set(ids)) != len(ids):
        problems.append("selection produced duplicate genome IDs")
    return chosen, problems, selection_rule


def write_manifest(rows, output, selection_rule, metadata=None):
    frame = pd.DataFrame(rows).sort_values("genome_id")
    if frame.genome_id.duplicated().any():
        raise ValueError("Manifest genome IDs must be unique")
    save_table(frame, output)
    write_json(Path(output).with_suffix(".json"), dict(
        n_genomes=len(frame), selection_rule=selection_rule,
        lengths=dict(min=int(frame.length_bp.min()) if len(frame) else None,
                     max=int(frame.length_bp.max()) if len(frame) else None,
                     total=int(frame.length_bp.sum()) if len(frame) else None),
        gc_span=metadata or {}, sha256=sha256(output)))
    return frame


def link_reference(source, target, genome_id):
    """Reuse a suitable archived reference by symlink, never by copying a collection."""
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() or target.is_symlink():
        if target.resolve() != Path(source).resolve():
            raise ValueError(f"{target} already exists and points elsewhere")
        return target
    target.symlink_to(Path(source).resolve())
    return target