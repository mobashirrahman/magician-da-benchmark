"""Small streaming IO helpers; no copies of external genome libraries."""
from pathlib import Path
import gzip
import hashlib
import json
import os
import tempfile
import pandas as pd


def write_json(path, obj):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=p.parent, delete=False) as handle:
        json.dump(obj, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        name = handle.name
    os.replace(name, p)


def read_json(path):
    return json.loads(Path(path).read_text())


def table(path):
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def save_table(df, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")


def fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    name, parts = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            elif line.strip():
                if name is None:
                    raise ValueError(f"Invalid FASTA: {path}")
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts).upper()


def write_fasta(records, path):
    with Path(path).open("w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for i in range(0, len(sequence), 80):
                handle.write(sequence[i:i + 80] + "\n")


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def tree_bytes(root):
    root = Path(root)
    if not root.exists():
        return 0
    total = 0
    for base, dirs, files in os.walk(root, followlinks=False):
        dirs[:] = [d for d in dirs if not (Path(base) / d).is_symlink()]
        for name in files:
            p = Path(base) / name
            try:
                if not p.is_symlink():
                    total += p.stat().st_size
            except FileNotFoundError:
                pass
    return total
