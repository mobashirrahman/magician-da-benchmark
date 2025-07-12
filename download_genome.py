#!/usr/bin/env python3
"""
Download a random subset of genomes from the HOMD GTDB taxonomy table.

This script picks *n* random Genome IDs (assembly accessions) from the GTDB taxonomy
CSV downloaded by HOMD and downloads the corresponding genomic FASTA files from the
NCBI Assembly FTP server.

Example
-------
$ python download_genome.py \
    --csv data/HOMD/gtdb_taxonomy/GTDB_taxonomy20250406.csv \
    --n 200 \
    --out HOMD \
    --email you@example.com
"""

from __future__ import annotations

import pathlib, random, time, sys, os
from typing import Optional, List

import pandas as pd
import requests
from Bio import Entrez
import typer

app = typer.Typer(add_completion=False)

# ──────────────────────────────────────────────────────────────
# Helper functions
# ──────────────────────────────────────────────────────────────

def _assembly_uid(accession: str) -> Optional[str]:
    """Return the numeric UID for an assembly accession, or *None* if not found."""
    search_term = f"{accession}[Assembly Accession]"
    handle = Entrez.esearch(db="assembly", term=search_term, retmax=1)
    record = Entrez.read(handle)
    handle.close()
    ids: List[str] = record.get("IdList", [])
    return ids[0] if ids else None


def _ftp_path(uid: str) -> Optional[str]:
    """Given an NCBI Assembly UID, return the GenBank/RefSeq FTP path."""
    handle = Entrez.esummary(db="assembly", id=uid, report="full")
    summary = Entrez.read(handle)
    handle.close()
    docsum = summary["DocumentSummarySet"]["DocumentSummary"][0]
    ftp = docsum.get("FtpPath_GenBank") or docsum.get("FtpPath_RefSeq")
    return ftp if ftp else None


def download_genome(accession: str, out_dir: pathlib.Path, delay: float = 1.0) -> None:
    """Download the *_genomic.fna.gz* file for *accession* into *out_dir*.

    Skips the download if the destination file already exists.
    """
    dest = out_dir / f"{accession}.fna.gz"
    if dest.exists():
        typer.echo(f"· {accession}: already present, skipping")
        return

    uid = _assembly_uid(accession)
    if uid is None:
        typer.echo(f"! {accession}: accession not found on NCBI", err=True)
        return

    ftp = _ftp_path(uid)
    if not ftp:
        typer.echo(f"! {accession}: no public FTP path", err=True)
        return

    base = os.path.basename(ftp)
    # Convert the FTP address to HTTPS because `requests` does not handle the
    # ftp:// scheme.  NCBI serves identical content over HTTPS.
    http_path = ftp.replace("ftp://", "https://")
    url = f"{http_path}/{base}_genomic.fna.gz"

    try:
        with requests.get(url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(dest, "wb") as fh:
                for chunk in r.iter_content(chunk_size=8192):
                    fh.write(chunk)
        typer.echo(f"✓ {accession}")
    except requests.HTTPError as e:
        typer.echo(f"! {accession}: HTTP {e.response.status_code} – {e}", err=True)
    except Exception as e:
        typer.echo(f"! {accession}: {e}", err=True)
    finally:
        time.sleep(delay)


# ──────────────────────────────────────────────────────────────
# CLI entry point
# ──────────────────────────────────────────────────────────────

@app.command()
def main(
    csv: pathlib.Path = typer.Option(
        "data/HOMD/gtdb_taxonomy/GTDB_taxonomy20250406.csv",
        help="Path to GTDB taxonomy CSV (tab-delimited).",
    ),
    n: int = typer.Option(200, help="Number of genomes to sample."),
    out: pathlib.Path = typer.Option(
        "HOMD", help="Output directory where genomes will be stored."
    ),
    email: str = typer.Option(..., help="Email address for NCBI Entrez."),
    api_key: Optional[str] = typer.Option(
        None, help="NCBI API key to increase request limits."
    ),
    seed: Optional[int] = typer.Option(None, help="Random seed for reproducibility."),
):
    """Download *n* random genomes listed in *csv* to *out*."""

    # 0 ▸ setup
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
        delay = 0.34  # ~3 req/s allowed with API key
    else:
        delay = 1.0   # 1 req/s to be polite

    out.mkdir(parents=True, exist_ok=True)

    # 1 ▸ read CSV
    try:
        df = pd.read_csv(csv, sep="\t", comment="#", skiprows=[0])
    except Exception as e:
        typer.echo(f"Error reading {csv}: {e}", err=True)
        sys.exit(1)

    if "Genome-ID" not in df.columns:
        typer.echo("CSV must contain a 'Genome-ID' column.", err=True)
        sys.exit(1)

    accessions = df["Genome-ID"].dropna().unique().tolist()
    if len(accessions) < n:
        typer.echo(f"CSV only contains {len(accessions)} genomes (< {n}).", err=True)
        sys.exit(1)

    rng = random.Random(seed)
    sample = rng.sample(accessions, n)
    typer.echo(f"Selecting {n} genomes out of {len(accessions)} candidates…")

    # 2 ▸ download
    for acc in sample:
        download_genome(acc, out_dir=out, delay=delay)

    typer.echo(f"Done. Genomes saved to {out.resolve()}")


if __name__ == "__main__":
    app() 