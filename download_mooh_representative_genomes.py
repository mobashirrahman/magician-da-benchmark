#!/usr/bin/env python3
"""
Download representative genomes for taxa from MetaPhlAn4 97% representative species.

This script reads the 97_representative_species.csv file containing MetaPhlAn4 taxonomic
names and downloads representative genomes from NCBI for each species. Uses a fallback
strategy: first tries complete genomes, then reference/representative genomes, finally
any available genome for maximum coverage.

Example
-------
$ python download_mooh_representative_genomes.py \
    --csv data/97_representative_species.csv \
    --out /local/mrahman/magician/data/mooh_representative \
    --email you@example.com
"""

from __future__ import annotations

import pathlib, time, sys, os, re
from typing import Optional, List, Dict, Tuple

import pandas as pd
import requests
from Bio import Entrez
import typer

app = typer.Typer(add_completion=False)

# ──────────────────────────────────────────────────────────────
# Helper functions
# ──────────────────────────────────────────────────────────────

def parse_metaphlan_taxon(taxon_string: str) -> Dict[str, str]:
    """Parse MetaPhlAn4 taxonomic string into components."""
    taxonomy = {}
    parts = taxon_string.split("|")
    
    for part in parts:
        if part.startswith("k__"):
            taxonomy["kingdom"] = part[3:]
        elif part.startswith("p__"):
            taxonomy["phylum"] = part[3:]
        elif part.startswith("c__"):
            taxonomy["class"] = part[3:]
        elif part.startswith("o__"):
            taxonomy["order"] = part[3:]
        elif part.startswith("f__"):
            taxonomy["family"] = part[3:]
        elif part.startswith("g__"):
            taxonomy["genus"] = part[3:]
        elif part.startswith("s__"):
            taxonomy["species"] = part[3:]
        elif part.startswith("t__"):
            taxonomy["strain"] = part[3:]
    
    return taxonomy

def get_species_name(taxon_string: str) -> Optional[str]:
    """Extract species name from MetaPhlAn4 taxonomic string."""
    taxonomy = parse_metaphlan_taxon(taxon_string)
    
    if "genus" in taxonomy and "species" in taxonomy:
        genus = taxonomy["genus"]
        species = taxonomy["species"]
        
        # Handle cases where species contains strain info or SGB codes
        species_clean = species.replace(f"{genus}_", "")
        species_clean = re.sub(r"_SGB\d+.*", "", species_clean)
        
        # Skip if species is just an SGB code or unclear
        if species_clean.startswith("SGB") or not species_clean:
            return None
            
        return f"{genus} {species_clean}"
    
    return None

def search_genomes_for_species(species: str, max_results: int = 5) -> List[str]:
    """Search for genomes for a given species, preferring complete genomes."""
    # First try to find complete genomes
    search_term = f'"{species}"[Organism] AND "complete genome"[Assembly Level]'
    
    try:
        handle = Entrez.esearch(
            db="assembly", 
            term=search_term, 
            retmax=max_results,
            sort="relevance"
        )
        record = Entrez.read(handle)
        handle.close()
        
        assembly_ids = record.get("IdList", [])
        
        # If no complete genomes found, search for any representative genome
        if not assembly_ids:
            typer.echo(f"  No complete genomes found for {species}, searching for any representative genome...")
            search_term = f'"{species}"[Organism] AND ("reference genome"[Assembly Level] OR "representative genome"[Filter] OR "latest"[Filter])'
            
            handle = Entrez.esearch(
                db="assembly", 
                term=search_term, 
                retmax=max_results,
                sort="relevance"
            )
            record = Entrez.read(handle)
            handle.close()
            
            assembly_ids = record.get("IdList", [])
            
            # If still no results, try broader search without assembly level filter
            if not assembly_ids:
                typer.echo(f"  No representative genomes found, searching for any genome...")
                search_term = f'"{species}"[Organism]'
                
                handle = Entrez.esearch(
                    db="assembly", 
                    term=search_term, 
                    retmax=max_results,
                    sort="relevance"
                )
                record = Entrez.read(handle)
                handle.close()
                
                assembly_ids = record.get("IdList", [])
        
        return assembly_ids
    except Exception as e:
        typer.echo(f"Error searching for {species}: {e}", err=True)
        return []



def _ftp_path(uid: str) -> Optional[str]:
    """Given an NCBI Assembly UID, return the GenBank/RefSeq FTP path."""
    try:
        handle = Entrez.esummary(db="assembly", id=uid, report="full")
        summary = Entrez.read(handle)
        handle.close()
        docsum = summary["DocumentSummarySet"]["DocumentSummary"][0]
        ftp = docsum.get("FtpPath_GenBank") or docsum.get("FtpPath_RefSeq")
        return ftp if ftp else None
    except Exception as e:
        typer.echo(f"Error getting FTP path for UID {uid}: {e}", err=True)
        return None

def get_assembly_info(uid: str) -> Tuple[Optional[str], Optional[str]]:
    """Get assembly accession and assembly level from UID."""
    try:
        handle = Entrez.esummary(db="assembly", id=uid, report="full")
        summary = Entrez.read(handle)
        handle.close()
        
        docsum = summary["DocumentSummarySet"]["DocumentSummary"][0]
        accession = docsum.get("AssemblyAccession")
        assembly_level = docsum.get("AssemblyLevel", "Unknown")
        return accession, assembly_level
    except Exception as e:
        typer.echo(f"Error getting assembly info for UID {uid}: {e}", err=True)
        return None, None

def download_genome_by_uid(uid: str, species: str, out_dir: pathlib.Path, delay: float = 1.0) -> bool:
    """Download genome using NCBI Assembly UID."""
    
    # Get accession and assembly level for filename and info
    accession, assembly_level = get_assembly_info(uid)
    if not accession:
        typer.echo(f"! {species}: Could not get accession for UID {uid}", err=True)
        return False
    
    # Create safe filename from species name
    safe_species = re.sub(r'[^\w\s-]', '', species).strip()
    safe_species = re.sub(r'[-\s]+', '_', safe_species)
    dest = out_dir / f"{safe_species}_{accession}.fna.gz"
    
    if dest.exists():
        typer.echo(f"· {species}: already present, skipping")
        return True

    ftp = _ftp_path(uid)
    if not ftp:
        typer.echo(f"! {species}: no public FTP path for {accession}", err=True)
        return False

    base = os.path.basename(ftp)
    # Convert FTP to HTTPS
    http_path = ftp.replace("ftp://", "https://")
    url = f"{http_path}/{base}_genomic.fna.gz"

    try:
        with requests.get(url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(dest, "wb") as fh:
                for chunk in r.iter_content(chunk_size=8192):
                    fh.write(chunk)
        
        # Show assembly level in output
        level_indicator = "●" if assembly_level == "Complete Genome" else "○"
        typer.echo(f"✓ {species}: {accession} [{level_indicator} {assembly_level}]")
        return True
    except requests.HTTPError as e:
        typer.echo(f"! {species}: HTTP {e.response.status_code} – {e}", err=True)
        return False
    except Exception as e:
        typer.echo(f"! {species}: {e}", err=True)
        return False
    finally:
        time.sleep(delay)

# ──────────────────────────────────────────────────────────────
# CLI entry point
# ──────────────────────────────────────────────────────────────

@app.command()
def main(
    csv: pathlib.Path = typer.Option(
        "data/97_representative_species.csv",
        help="Path to MetaPhlAn4 representative species CSV.",
    ),
    out: pathlib.Path = typer.Option(
        "/local/mrahman/magician/data/mooh_representative",
        help="Output directory where genomes will be stored."
    ),
    email: str = typer.Option(..., help="Email address for NCBI Entrez."),
    api_key: Optional[str] = typer.Option(
        None, help="NCBI API key to increase request limits."
    ),
    max_per_species: int = typer.Option(
        1, help="Maximum number of genomes to download per species."
    ),
    min_abundance: float = typer.Option(
        0.0, help="Minimum relative abundance threshold (%)."
    ),
):
    """Download representative genomes for species from MetaPhlAn4 data."""

    # Setup
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
        delay = 0.34  # ~3 req/s allowed with API key
    else:
        delay = 1.0   # 1 req/s to be polite

    out.mkdir(parents=True, exist_ok=True)

    # Read CSV
    try:
        df = pd.read_csv(csv)
    except Exception as e:
        typer.echo(f"Error reading {csv}: {e}", err=True)
        sys.exit(1)

    if "taxon" not in df.columns or "relative_abundance" not in df.columns:
        typer.echo("CSV must contain 'taxon' and 'relative_abundance' columns.", err=True)
        sys.exit(1)

    # Filter by abundance threshold
    df_filtered = df[df["relative_abundance"] >= min_abundance]
    typer.echo(f"Found {len(df_filtered)} taxa above {min_abundance}% abundance threshold")

    # Extract species names
    species_data = []
    for _, row in df_filtered.iterrows():
        species_name = get_species_name(row["taxon"])
        if species_name:
            species_data.append({
                "species": species_name,
                "taxon": row["taxon"],
                "abundance": row["relative_abundance"]
            })

    typer.echo(f"Extracted {len(species_data)} valid species names")

    # Download genomes
    successful_downloads = 0
    failed_species = []

    for i, data in enumerate(species_data, 1):
        species = data["species"]
        abundance = data["abundance"]
        
        typer.echo(f"[{i}/{len(species_data)}] Processing {species} ({abundance:.2f}%)")
        
        # Search for genomes
        assembly_ids = search_genomes_for_species(species, max_results=max_per_species * 2)
        
        if not assembly_ids:
            typer.echo(f"! {species}: No complete genomes found")
            failed_species.append(species)
            time.sleep(delay)
            continue
        
        # Try to download up to max_per_species genomes
        downloaded_count = 0
        for uid in assembly_ids:
            if downloaded_count >= max_per_species:
                break
                
            success = download_genome_by_uid(uid, species, out, delay)
            if success:
                downloaded_count += 1
        
        if downloaded_count > 0:
            successful_downloads += 1
        else:
            failed_species.append(species)

    # Summary
    typer.echo(f"\n=== Download Summary ===")
    typer.echo(f"Species processed: {len(species_data)}")
    typer.echo(f"Successful downloads: {successful_downloads}")
    typer.echo(f"Failed species: {len(failed_species)}")
    
    if failed_species:
        typer.echo(f"\nFailed species:")
        for species in failed_species[:10]:  # Show first 10
            typer.echo(f"  - {species}")
        if len(failed_species) > 10:
            typer.echo(f"  ... and {len(failed_species) - 10} more")

    typer.echo(f"\nGenomes saved to: {out.resolve()}")


if __name__ == "__main__":
    app()
