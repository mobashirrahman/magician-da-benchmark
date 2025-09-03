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


def _is_valid_strain_token(token: str) -> bool:
    """Heuristic to decide if MetaPhlAn t__ token likely represents a real strain identifier."""
    if not token:
        return False
    t = token.strip()
    # Exclude common placeholders/codes
    bad_markers = ["SGB", "group", "unclassified", "EUK", "GGB", "_SGB", "_group"]
    if any(m.lower() in t.lower() for m in bad_markers):
        return False
    # Avoid generic sp_ tokens
    if re.match(r"^sp[_-]", t, flags=re.IGNORECASE):
        return False
    return True


def get_strain_candidate(taxon_string: str) -> Optional[str]:
    """Extract strain-like token from MetaPhlAn taxon (t__) if it appears valid."""
    taxonomy = parse_metaphlan_taxon(taxon_string)
    strain = taxonomy.get("strain")
    if strain and _is_valid_strain_token(strain):
        return strain
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


def search_genomes_for_strain(species: str, strain: str, max_results: int = 10) -> List[str]:
    """Search for strain-matching assemblies within a species (best-effort via All Fields)."""
    search_term = f'"{species}"[Organism] AND ("{strain}"[All Fields] OR "{species} {strain}"[All Fields])'
    try:
        handle = Entrez.esearch(
            db="assembly",
            term=search_term,
            retmax=max_results,
            sort="relevance",
        )
        record = Entrez.read(handle)
        handle.close()
        return record.get("IdList", [])
    except Exception as e:
        typer.echo(f"Error searching strain {species} / {strain}: {e}", err=True)
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

def _assembly_rank(accession: Optional[str], level: Optional[str]) -> tuple[int, int]:
    """Ranking helper: prefer RefSeq over GenBank, and higher assembly levels."""
    refseq_score = 1 if accession and str(accession).startswith("GCF_") else 0
    level_rank = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}
    return (refseq_score, level_rank.get(level or "", 0))


def pick_best_assembly(assembly_ids: List[str]) -> Optional[Tuple[str, str, str]]:
    """Return (uid, accession, level) for best assembly among provided IDs."""
    best = None
    best_score = (-1, -1)
    for uid in assembly_ids:
        acc, lvl = get_assembly_info(uid)
        score = _assembly_rank(acc, lvl)
        if acc and score > best_score:
            best = (uid, acc, lvl or "Unknown")
            best_score = score
    return best


def download_genome_by_uid(uid: str, species: str, out_dir: pathlib.Path, delay: float = 1.0, strain: Optional[str] = None) -> bool:
    """Download genome using NCBI Assembly UID."""
    
    # Get accession and assembly level for filename and info
    accession, assembly_level = get_assembly_info(uid)
    if not accession:
        typer.echo(f"! {species}: Could not get accession for UID {uid}", err=True)
        return False
    
    # Create safe filename from species (and optional strain) name
    safe_species = re.sub(r'[^\w\s-]', '', species).strip()
    safe_species = re.sub(r'[-\s]+', '_', safe_species)
    safe_strain = None
    if strain:
        st = re.sub(r'[^\w\s-]', '', strain).strip()
        st = re.sub(r'[-\s]+', '_', st)
        safe_strain = st if st else None
    base_name = f"{safe_species}{('_' + safe_strain) if safe_strain else ''}_{accession}"
    dest = out_dir / f"{base_name}.fna.gz"
    
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
        label = f"{species}{(' ' + strain) if strain else ''}"
        typer.echo(f"✓ {label}: {accession} [{level_indicator} {assembly_level}]")
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
    out_mapping: Optional[pathlib.Path] = typer.Option(
        pathlib.Path("data/metaphlan4_genome_mapping_strain.csv"),
        help="Optional path to write mapping CSV with chosen assemblies",
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

    # Extract species + strain candidates
    species_data = []
    for _, row in df_filtered.iterrows():
        taxon = row["taxon"]
        species_name = get_species_name(taxon)
        if species_name:
            species_data.append({
                "species": species_name,
                "strain": get_strain_candidate(taxon),
                "taxon": taxon,
                "abundance": row["relative_abundance"]
            })

    typer.echo(f"Extracted {len(species_data)} valid species names")

    # Download genomes with strain preference; ensure uniqueness by (species, strain)
    successful_downloads = 0
    failed_species = []
    seen_keys = set()
    mapping_rows: List[Dict[str, str]] = []

    for i, data in enumerate(species_data, 1):
        species = data["species"]
        strain = data.get("strain")
        abundance = data["abundance"]

        key = (species, strain or None)
        if key in seen_keys:
            typer.echo(f"· {species}{(' ' + strain) if strain else ''}: already processed, skipping")
            continue

        typer.echo(f"[{i}/{len(species_data)}] Processing {species}{(' ' + strain) if strain else ''} ({abundance:.2f}%)")

        # Search for genomes (strain-first)
        assembly_ids: List[str] = []
        chosen_level = "species"
        if strain:
            assembly_ids = search_genomes_for_strain(species, strain, max_results=max_per_species * 3)
            if assembly_ids:
                chosen_level = "strain"
        if not assembly_ids:
            assembly_ids = search_genomes_for_species(species, max_results=max_per_species * 3)
        
        if not assembly_ids:
            typer.echo(f"! {species}{(' ' + strain) if strain else ''}: No assemblies found")
            failed_species.append(species)
            time.sleep(delay)
            continue

        best = pick_best_assembly(assembly_ids)
        if not best:
            typer.echo(f"! {species}{(' ' + strain) if strain else ''}: No suitable assemblies after ranking")
            failed_species.append(species)
            continue

        uid, accession, level = best
        success = download_genome_by_uid(uid, species, out, delay, strain=strain)
        if success:
            successful_downloads += 1
            seen_keys.add(key)
            # Build output filename consistently with download_genome_by_uid
            safe_species = re.sub(r'[^\w\s-]', '', species).strip()
            safe_species = re.sub(r'[-\s]+', '_', safe_species)
            safe_strain = None
            if strain:
                st = re.sub(r'[^\w\s-]', '', strain).strip()
                st = re.sub(r'[-\s]+', '_', st)
                safe_strain = st if st else None
            filename = f"{safe_species}{('_' + safe_strain) if safe_strain else ''}_{accession}.fna.gz"
            genome_status = "RefSeq" if accession.startswith("GCF_") else "GenBank"
            mapping_rows.append({
                "original_metaphlan_taxon": data["taxon"],
                "species_search_term": species,
                "strain_candidate": strain or "",
                "chosen_level": chosen_level,
                "assembly_level": level or "Unknown",
                "abundance": data["abundance"],
                "relative_abundance": data["abundance"],
                "genome_status": genome_status,
                "genome_accession": accession,
                "genome_filename": filename,
            })
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

    # Write updated mapping if requested
    if out_mapping:
        try:
            out_mapping.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(mapping_rows).to_csv(out_mapping, index=False)
            typer.echo(f"\nUpdated mapping written to: {out_mapping.resolve()}")
        except Exception as e:
            typer.echo(f"Error writing mapping CSV: {e}", err=True)

    typer.echo(f"\nGenomes saved to: {out.resolve()}")


if __name__ == "__main__":
    app()
