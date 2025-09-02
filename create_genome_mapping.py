#!/usr/bin/env python3
"""
Create a comprehensive mapping CSV of MetaPhlAn4 taxa to downloaded genomes.

This script combines:
- Original taxonomic names from the MetaPhlAn4 RDS file
- Species search terms used for genome downloads
- Genome status (complete/incomplete)
- Genome accessions (GCA/GCF IDs)
"""

import pandas as pd
import re
import pathlib
import sys
from typing import Optional, Dict, List

# Import our genome download functions
sys.path.append('/home/mrahman/magician')
from download_mooh_representative_genomes import get_species_name, parse_metaphlan_taxon

def extract_genome_info(filename: str) -> tuple:
    """Extract species name and accession from filename."""
    # Filename format: Species_name_ACCESSION.fna.gz
    # Extract accession (GCA_XXXXXX.X or GCF_XXXXXX.X)
    accession_match = re.search(r'(GC[AF]_\d+\.\d+)', filename)
    if accession_match:
        accession = accession_match.group(1)
        # Extract species name (everything before the accession)
        species_part = filename.replace(f'_{accession}.fna.gz', '')
        species_name = species_part.replace('_', ' ')
        return species_name, accession
    return None, None

def determine_genome_status(accession: str) -> str:
    """Determine if genome is complete or incomplete based on accession type."""
    # GCF = RefSeq (curated, often complete genomes)
    # GCA = GenBank (submitted genomes, mixed quality)
    # This is a reasonable heuristic based on NCBI practices
    if accession.startswith('GCF_'):
        return "RefSeq"  # RefSeq genomes (usually higher quality)
    else:
        return "GenBank"  # GenBank genomes (variable quality)

def main():
    # Read the 97% representative species CSV
    print("Reading 97% representative species CSV...")
    rep_species_df = pd.read_csv('/home/mrahman/magician/data/97_representative_species.csv')
    
    # Get list of downloaded genome files
    print("Scanning downloaded genomes...")
    genome_dir = pathlib.Path('/local/mrahman/magician/data/mooh_representative')
    genome_files = list(genome_dir.glob('*.fna.gz'))
    
    # Create a mapping of species names to genome info
    downloaded_genomes = {}
    for genome_file in genome_files:
        species_name, accession = extract_genome_info(genome_file.name)
        if species_name and accession:
            # Normalize species name for matching
            species_key = species_name.lower().replace('_', ' ').strip()
            downloaded_genomes[species_key] = {
                'accession': accession,
                'filename': genome_file.name,
                'status': determine_genome_status(accession)
            }
    
    print(f"Found {len(downloaded_genomes)} downloaded genomes")
    
    # Create comprehensive mapping
    mapping_data = []
    
    for _, row in rep_species_df.iterrows():
        original_taxon = row['taxon']
        abundance = row['abundance']
        relative_abundance = row['relative_abundance']
        
        # Extract species name using our function
        species_search_term = get_species_name(original_taxon)
        
        # Check if we have a downloaded genome for this species
        genome_accession = "Not Available"
        genome_status = "Not Downloaded"
        genome_filename = "Not Available"
        
        if species_search_term:
            species_key = species_search_term.lower().replace('_', ' ').strip()
            if species_key in downloaded_genomes:
                genome_info = downloaded_genomes[species_key]
                genome_accession = genome_info['accession']
                genome_status = genome_info['status']
                genome_filename = genome_info['filename']
            else:
                genome_status = "Download Failed"
        else:
            species_search_term = "Not Extractable"
            genome_status = "Species Not Extractable"
        
        mapping_data.append({
            'original_metaphlan_taxon': original_taxon,
            'species_search_term': species_search_term if species_search_term else "Not Extractable",
            'abundance': abundance,
            'relative_abundance': relative_abundance,
            'genome_status': genome_status,
            'genome_accession': genome_accession,
            'genome_filename': genome_filename
        })
    
    # Create DataFrame and save
    mapping_df = pd.DataFrame(mapping_data)
    
    # Sort by relative abundance (descending)
    mapping_df = mapping_df.sort_values('relative_abundance', ascending=False)
    
    output_file = '/home/mrahman/magician/data/metaphlan4_genome_mapping.csv'
    mapping_df.to_csv(output_file, index=False)
    
    print(f"\n=== Mapping Summary ===")
    print(f"Total taxa processed: {len(mapping_df)}")
    print(f"Species extractable: {len(mapping_df[mapping_df['species_search_term'] != 'Not Extractable'])}")
    print(f"Genomes downloaded: {len(mapping_df[mapping_df['genome_status'].isin(['RefSeq', 'GenBank'])])}")
    print(f"Download failures: {len(mapping_df[mapping_df['genome_status'] == 'Download Failed'])}")
    print(f"Non-extractable species: {len(mapping_df[mapping_df['genome_status'] == 'Species Not Extractable'])}")
    
    print(f"\nGenome Status Distribution:")
    status_counts = mapping_df['genome_status'].value_counts()
    for status, count in status_counts.items():
        print(f"  {status}: {count}")
    
    print(f"\nMapping file saved: {output_file}")
    
    # Show a preview
    print(f"\nPreview of mapping (top 10 taxa):")
    preview_cols = ['original_metaphlan_taxon', 'species_search_term', 'relative_abundance', 'genome_status', 'genome_accession']
    print(mapping_df[preview_cols].head(10).to_string(index=False))

if __name__ == "__main__":
    main()
