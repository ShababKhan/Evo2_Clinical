#!/usr/bin/env python
"""
API Integration Example for Evo2_Clinical

This script demonstrates how to use the API integration module to fetch data
from various biomedical databases.
"""

import os
import sys
import json
import pandas as pd
from pathlib import Path

# Add the project root to Python path to ensure imports work correctly
project_root = Path(__file__).parent.parent.absolute()
sys.path.insert(0, str(project_root))

from evo2_clinical.utils.api_integration import (
    EncodeClient,
    GWASCatalogClient,
    ThousandGenomesClient,
    OpenTargetsClient
)
from evo2_clinical.config import config


def encode_api_example():
    """Example of using the ENCODE API client."""
    print("\n=== ENCODE API Example ===")
    
    # Initialize client
    encode_client = EncodeClient()
    
    # Search for ChIP-seq experiments in endothelial cells
    print("Searching for ChIP-seq experiments in endothelial cells...")
    experiments = encode_client.search_experiments(
        assay_type="ChIP-seq",
        cell_type="endothelial",
        limit=5
    )
    
    print(f"Found {len(experiments)} experiments")
    
    # Display experiment details
    if experiments:
        print("\nExperiment details:")
        experiment = experiments[0]
        print(f"  ID: {experiment.get('accession', 'N/A')}")
        print(f"  Assay: {experiment.get('assay_title', 'N/A')}")
        print(f"  Target: {experiment.get('target', {}).get('label', 'N/A')}")
        print(f"  Description: {experiment.get('description', 'N/A')}")
        
        # Get files for the experiment
        experiment_id = experiment.get('accession')
        if experiment_id:
            print("\nFinding files for experiment...")
            files = encode_client.get_files(
                experiment_id=experiment_id,
                file_type="bigWig"
            )
            
            print(f"Found {len(files)} files")
            
            # Display file details
            if files:
                print("\nFile details:")
                file = files[0]
                print(f"  ID: {file.get('accession', 'N/A')}")
                print(f"  Type: {file.get('file_type', 'N/A')}")
                print(f"  Format: {file.get('file_format', 'N/A')}")
                print(f"  Size: {file.get('file_size', 'N/A')}")
                
                # Note about file download
                print("\nTo download this file:")
                print(f"  encode_client.download_file('{file.get('accession')}', 'output/encode')")


def gwas_catalog_example():
    """Example of using the GWAS Catalog API client."""
    print("\n=== GWAS Catalog API Example ===")
    
    # Initialize client
    gwas_client = GWASCatalogClient()
    
    # Search for studies on lung cancer
    print("Searching for GWAS studies on lung cancer...")
    studies = gwas_client.search_studies(trait="lung cancer", limit=5)
    
    print(f"Found {len(studies)} studies")
    
    # Display study details
    if studies:
        print("\nStudy details:")
        study = studies[0]
        print(f"  Accession: {study.get('accessionId', 'N/A')}")
        print(f"  Title: {study.get('publicationInfo', {}).get('title', 'N/A')}")
        print(f"  Author: {study.get('publicationInfo', {}).get('author', {}).get('fullname', 'N/A')}")
        print(f"  Journal: {study.get('publicationInfo', {}).get('publication', 'N/A')}")
        
        # Get associations for the study
        study_id = study.get('accessionId')
        if study_id:
            print("\nFinding associations for study...")
            associations = gwas_client.get_associations(study_id)
            
            print(f"Found {len(associations)} associations")
            
            # Display association details
            if associations:
                print("\nAssociation details:")
                assoc = associations[0]
                print(f"  Variant: {assoc.get('_links', {}).get('snp', {}).get('href', 'N/A').split('/')[-1]}")
                print(f"  P-value: {assoc.get('pvalue', 'N/A')}")
                print(f"  OR/Beta: {assoc.get('orPerCopyNum', 'N/A')}")
                
                # Search for a specific variant
                print("\nSearching for variant rs17085007...")
                variant_data = gwas_client.search_variants("rs17085007")
                print(f"Found {len(variant_data)} associations for rs17085007")


def thousand_genomes_example():
    """Example of using the 1000 Genomes Project API client."""
    print("\n=== 1000 Genomes Project API Example ===")
    
    # Initialize client
    tg_client = ThousandGenomesClient()
    
    # Get population data
    print("Fetching population data...")
    populations = tg_client.get_populations()
    
    print(f"Found {len(populations)} populations")
    
    # Display population details
    if not populations.empty:
        print("\nPopulation details:")
        for i, (_, pop) in enumerate(populations.head(5).iterrows()):
            print(f"  {pop.get('code', 'N/A')}: {pop.get('text', 'N/A')}")
        
        # Get samples for a specific population
        print("\nFetching samples for GBR (British) population...")
        samples = tg_client.get_samples(population="GBR")
        
        print(f"Found {len(samples)} samples")
        
        # Display sample details
        if not samples.empty:
            print("\nSample details:")
            for i, (_, sample) in enumerate(samples.head(3).iterrows()):
                print(f"  Sample: {sample.get('name', 'N/A')}")
                print(f"  Sex: {sample.get('sex', 'N/A')}")
                print(f"  Population: {sample.get('population', {}).get('code', 'N/A')}")
        
        # Note about variant data
        print("\nTo fetch variant data for a region:")
        print("  tg_client.get_variant_data('1', 10000, 20000)")


def open_targets_example():
    """Example of using the Open Targets Platform API client."""
    print("\n=== Open Targets Platform API Example ===")
    
    # Initialize client
    ot_client = OpenTargetsClient()
    
    # Search for a disease
    print("Searching for lung cancer in Open Targets...")
    diseases = ot_client.search_disease("lung cancer", limit=5)
    
    print(f"Found {len(diseases)} diseases")
    
    # Display disease details
    if diseases:
        print("\nDisease details:")
        disease = diseases[0]
        disease_entity = disease.get('entity', {})
        print(f"  ID: {disease_entity.get('id', 'N/A')}")
        print(f"  Name: {disease_entity.get('name', 'N/A')}")
        
        # Get disease associations
        disease_id = disease_entity.get('id')
        if disease_id:
            print(f"\nFetching target associations for {disease_id}...")
            associations = ot_client.get_disease_associations(disease_id, limit=5)
            
            print(f"Found {len(associations)} target associations")
            
            # Display association details
            if associations:
                print("\nTarget association details:")
                assoc = associations[0]
                target = assoc.get('target', {})
                print(f"  Target ID: {target.get('id', 'N/A')}")
                print(f"  Symbol: {target.get('approvedSymbol', 'N/A')}")
                print(f"  Name: {target.get('approvedName', 'N/A')}")
                print(f"  Association score: {assoc.get('score', 'N/A')}")
                
                # Get target details
                target_id = target.get('id')
                if target_id:
                    print(f"\nFetching details for target {target_id}...")
                    target_details = ot_client.get_target_details(target_id)
                    
                    if target_details:
                        print("\nTarget details:")
                        print(f"  Symbol: {target_details.get('approvedSymbol', 'N/A')}")
                        print(f"  Name: {target_details.get('approvedName', 'N/A')}")
                        print(f"  Biotype: {target_details.get('biotype', 'N/A')}")
                        
                        # Get genomic location
                        location = target_details.get('genomicLocation', {})
                        print(f"\n  Genomic Location:")
                        print(f"    Chromosome: {location.get('chromosome', 'N/A')}")
                        print(f"    Start: {location.get('start', 'N/A')}")
                        print(f"    End: {location.get('end', 'N/A')}")


def main():
    """Main function to run all API examples."""
    print("Evo2_Clinical API Integration Examples")
    print("=====================================")
    
    # Run examples
    try:
        encode_api_example()
    except Exception as e:
        print(f"Error in ENCODE API example: {e}")
    
    try:
        gwas_catalog_example()
    except Exception as e:
        print(f"Error in GWAS Catalog API example: {e}")
    
    try:
        thousand_genomes_example()
    except Exception as e:
        print(f"Error in 1000 Genomes API example: {e}")
    
    try:
        open_targets_example()
    except Exception as e:
        print(f"Error in Open Targets API example: {e}")
    
    print("\nAll examples completed.")


if __name__ == "__main__":
    main()