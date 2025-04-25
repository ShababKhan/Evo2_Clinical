#!/usr/bin/env python
"""
Basic example of using the Evo2_Clinical package.

This script demonstrates the basic workflow of loading data,
analyzing variants, and storing results.
"""

import os
import sys
import logging
from pathlib import Path

# Add the parent directory to the system path to import the package
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from evo2_clinical.utils.helpers import setup_logging
from evo2_clinical.config import Config
from evo2_clinical.data.io import load_vcf_file, load_gene_list
from evo2_clinical.analysis.variant_scoring import run_evo2_variant_scoring
from evo2_clinical.database.manager import DatabaseManager, create_endothelial_variants_database


def main():
    # Set up logging
    log_file = setup_logging(os.path.join(os.path.dirname(__file__), 'logs'))
    logger = logging.getLogger(__name__)
    logger.info("Starting Evo2_Clinical basic example")

    # Create a configuration
    config = Config()
    
    # Paths to example data files
    data_dir = os.path.join(os.path.dirname(__file__), '..', 'data')
    vcf_file = os.path.join(data_dir, '1000genomes', 'chr19_variants.vcf')
    
    # Load gene list (creating a simple one for this example)
    endothelial_genes = ['GATA2', 'KLF2', 'KLF4', 'FOXC1', 'ERG']
    logger.info(f"Using example endothelial genes: {', '.join(endothelial_genes)}")
    
    # Load VCF data if file exists
    if os.path.exists(vcf_file):
        try:
            logger.info(f"Loading VCF data from {vcf_file}")
            variants_df = load_vcf_file(vcf_file)
            logger.info(f"Loaded {len(variants_df)} variants")
            
            # Filter variants (this is a simplified example)
            # In a real scenario, you would filter based on gene coordinates
            filtered_variants = variants_df.head(100)  # Just take first 100 for example
            logger.info(f"Using {len(filtered_variants)} variants for demonstration")
            
            # Score variants
            logger.info("Scoring variants (simulated - no actual Evo2 executable)")
            # Note: The actual scoring would require the Evo2 executable
            # Here we'll add a mock functional score column
            scored_variants = filtered_variants.copy()
            scored_variants['functional_score'] = pd.Series([0.1 * i for i in range(len(scored_variants))])
            
            # Store in database
            logger.info("Creating and populating database")
            db_dir = os.path.join(os.path.dirname(__file__), 'output', 'databases')
            os.makedirs(db_dir, exist_ok=True)
            
            db_manager = DatabaseManager(db_dir)
            create_endothelial_variants_database(db_manager)
            
            # Modify the DataFrame to match expected schema
            if 'CHROM' not in scored_variants.columns and '#CHROM' in scored_variants.columns:
                scored_variants = scored_variants.rename(columns={'#CHROM': 'chrom'})
            else:
                scored_variants['chrom'] = scored_variants.get('#CHROM', 'chr1')
            
            scored_variants['gene'] = 'EXAMPLE'
            scored_variants['is_common'] = 1
            scored_variants['is_rare'] = 0
            scored_variants['allele_freq'] = 0.1
            
            selected_cols = [
                'chrom', 'POS', 'REF', 'ALT', 'gene', 
                'is_common', 'is_rare', 'allele_freq', 'functional_score'
            ]
            
            # Remove id column if it exists to avoid conflict with database schema
            if 'id' in scored_variants.columns:
                scored_variants = scored_variants.drop(columns=['id'])
            
            # Only include required columns if they exist
            valid_cols = [col for col in selected_cols if col in scored_variants.columns]
            db_variants = scored_variants[valid_cols].rename(columns={'POS': 'pos', 'REF': 'ref', 'ALT': 'alt'})
            
            db_manager.populate_database('endothelial_variants_db', db_variants, 'variants')
            logger.info("Database populated successfully")
            
            # Query the database to verify
            query = "SELECT * FROM variants LIMIT 5"
            results = db_manager.query_database('endothelial_variants_db', query)
            logger.info(f"Database query result: {results.shape[0]} rows")
            logger.info(f"Columns: {', '.join(results.columns)}")
            
        except Exception as e:
            logger.error(f"Error during example execution: {e}", exc_info=True)
    else:
        logger.warning(f"VCF file not found: {vcf_file}")
        logger.info("This is just an example script. You need to provide actual data files to run the full pipeline.")
    
    logger.info("Example completed")


if __name__ == '__main__':
    main()