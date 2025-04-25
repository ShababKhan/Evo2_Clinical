"""
Main pipeline module for the Evo2_Clinical package.

This module orchestrates the execution of the complete Evo2 computational pipeline,
managing data loading, analysis, and database generation.
"""

import os
import pandas as pd
import logging
from typing import List, Dict, Union, Optional, Tuple, Any
import time

from ..config import config
from ..data import io
from ..analysis import variant_scoring, lncrna_scoring
from ..database.manager import DatabaseManager, create_all_databases

logger = logging.getLogger(__name__)


class Pipeline:
    """
    Main class for orchestrating the Evo2 computational pipeline.
    """
    
    def __init__(self, config_file: Optional[str] = None, db_manager: Optional[DatabaseManager] = None):
        """
        Initialize the pipeline with configuration and database manager.
        
        Args:
            config_file: Optional path to a configuration file.
            db_manager: Optional DatabaseManager instance.
        """
        # Load configuration if specified
        if config_file:
            config.load_config(config_file)
        
        # Initialize database manager
        self.db_manager = db_manager or DatabaseManager()
        
        # Ensure output directories exist
        config.ensure_output_dirs()
        
        # Initialize logging
        self._setup_logging()
    
    def run_pipeline(self) -> bool:
        """
        Execute the complete Evo2 computational pipeline.
        
        Returns:
            bool: True if successful, False if there were errors.
        """
        logger.info("Starting Evo2 computational pipeline execution...")
        start_time = time.time()
        
        try:
            # Step 1: Data Loading and Integration
            data_loading_success = self._data_loading_step()
            if not data_loading_success:
                logger.error("Data loading step failed, stopping pipeline")
                return False
            
            # Step 2: Variant and lncRNA Scoring
            scoring_success = self._scoring_step()
            if not scoring_success:
                logger.error("Scoring step failed, stopping pipeline")
                return False
            
            # Step 3: Database Generation
            database_success = self._database_generation_step()
            if not database_success:
                logger.error("Database generation step failed")
                # Continue despite errors, as partial results may still be useful
            
            # Step 4: Specific Analyses
            analysis_success = self._specific_analyses_step()
            if not analysis_success:
                logger.error("Specific analyses step failed")
                # Continue despite errors, as partial results may still be useful
            
            execution_time = time.time() - start_time
            logger.info(f"Evo2 computational pipeline execution completed in {execution_time:.2f} seconds")
            return True
            
        except Exception as e:
            logger.error(f"Pipeline execution failed with error: {e}", exc_info=True)
            return False
    
    def _setup_logging(self) -> None:
        """Setup logging configuration."""
        log_dir = config.get_absolute_path(config.OUTPUT_DIRS.get('logs', 'output/logs'))
        os.makedirs(log_dir, exist_ok=True)
        
        log_file = os.path.join(log_dir, f"pipeline_{time.strftime('%Y%m%d_%H%M%S')}.log")
        
        # Configure root logger
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.INFO)
        
        # File handler
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)
        file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(file_formatter)
        root_logger.addHandler(file_handler)
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_formatter = logging.Formatter('%(levelname)s - %(message)s')
        console_handler.setFormatter(console_formatter)
        root_logger.addHandler(console_handler)
        
        logger.info(f"Logging to {log_file}")
    
    def _data_loading_step(self) -> bool:
        """
        Load and process data from various sources.
        
        Returns:
            bool: True if successful, False otherwise.
        """
        logger.info("Starting data loading and integration step")
        
        try:
            # Load gene lists
            self.endothelial_genes = io.load_gene_list(
                config.DATA_PATHS.get('endothelial_genes_list', '')
            )
            logger.info(f"Loaded {len(self.endothelial_genes)} endothelial genes")
            
            self.lncrna_list = io.load_lncrna_list(
                config.DATA_PATHS.get('lncrna_list', '')
            )
            logger.info(f"Loaded {len(self.lncrna_list)} lncRNAs")
            
            self.emt_genes = io.load_emt_genes(
                config.DATA_PATHS.get('emt_genes_list', '')
            )
            logger.info(f"Loaded {len(self.emt_genes)} EMT genes")
            
            # Load GWAS data
            self.gwas_data = io.load_gwas_catalog_data(
                config.DATA_PATHS.get('gwas_catalog', '')
            )
            logger.info(f"Loaded {len(self.gwas_data)} GWAS catalog entries")
            
            # Load and filter 1000 Genomes data
            self.all_1000_genomes_variants = io.load_1000_genomes_data(
                config.DATA_PATHS.get('1000_genomes_vcf', ''),
                self.endothelial_genes
            )
            logger.info(f"Loaded {len(self.all_1000_genomes_variants)} variants from 1000 Genomes")
            
            # Filter variants by ENCODE data for endothelial cells
            self.endothelial_active_variants = io.filter_by_encode_data(
                self.all_1000_genomes_variants,
                config.DATA_PATHS.get('encode_data', ''),
                cell_type="ENDOS"
            )
            logger.info(f"Filtered to {len(self.endothelial_active_variants)} variants active in endothelial cells")
            
            # Define specialized variant subsets
            # These will be populated in later steps
            self.gata2_as1_variants = pd.DataFrame()
            self.epigenetic_mediator_variants = pd.DataFrame()
            self.emt_variants = pd.DataFrame()
            self.cteph_gwas_variants = pd.DataFrame()
            self.pah_known_variants = pd.DataFrame()
            self.mesothelioma_variants = pd.DataFrame()
            
            logger.info("Data loading and integration completed successfully")
            return True
            
        except Exception as e:
            logger.error(f"Error in data loading step: {e}", exc_info=True)
            return False
    
    def _scoring_step(self) -> bool:
        """
        Score variants and lncRNAs using Evo2.
        
        Returns:
            bool: True if successful, False otherwise.
        """
        logger.info("Starting variant and lncRNA scoring step")
        
        try:
            # Score endothelial variants
            self.scored_endothelial_variants = variant_scoring.run_evo2_variant_scoring(
                self.endothelial_active_variants,
                evo2_path=config.TOOL_PATHS.get('evo2_executable', '')
            )
            logger.info(f"Scored {len(self.scored_endothelial_variants)} endothelial variants")
            
            # Score lncRNAs
            self.scored_lncrnas = lncrna_scoring.run_evo2_lncrna_scoring(
                self.lncrna_list,
                evo2_path=config.TOOL_PATHS.get('evo2_executable', '')
            )
            logger.info(f"Scored {len(self.scored_lncrnas)} lncRNAs")
            
            # Extract GATA2-AS1 variants for specialized analysis
            # This is a placeholder implementation - in reality, you would need to identify these variants
            # based on genomic coordinates or annotations
            self.gata2_as1_variants = self.endothelial_active_variants[
                self.endothelial_active_variants['INFO'].str.contains('GATA2-AS1', na=False)
            ]
            logger.info(f"Identified {len(self.gata2_as1_variants)} variants in GATA2-AS1")
            
            logger.info("Variant and lncRNA scoring completed successfully")
            return True
            
        except Exception as e:
            logger.error(f"Error in scoring step: {e}", exc_info=True)
            return False
    
    def _database_generation_step(self) -> bool:
        """
        Create and populate databases with analysis results.
        
        Returns:
            bool: True if successful, False otherwise.
        """
        logger.info("Starting database generation step")
        
        try:
            # Create all database schemas
            create_all_databases(self.db_manager)
            
            # Populate endothelial variants database
            self.db_manager.populate_database(
                'endothelial_variants_db',
                self.scored_endothelial_variants,
                table_name='variants'
            )
            
            # Populate lncRNA functionality database
            self.db_manager.populate_database(
                'lncrna_functionality_db',
                self.scored_lncrnas,
                table_name='lncrnas'
            )
            
            logger.info("Database generation completed successfully")
            return True
            
        except Exception as e:
            logger.error(f"Error in database generation step: {e}", exc_info=True)
            return False
    
    def _specific_analyses_step(self) -> bool:
        """
        Perform specialized analyses for specific targets.
        
        Returns:
            bool: True if successful, False otherwise.
        """
        logger.info("Starting specific analyses step")
        
        try:
            # GATA2-AS1 variant effect predictions
            self.gata2_as1_predictions = variant_scoring.predict_gata2_as1_variant_effects(
                self.gata2_as1_variants,
                evo2_path=config.TOOL_PATHS.get('evo2_executable', '')
            )
            logger.info(f"Generated effect predictions for {len(self.gata2_as1_predictions)} GATA2-AS1 variants")
            
            # Populate GATA2-AS1 predictions database
            self.db_manager.populate_database(
                'gata2_as1_predictions_db',
                self.gata2_as1_predictions,
                table_name='predictions'
            )
            
            # Analyze EMT pathway genes
            self.emt_scores = variant_scoring.analyze_emt_genes(
                self.all_1000_genomes_variants,
                self.emt_genes,
                evo2_path=config.TOOL_PATHS.get('evo2_executable', '')
            )
            logger.info(f"Analyzed {len(self.emt_scores)} variants in EMT genes")
            
            # Populate EMT variants in functional variants database
            self.db_manager.populate_database(
                'functional_variants_db',
                self.emt_scores,
                table_name='emt_variants'
            )
            
            # Analyze specific lncRNAs (e.g., ANRIL, PIRAT, LRAC, GATA2-AS1)
            target_lncrnas = ['ANRIL', 'PIRAT', 'LRAC', 'GATA2-AS1']
            self.specific_lncrna_analysis = lncrna_scoring.analyze_specific_lncrnas(
                self.lncrna_list,
                target_lncrnas,
                evo2_path=config.TOOL_PATHS.get('evo2_executable', '')
            )
            logger.info(f"Completed specific analysis for {len(self.specific_lncrna_analysis)} lncRNAs")
            
            logger.info("Specific analyses completed successfully")
            return True
            
        except Exception as e:
            logger.error(f"Error in specific analyses step: {e}", exc_info=True)
            return False


def run_evo2_pipeline(config_file: Optional[str] = None) -> bool:
    """
    Wrapper function to execute the complete Evo2 computational pipeline.
    
    Args:
        config_file: Optional path to a configuration file.
    
    Returns:
        bool: True if successful, False if there were errors.
    """
    pipeline = Pipeline(config_file)
    return pipeline.run_pipeline()