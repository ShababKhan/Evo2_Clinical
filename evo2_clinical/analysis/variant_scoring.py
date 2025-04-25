"""
Variant scoring module for the Evo2_Clinical package.

This module provides functions for scoring genetic variants for functional impact
using the Evo2 tool.
"""

import os
import pandas as pd
import numpy as np
import subprocess
import logging
import tempfile
import json
from typing import List, Dict, Union, Optional, Tuple, Any

from ..config import config

logger = logging.getLogger(__name__)


class Evo2Runner:
    """
    Class for interfacing with the Evo2 tool to score genetic variants.
    """
    
    def __init__(self, evo2_path: Optional[str] = None, options: Optional[Dict[str, Any]] = None):
        """
        Initialize the Evo2Runner with path to the Evo2 executable and options.
        
        Args:
            evo2_path: Path to the Evo2 executable or interface.
            options: Optional dictionary of additional options for Evo2.
        """
        self.evo2_path = evo2_path or config.TOOL_PATHS.get('evo2_executable', '')
        self.options = options or {}
        
        if not self.evo2_path:
            logger.warning("Evo2 executable path not provided. Set it with evo2_path parameter or in config.")
    
    def score_variants(self, variants_df: pd.DataFrame, context_window: int = 1000000) -> pd.DataFrame:
        """
        Score genetic variants for functional impact using Evo2.
        
        Args:
            variants_df: DataFrame of genetic variants to score.
            context_window: Size of context window in base pairs (default: 1,000,000).
        
        Returns:
            DataFrame with added functional scores for each variant.
            
        Raises:
            FileNotFoundError: If the Evo2 executable is not found.
            RuntimeError: If Evo2 execution fails.
        """
        if not self.evo2_path:
            raise FileNotFoundError("Evo2 executable path not set")
        
        if not os.path.exists(self.evo2_path):
            raise FileNotFoundError(f"Evo2 executable not found at {self.evo2_path}")
        
        logger.info(f"Running Evo2 for variant scoring with context window size: {context_window}")
        
        # Create a copy of the dataframe to avoid modifying the original
        scored_variants_df = variants_df.copy()
        
        try:
            # Create a temporary file for the input variants
            with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as temp_input:
                temp_input_path = temp_input.name
                # Write variants to the temporary file in VCF format
                _write_vcf_format(scored_variants_df, temp_input)
            
            # Create a temporary file for the output
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as temp_output:
                temp_output_path = temp_output.name
            
            # Build the command to run Evo2
            cmd = [
                self.evo2_path,
                'score',
                '--input', temp_input_path,
                '--output', temp_output_path,
                '--context-window', str(context_window)
            ]
            
            # Add any additional options
            for opt_name, opt_value in self.options.items():
                cmd.extend([f'--{opt_name}', str(opt_value)])
            
            # Execute the command
            logger.debug(f"Executing command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Log any warning messages
            if result.stderr:
                logger.warning(f"Evo2 warnings: {result.stderr}")
            
            # Read the output
            with open(temp_output_path, 'r') as f:
                scores_data = json.load(f)
            
            # Merge scores with the input variants
            scored_variants_df = _merge_scores_with_variants(scored_variants_df, scores_data)
            
            logger.info("Evo2 variant scoring completed successfully")
            return scored_variants_df
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Evo2 execution failed: {e.stderr}")
            raise RuntimeError(f"Evo2 execution failed: {e.stderr}")
        
        except Exception as e:
            logger.error(f"Error during variant scoring: {str(e)}")
            raise
        
        finally:
            # Clean up temporary files
            for path in [temp_input_path, temp_output_path]:
                if os.path.exists(path):
                    os.remove(path)

    def predict_variant_effects(self, variants_df: pd.DataFrame, target_gene: str) -> pd.DataFrame:
        """
        Predict effects of variants on a specific gene using Evo2.
        
        Args:
            variants_df: DataFrame of genetic variants to predict effects for.
            target_gene: Target gene name to analyze effects on.
        
        Returns:
            DataFrame with predicted effects for each variant.
            
        Raises:
            FileNotFoundError: If the Evo2 executable is not found.
            RuntimeError: If Evo2 execution fails.
        """
        if not self.evo2_path:
            raise FileNotFoundError("Evo2 executable path not set")
        
        logger.info(f"Predicting variant effects on gene {target_gene}")
        
        # Create a copy of the dataframe to avoid modifying the original
        effects_df = variants_df.copy()
        
        try:
            # Create a temporary file for the input variants
            with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as temp_input:
                temp_input_path = temp_input.name
                # Write variants to the temporary file in VCF format
                _write_vcf_format(effects_df, temp_input)
            
            # Create a temporary file for the output
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as temp_output:
                temp_output_path = temp_output.name
            
            # Build the command to run Evo2
            cmd = [
                self.evo2_path,
                'predict-effects',
                '--input', temp_input_path,
                '--output', temp_output_path,
                '--target-gene', target_gene
            ]
            
            # Add any additional options
            for opt_name, opt_value in self.options.items():
                cmd.extend([f'--{opt_name}', str(opt_value)])
            
            # Execute the command
            logger.debug(f"Executing command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Log any warning messages
            if result.stderr:
                logger.warning(f"Evo2 warnings: {result.stderr}")
            
            # Read the output
            with open(temp_output_path, 'r') as f:
                effects_data = json.load(f)
            
            # Merge effects predictions with the input variants
            effects_df = _merge_effects_with_variants(effects_df, effects_data)
            
            logger.info(f"Evo2 variant effect prediction for {target_gene} completed successfully")
            return effects_df
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Evo2 execution failed: {e.stderr}")
            raise RuntimeError(f"Evo2 execution failed: {e.stderr}")
        
        except Exception as e:
            logger.error(f"Error during variant effect prediction: {str(e)}")
            raise
        
        finally:
            # Clean up temporary files
            for path in [temp_input_path, temp_output_path]:
                if os.path.exists(path):
                    os.remove(path)


def run_evo2_variant_scoring(variants_df: pd.DataFrame, evo2_path: Optional[str] = None,
                            context_window: int = 1000000) -> pd.DataFrame:
    """
    Interface function to score genetic variants for functional impact using Evo2.
    
    Args:
        variants_df: DataFrame of genetic variants to score.
        evo2_path: Path to the Evo2 executable or interface (if None, uses config value).
        context_window: Size of context window in base pairs (default: 1,000,000).
    
    Returns:
        DataFrame with added functional scores for each variant.
    """
    evo2_runner = Evo2Runner(evo2_path=evo2_path)
    return evo2_runner.score_variants(variants_df, context_window=context_window)


def analyze_emt_genes(variants_df: pd.DataFrame, emt_genes_list: List[str], 
                     evo2_path: Optional[str] = None) -> pd.DataFrame:
    """
    Analyzes genetic variants in the Endothelial-Mesenchymal Transition (EMT) pathway
    using Evo2.
    
    Args:
        variants_df: DataFrame of variants to analyze.
        emt_genes_list: List of EMT gene names.
        evo2_path: Path to the Evo2 executable or interface (if None, uses config value).
    
    Returns:
        DataFrame with scores/analysis results for EMT gene variants.
    """
    logger.info(f"Analyzing genetic variants in {len(emt_genes_list)} EMT genes")
    
    # Filter variants to only include those in EMT genes
    # This is a placeholder implementation - actual filtering would depend on
    # how genes are identified in the variants dataframe
    emt_variants = variants_df[variants_df['INFO'].apply(
        lambda x: any(gene in x for gene in emt_genes_list)
    )]
    
    if len(emt_variants) == 0:
        logger.warning("No variants found in EMT genes")
        return pd.DataFrame(columns=variants_df.columns + ['functional_score'])
    
    logger.info(f"Found {len(emt_variants)} variants in EMT genes")
    
    # Score the variants
    emt_scores_df = run_evo2_variant_scoring(emt_variants, evo2_path)
    
    logger.info("EMT gene analysis completed")
    return emt_scores_df


def predict_gata2_as1_variant_effects(variants_df: pd.DataFrame, evo2_path: Optional[str] = None) -> pd.DataFrame:
    """
    Performs specific variant effect predictions for the lncRNA GATA2-AS1 using Evo2.
    
    Args:
        variants_df: DataFrame of variants for analysis.
        evo2_path: Path to the Evo2 executable or interface (if None, uses config value).
    
    Returns:
        DataFrame with predicted effects for GATA2-AS1 variants.
    """
    logger.info("Predicting variant effects for GATA2-AS1")
    
    evo2_runner = Evo2Runner(evo2_path=evo2_path)
    gata2_as1_predictions = evo2_runner.predict_variant_effects(variants_df, target_gene="GATA2-AS1")
    
    logger.info("GATA2-AS1 variant effect prediction completed")
    return gata2_as1_predictions


def analyze_epigenetic_mediators(variants_df: pd.DataFrame, epigenetic_genes: List[str], 
                               evo2_path: Optional[str] = None) -> pd.DataFrame:
    """
    Analyzes epigenetic mediators with functionally relevant SNPs using Evo2.
    
    Args:
        variants_df: DataFrame of variants for analysis.
        epigenetic_genes: List of epigenetic mediator gene names.
        evo2_path: Path to the Evo2 executable or interface (if None, uses config value).
    
    Returns:
        DataFrame with scores for epigenetic mediator variants.
    """
    logger.info(f"Analyzing {len(epigenetic_genes)} epigenetic mediators with functional SNPs")
    
    # Filter variants to only include those in epigenetic mediator genes
    # This is a placeholder implementation
    epigenetic_variants = variants_df[variants_df['INFO'].apply(
        lambda x: any(gene in x for gene in epigenetic_genes)
    )]
    
    if len(epigenetic_variants) == 0:
        logger.warning("No variants found in epigenetic mediator genes")
        return pd.DataFrame(columns=variants_df.columns + ['functional_score'])
    
    logger.info(f"Found {len(epigenetic_variants)} variants in epigenetic mediator genes")
    
    # Score the variants
    epigenetic_scores_df = run_evo2_variant_scoring(epigenetic_variants, evo2_path)
    
    logger.info("Epigenetic mediator analysis completed")
    return epigenetic_scores_df


def run_swap_snp_pipeline(human_variants_df: pd.DataFrame, animal_model_data: pd.DataFrame,
                        evo2_path: Optional[str] = None, context_window: int = 1000000) -> pd.DataFrame:
    """
    Runs the "swap SNP pipeline" with Evo2 for translational bridging between human and animal model data.
    
    Args:
        human_variants_df: Human variant data.
        animal_model_data: Variant data from the animal model.
        evo2_path: Path to the Evo2 executable or interface (if None, uses config value).
        context_window: Size of context window in base pairs (default: 1,000,000).
    
    Returns:
        DataFrame with results of the translational bridging analysis.
    """
    logger.info("Running swap SNP pipeline for translational bridging")
    
    evo2_path = evo2_path or config.TOOL_PATHS.get('evo2_executable', '')
    
    if not evo2_path:
        raise FileNotFoundError("Evo2 executable path not set")
    
    # This is a placeholder implementation - actual implementation would depend on
    # the specific Evo2 swap SNP pipeline functionality
    
    # For demonstration, we'll combine human and animal variants and score them
    combined_df = pd.concat([
        human_variants_df.assign(source='human'),
        animal_model_data.assign(source='animal')
    ]).reset_index(drop=True)
    
    # Score all variants
    scored_variants = run_evo2_variant_scoring(combined_df, evo2_path, context_window)
    
    # Create a result dataframe with comparative analysis
    # In a real implementation, this would involve more sophisticated analysis
    translational_results_df = scored_variants.copy()
    
    logger.info("Swap SNP pipeline completed")
    return translational_results_df


# --- Helper functions ---

def _write_vcf_format(variants_df: pd.DataFrame, file_handle) -> None:
    """Write variants DataFrame to a file handle in VCF format."""
    # Write VCF header
    file_handle.write("##fileformat=VCFv4.2\n")
    file_handle.write("##source=Evo2_Clinical\n")
    
    # Write column headers
    file_handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    
    # Write data rows
    for _, row in variants_df.iterrows():
        chrom = row.get('#CHROM', row.get('CHROM', '.'))
        pos = row.get('POS', '.')
        id_val = row.get('ID', '.')
        ref = row.get('REF', '.')
        alt = row.get('ALT', '.')
        qual = row.get('QUAL', '.')
        filter_val = row.get('FILTER', '.')
        info = row.get('INFO', '.')
        
        file_handle.write(f"{chrom}\t{pos}\t{id_val}\t{ref}\t{alt}\t{qual}\t{filter_val}\t{info}\n")


def _merge_scores_with_variants(variants_df: pd.DataFrame, scores_data: Dict) -> pd.DataFrame:
    """Merge Evo2 scores with input variants DataFrame."""
    # This is a placeholder implementation - actual implementation would depend on
    # the specific output format of the Evo2 tool
    
    result_df = variants_df.copy()
    
    # Assuming scores_data is a dictionary with variant IDs as keys
    # and scores as values
    score_dict = {}
    
    # Extract scores from the scores_data
    for variant_id, data in scores_data.get('variants', {}).items():
        score_dict[variant_id] = data.get('functional_score', 0.0)
    
    # Create a functional score column
    def get_score(row):
        # Try to match by ID if available
        if row.get('ID', '.') != '.' and row.get('ID') in score_dict:
            return score_dict[row['ID']]
        
        # Otherwise, try to create a key from position info
        key = f"{row.get('#CHROM', row.get('CHROM', '.'))}:{row.get('POS')}_{row.get('REF')}_{row.get('ALT')}"
        return score_dict.get(key, np.nan)
    
    result_df['functional_score'] = result_df.apply(get_score, axis=1)
    
    # Fill any missing scores with NaN
    result_df['functional_score'] = pd.to_numeric(result_df['functional_score'], errors='coerce')
    
    return result_df


def _merge_effects_with_variants(variants_df: pd.DataFrame, effects_data: Dict) -> pd.DataFrame:
    """Merge Evo2 effect predictions with input variants DataFrame."""
    # This is a placeholder implementation - actual implementation would depend on
    # the specific output format of the Evo2 tool
    
    result_df = variants_df.copy()
    
    # Assuming effects_data contains variant IDs and their predicted effects
    effect_dict = {}
    confidence_dict = {}
    
    # Extract effects from the effects_data
    for variant_id, data in effects_data.get('variants', {}).items():
        effect_dict[variant_id] = data.get('predicted_effect', 'unknown')
        confidence_dict[variant_id] = data.get('confidence', 0.0)
    
    # Create effect and confidence columns
    def get_effect_data(row):
        # Try to match by ID if available
        variant_id = row.get('ID', '.')
        if variant_id != '.' and variant_id in effect_dict:
            return effect_dict[variant_id], confidence_dict[variant_id]
        
        # Otherwise, try to create a key from position info
        key = f"{row.get('#CHROM', row.get('CHROM', '.'))}:{row.get('POS')}_{row.get('REF')}_{row.get('ALT')}"
        return effect_dict.get(key, 'unknown'), confidence_dict.get(key, 0.0)
    
    # Apply the function and create new columns
    effects_and_confidence = result_df.apply(get_effect_data, axis=1, result_type='expand')
    result_df['predicted_effect'] = effects_and_confidence[0]
    result_df['confidence_score'] = pd.to_numeric(effects_and_confidence[1], errors='coerce')
    
    return result_df