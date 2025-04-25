"""
lncRNA scoring module for the Evo2_Clinical package.

This module provides functions for scoring lncRNA functionality using the Evo2 tool.
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


class lncRNAScorer:
    """
    Class for scoring lncRNA functionality using the Evo2 tool.
    """
    
    def __init__(self, evo2_path: Optional[str] = None, options: Optional[Dict[str, Any]] = None):
        """
        Initialize the lncRNAScorer with path to the Evo2 executable and options.
        
        Args:
            evo2_path: Path to the Evo2 executable or interface.
            options: Optional dictionary of additional options for Evo2.
        """
        self.evo2_path = evo2_path or config.TOOL_PATHS.get('evo2_executable', '')
        self.options = options or {}
        
        if not self.evo2_path:
            logger.warning("Evo2 executable path not provided. Set it with evo2_path parameter or in config.")
    
    def score_lncrnas(self, lncrnas: Union[List[str], pd.DataFrame]) -> pd.DataFrame:
        """
        Score lncRNAs for functionality using Evo2.
        
        Args:
            lncrnas: Either a list of lncRNA names or a DataFrame containing lncRNA data.
                     If DataFrame, must contain 'lncrna' column with names.
        
        Returns:
            DataFrame with lncRNA names and their functionality scores.
            
        Raises:
            FileNotFoundError: If the Evo2 executable is not found.
            RuntimeError: If Evo2 execution fails.
            ValueError: If invalid input format is provided.
        """
        if not self.evo2_path:
            raise FileNotFoundError("Evo2 executable path not set")
            
        if not os.path.exists(self.evo2_path):
            raise FileNotFoundError(f"Evo2 executable not found at {self.evo2_path}")
        
        # Convert input to a standard format
        if isinstance(lncrnas, list):
            lncrnas_df = pd.DataFrame({'lncrna': lncrnas})
        elif isinstance(lncrnas, pd.DataFrame):
            if 'lncrna' not in lncrnas.columns:
                raise ValueError("Input DataFrame must contain 'lncrna' column")
            lncrnas_df = lncrnas.copy()
        else:
            raise ValueError("Input must be either a list of lncRNA names or a DataFrame")
        
        logger.info(f"Scoring functionality for {len(lncrnas_df)} lncRNAs")
        
        try:
            # Create a temporary file for the input lncRNAs
            with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_input:
                temp_input_path = temp_input.name
                # Write lncRNA names to the temporary file
                for lncrna in lncrnas_df['lncrna']:
                    temp_input.write(f"{lncrna}\n")
            
            # Create a temporary file for the output
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as temp_output:
                temp_output_path = temp_output.name
            
            # Build the command to run Evo2
            cmd = [
                self.evo2_path,
                'score-lncrna',  # Assuming this is the lncRNA scoring command
                '--input', temp_input_path,
                '--output', temp_output_path
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
            
            # Process the scores
            scored_lncrnas_df = _process_lncrna_scores(lncrnas_df, scores_data)
            
            logger.info("lncRNA scoring completed successfully")
            return scored_lncrnas_df
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Evo2 execution failed: {e.stderr}")
            raise RuntimeError(f"Evo2 execution failed: {e.stderr}")
        
        except Exception as e:
            logger.error(f"Error during lncRNA scoring: {str(e)}")
            raise
        
        finally:
            # Clean up temporary files
            for path in [temp_input_path, temp_output_path]:
                if os.path.exists(path):
                    os.remove(path)
    
    def analyze_lncrna_structure(self, lncrna_name: str, sequence: Optional[str] = None) -> Dict[str, Any]:
        """
        Analyze the structure of a specific lncRNA.
        
        Args:
            lncrna_name: Name of the lncRNA to analyze.
            sequence: Optional RNA sequence. If None, the sequence will be retrieved from databases.
        
        Returns:
            Dictionary containing structural analysis results.
            
        Raises:
            FileNotFoundError: If the Evo2 executable is not found.
            RuntimeError: If Evo2 execution fails.
        """
        if not self.evo2_path:
            raise FileNotFoundError("Evo2 executable path not set")
            
        logger.info(f"Analyzing structure of lncRNA {lncrna_name}")
        
        try:
            # Create a temporary file for the output
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as temp_output:
                temp_output_path = temp_output.name
            
            # Build the command
            cmd = [
                self.evo2_path,
                'analyze-lncrna',  # Assuming this is the lncRNA analysis command
                '--name', lncrna_name,
                '--output', temp_output_path
            ]
            
            # Add sequence if provided
            if sequence:
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_seq:
                    temp_seq_path = temp_seq.name
                    temp_seq.write(f">{lncrna_name}\n{sequence}\n")
                cmd.extend(['--sequence', temp_seq_path])
            
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
                analysis_data = json.load(f)
            
            logger.info(f"lncRNA {lncrna_name} structure analysis completed")
            return analysis_data
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Evo2 execution failed: {e.stderr}")
            raise RuntimeError(f"Evo2 execution failed: {e.stderr}")
        
        except Exception as e:
            logger.error(f"Error during lncRNA structure analysis: {str(e)}")
            raise
        
        finally:
            # Clean up temporary files
            for path in [temp_output_path]:
                if os.path.exists(path):
                    os.remove(path)
            if sequence and 'temp_seq_path' in locals():
                if os.path.exists(temp_seq_path):
                    os.remove(temp_seq_path)


def run_evo2_lncrna_scoring(lncrnas: Union[List[str], pd.DataFrame], 
                           evo2_path: Optional[str] = None) -> pd.DataFrame:
    """
    Interface function to score lncRNAs for functionality using Evo2.
    
    Args:
        lncrnas: Either a list of lncRNA names or a DataFrame containing lncRNA data.
        evo2_path: Path to the Evo2 executable or interface (if None, uses config value).
    
    Returns:
        DataFrame with lncRNA names and their functionality scores.
    """
    scorer = lncRNAScorer(evo2_path=evo2_path)
    return scorer.score_lncrnas(lncrnas)


def analyze_specific_lncrnas(lncrna_list: List[str], target_lncrnas: List[str], 
                            evo2_path: Optional[str] = None) -> pd.DataFrame:
    """
    Analyze specific lncRNAs (e.g., ANRIL, PIRAT, LRAC, GATA2-AS1) mentioned in the project.
    
    Args:
        lncrna_list: List of all lncRNA names.
        target_lncrnas: List of specific lncRNAs to focus on.
        evo2_path: Path to the Evo2 executable or interface (if None, uses config value).
    
    Returns:
        DataFrame with enriched analysis for the target lncRNAs.
    """
    logger.info(f"Performing specialized analysis for specific lncRNAs: {', '.join(target_lncrnas)}")
    
    # Score all lncRNAs first
    all_scores_df = run_evo2_lncrna_scoring(lncrna_list, evo2_path)
    
    # Filter for the target lncRNAs
    target_scores_df = all_scores_df[all_scores_df['lncrna'].isin(target_lncrnas)].copy()
    
    if len(target_scores_df) == 0:
        logger.warning(f"None of the target lncRNAs {target_lncrnas} found in the results")
        return pd.DataFrame(columns=all_scores_df.columns)
    
    # Enrich with additional analyses
    scorer = lncRNAScorer(evo2_path=evo2_path)
    
    # Perform detailed structure analysis for each target lncRNA
    structure_data = []
    for lncrna in target_scores_df['lncrna']:
        try:
            analysis = scorer.analyze_lncrna_structure(lncrna)
            structure_data.append({
                'lncrna': lncrna,
                'structure_score': analysis.get('structure_score', np.nan),
                'secondary_structure': analysis.get('secondary_structure', ''),
                'binding_domains': len(analysis.get('binding_domains', [])),
                'conservation_score': analysis.get('conservation_score', np.nan)
            })
        except Exception as e:
            logger.error(f"Error analyzing structure for {lncrna}: {e}")
            structure_data.append({
                'lncrna': lncrna,
                'structure_score': np.nan,
                'secondary_structure': '',
                'binding_domains': 0,
                'conservation_score': np.nan
            })
    
    # Create DataFrame from structure analysis
    structure_df = pd.DataFrame(structure_data)
    
    # Merge with scores
    enriched_df = pd.merge(target_scores_df, structure_df, on='lncrna', how='left')
    
    logger.info(f"Completed specialized analysis for {len(enriched_df)} lncRNAs")
    return enriched_df


# --- Helper functions ---

def _process_lncrna_scores(lncrnas_df: pd.DataFrame, scores_data: Dict) -> pd.DataFrame:
    """Process and merge lncRNA scores from Evo2 output."""
    # This is a placeholder implementation - actual implementation would depend on
    # the specific output format of the Evo2 tool
    
    result_df = lncrnas_df.copy()
    
    # Check if the scores_data dictionary contains the expected structure
    if 'lncrnas' not in scores_data:
        logger.warning("Unexpected Evo2 output format: 'lncrnas' key missing")
        result_df['functionality_score'] = np.nan
        return result_df
    
    # Create a dictionary mapping lncRNA names to scores
    score_dict = {
        lnc_data.get('name'): lnc_data.get('functionality_score', np.nan)
        for lnc_data in scores_data['lncrnas']
    }
    
    # Add functionality scores to the dataframe
    result_df['functionality_score'] = result_df['lncrna'].map(score_dict)
    result_df['functionality_score'] = pd.to_numeric(result_df['functionality_score'], errors='coerce')
    
    # Check for missing scores
    missing_scores = result_df['functionality_score'].isna().sum()
    if missing_scores > 0:
        logger.warning(f"{missing_scores} lncRNAs have missing scores")
    
    return result_df