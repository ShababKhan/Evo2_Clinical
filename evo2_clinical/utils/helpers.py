"""
Utility functions for the Evo2_Clinical package.

This module provides helper functions used across the package.
"""

import os
import time
import logging
import pandas as pd
import numpy as np
from typing import List, Dict, Optional, Any, Union
import yaml
from datetime import datetime

logger = logging.getLogger(__name__)


def setup_logging(log_dir: str, log_level: int = logging.INFO) -> str:
    """
    Set up logging with file and console handlers.
    
    Args:
        log_dir: Directory where log files will be stored.
        log_level: Logging level (default: logging.INFO).
    
    Returns:
        str: Path to the created log file.
    """
    os.makedirs(log_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"evo2_clinical_{timestamp}.log")
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    
    # Remove any existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(log_level)
    file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)
    root_logger.addHandler(file_handler)
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_formatter = logging.Formatter('%(levelname)s - %(message)s')
    console_handler.setFormatter(console_formatter)
    root_logger.addHandler(console_handler)
    
    logger.info(f"Logging setup complete. Log file: {log_file}")
    return log_file


def load_yaml_config(config_file: str) -> Dict[str, Any]:
    """
    Load a YAML configuration file.
    
    Args:
        config_file: Path to the YAML configuration file.
    
    Returns:
        Dictionary containing the configuration.
        
    Raises:
        FileNotFoundError: If the file doesn't exist.
        ValueError: If the file is not a valid YAML file.
    """
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        return config
    except FileNotFoundError:
        logger.error(f"Configuration file not found: {config_file}")
        raise
    except yaml.YAMLError as e:
        logger.error(f"Error parsing YAML configuration: {e}")
        raise ValueError(f"Invalid YAML in configuration file: {e}")


def save_yaml_config(config: Dict[str, Any], config_file: str) -> None:
    """
    Save a dictionary to a YAML configuration file.
    
    Args:
        config: Dictionary to save.
        config_file: Path where to save the YAML file.
        
    Raises:
        IOError: If the file can't be written.
    """
    try:
        os.makedirs(os.path.dirname(os.path.abspath(config_file)), exist_ok=True)
        
        with open(config_file, 'w') as f:
            yaml.dump(config, f, default_flow_style=False)
            
        logger.info(f"Configuration saved to {config_file}")
    except Exception as e:
        logger.error(f"Error saving configuration to {config_file}: {e}")
        raise


def time_function(func):
    """
    Decorator to time a function's execution.
    
    Args:
        func: The function to time.
        
    Returns:
        Wrapped function that logs execution time.
    """
    def wrapper(*args, **kwargs):
        start_time = time.time()
        logger.info(f"Starting {func.__name__}")
        
        result = func(*args, **kwargs)
        
        end_time = time.time()
        execution_time = end_time - start_time
        logger.info(f"Completed {func.__name__} in {execution_time:.2f} seconds")
        
        return result
    
    return wrapper


def parse_variant_id(variant_id: str) -> Dict[str, str]:
    """
    Parse a variant ID string into its components.
    
    Args:
        variant_id: Variant ID string (format: "chr:pos_ref_alt").
        
    Returns:
        Dictionary with the components: {'chrom': str, 'pos': int, 'ref': str, 'alt': str}.
        
    Raises:
        ValueError: If the variant ID format is invalid.
    """
    try:
        # Split into chromosome:position and ref_alt parts
        location, alleles = variant_id.split('_', 1)
        chrom, pos = location.split(':', 1)
        ref, alt = alleles.split('_', 1)
        
        return {
            'chrom': chrom,
            'pos': int(pos),
            'ref': ref,
            'alt': alt
        }
    except Exception as e:
        raise ValueError(f"Invalid variant ID format '{variant_id}'. Expected 'chr:pos_ref_alt': {e}")


def create_variant_id(chrom: str, pos: Union[int, str], ref: str, alt: str) -> str:
    """
    Create a variant ID string from its components.
    
    Args:
        chrom: Chromosome name.
        pos: Position.
        ref: Reference allele.
        alt: Alternative allele.
        
    Returns:
        Variant ID string (format: "chr:pos_ref_alt").
    """
    return f"{chrom}:{pos}_{ref}_{alt}"


def extract_gene_name_from_vcf_info(info_field: str) -> Optional[str]:
    """
    Extract gene name from a VCF INFO field.
    
    Args:
        info_field: VCF INFO field string.
        
    Returns:
        Gene name if found, None otherwise.
    """
    if 'GENE=' in info_field:
        try:
            gene_part = info_field.split('GENE=')[1].split(';')[0]
            return gene_part
        except:
            pass
    
    return None


def summarize_dataframe(df: pd.DataFrame, title: str = "DataFrame Summary") -> str:
    """
    Create a string summary of a DataFrame.
    
    Args:
        df: DataFrame to summarize.
        title: Title for the summary.
        
    Returns:
        String containing the summary information.
    """
    if df is None or len(df) == 0:
        return f"{title}: Empty DataFrame"
    
    summary = [
        f"{title}:",
        f"  - Shape: {df.shape[0]} rows x {df.shape[1]} columns",
        f"  - Columns: {', '.join(df.columns)}",
        "  - Data types:"
    ]
    
    for col, dtype in df.dtypes.items():
        summary.append(f"    * {col}: {dtype}")
    
    summary.append("  - Missing values:")
    for col, missing in df.isna().sum().items():
        if missing > 0:
            summary.append(f"    * {col}: {missing} ({missing/len(df)*100:.1f}%)")
    
    return '\n'.join(summary)