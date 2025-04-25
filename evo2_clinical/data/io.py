"""
Data loading and processing module for the Evo2_Clinical package.

This module provides functions for loading and processing genomic data from various sources,
including GWAS catalogs, 1000 Genomes, ENCODE data, and custom gene lists.
"""

import os
import pandas as pd
import numpy as np
import gzip
import logging
from typing import List, Dict, Union, Optional, Tuple


logger = logging.getLogger(__name__)


def load_gwas_catalog_data(filepath: str, diseases: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Loads and preprocesses data from a GWAS catalog file.
    Filters for relevant diseases if provided.

    Args:
        filepath: Path to the GWAS catalog file.
        diseases: Optional list of disease names to filter for.

    Returns:
        DataFrame containing relevant GWAS data.
        
    Raises:
        FileNotFoundError: If the file doesn't exist.
        ValueError: If the file format is incorrect.
    """
    logger.info(f"Loading GWAS catalog data from {filepath}")
    
    try:
        # Assuming tab-separated file
        df = pd.read_csv(filepath, sep='\t')
        
        # Filter by disease if specified
        if diseases:
            disease_pattern = '|'.join(diseases)
            df = df[df['DISEASE/TRAIT'].str.contains(disease_pattern, case=False, na=False)]
            logger.info(f"Filtered GWAS data for diseases: {', '.join(diseases)}")
            
        logger.info(f"Loaded {len(df)} GWAS catalog entries")
        return df
    
    except FileNotFoundError:
        logger.error(f"GWAS catalog file not found: {filepath}")
        raise
    except Exception as e:
        logger.error(f"Error loading GWAS catalog data: {e}")
        raise ValueError(f"Error processing GWAS catalog file: {str(e)}")


def load_vcf_file(filepath: str, region: Optional[str] = None) -> pd.DataFrame:
    """
    Loads a VCF file into a pandas DataFrame.
    
    Args:
        filepath: Path to the VCF file (can be gzipped).
        region: Optional genomic region to filter (format: 'chr:start-end').
        
    Returns:
        DataFrame containing VCF data with standard VCF columns.
        
    Raises:
        FileNotFoundError: If the file doesn't exist.
        ValueError: If the file format is incorrect.
    """
    logger.info(f"Loading VCF data from {filepath}")
    
    try:
        # Check if file is gzipped
        open_func = gzip.open if filepath.endswith('.gz') else open
        
        # Parse VCF file
        header_lines = []
        header = None
        rows = []
        
        with open_func(filepath, 'rt') as f:
            for line in f:
                if line.startswith('##'):
                    header_lines.append(line.strip())
                elif line.startswith('#'):
                    header = line[1:].strip().split('\t')
                elif region is None:
                    rows.append(line.strip().split('\t'))
                elif region and _is_in_region(line, region):
                    rows.append(line.strip().split('\t'))
        
        if not header:
            raise ValueError("Invalid VCF format: no column header found")
        
        # Create DataFrame
        df = pd.DataFrame(rows, columns=header)
        logger.info(f"Loaded {len(df)} variants from VCF file")
        return df
    
    except FileNotFoundError:
        logger.error(f"VCF file not found: {filepath}")
        raise
    except Exception as e:
        logger.error(f"Error loading VCF file: {e}")
        raise ValueError(f"Error processing VCF file: {str(e)}")


def load_1000_genomes_data(filepath: str, endothelial_genes: Union[List[str], str], 
                           gene_coords_file: Optional[str] = None) -> pd.DataFrame:
    """
    Loads and processes 1000 Genomes VCF data.
    Filters variants to include only those within or near specified endothelial genes.
    Identifies common and rare variants.

    Args:
        filepath: Path to the 1000 Genomes VCF file.
        endothelial_genes: Either a list of endothelial gene names or a path to a file containing them.
        gene_coords_file: Optional path to a file with gene coordinates for filtering by position.

    Returns:
        DataFrame containing filtered 1000 Genomes variant data.
        
    Raises:
        FileNotFoundError: If any of the files don't exist.
        ValueError: If the file format is incorrect.
    """
    logger.info(f"Loading 1000 Genomes data and filtering by endothelial genes")
    
    # Load gene list if a string (assumed to be a file path)
    gene_list = _load_gene_list(endothelial_genes) if isinstance(endothelial_genes, str) else endothelial_genes
    
    # Load gene coordinates if provided
    gene_coords = None
    if gene_coords_file:
        logger.info(f"Loading gene coordinates from {gene_coords_file}")
        gene_coords = pd.read_csv(gene_coords_file, sep='\t')
    
    # Load 1000 Genomes data
    variants_df = load_vcf_file(filepath)
    
    # Filter variants by gene
    if gene_coords is not None:
        # Filter by coordinates using gene_coords dataframe
        filtered_df = _filter_variants_by_coords(variants_df, gene_coords, gene_list)
    else:
        # Basic filtering - assuming INFO field contains gene annotations
        filtered_df = variants_df[variants_df['INFO'].apply(lambda x: any(gene in x for gene in gene_list))]
    
    # Add common/rare variant annotations
    filtered_df = _annotate_variant_frequency(filtered_df)
    
    logger.info(f"Filtered to {len(filtered_df)} variants in endothelial genes")
    return filtered_df


def filter_by_encode_data(variants_df: pd.DataFrame, encode_filepath: str, 
                          cell_type: str = "ENDOS") -> pd.DataFrame:
    """
    Filters genetic variants based on ENCODE data, specifically for activity
    in the specified cell type (default: ENDOS - endothelial cells).

    Args:
        variants_df: DataFrame of genetic variants.
        encode_filepath: Path to the ENCODE data file or API endpoint.
        cell_type: The cell type to filter for (e.g., "ENDOS").

    Returns:
        Filtered DataFrame containing variants active in the specified cell type.
        
    Raises:
        FileNotFoundError: If the ENCODE file doesn't exist.
        ValueError: If there are issues with the file format.
    """
    logger.info(f"Filtering variants based on ENCODE data for {cell_type}")
    
    try:
        # Load ENCODE regulatory elements for the specified cell type
        encode_data = _load_encode_regulatory_elements(encode_filepath, cell_type)
        
        # Filter variants that overlap with regulatory elements
        filtered_df = _filter_variants_by_regulatory_elements(variants_df, encode_data)
        
        logger.info(f"Filtered to {len(filtered_df)} variants active in {cell_type} cells")
        return filtered_df
    
    except FileNotFoundError:
        logger.error(f"ENCODE data file not found: {encode_filepath}")
        raise
    except Exception as e:
        logger.error(f"Error filtering by ENCODE data: {e}")
        raise ValueError(f"Error processing ENCODE data: {str(e)}")


def load_gene_list(filepath: str) -> List[str]:
    """
    Loads a list of gene names from a file.

    Args:
        filepath: Path to the file containing gene names (one per line).

    Returns:
        List of gene names.
        
    Raises:
        FileNotFoundError: If the file doesn't exist.
    """
    return _load_gene_list(filepath)


def load_lncrna_list(filepath: str) -> List[str]:
    """
    Loads a list of lncRNA names from a file.

    Args:
        filepath: Path to the file containing lncRNA names (one per line).

    Returns:
        List of lncRNA names.
        
    Raises:
        FileNotFoundError: If the file doesn't exist.
    """
    logger.info(f"Loading lncRNA list from {filepath}")
    return _load_gene_list(filepath)


def load_emt_genes(filepath: str) -> List[str]:
    """
    Loads a list of EMT gene names from a file.

    Args:
        filepath: Path to the file containing EMT gene names (one per line).

    Returns:
        List of EMT gene names.
        
    Raises:
        FileNotFoundError: If the file doesn't exist.
    """
    logger.info(f"Loading EMT gene list from {filepath}")
    return _load_gene_list(filepath)


def load_gtf_file(filepath: str) -> pd.DataFrame:
    """
    Loads a GTF file into a pandas DataFrame.
    
    Args:
        filepath: Path to the GTF file (can be gzipped).
        
    Returns:
        DataFrame containing GTF data with standard columns.
        
    Raises:
        FileNotFoundError: If the file doesn't exist.
        ValueError: If the file format is incorrect.
    """
    logger.info(f"Loading GTF data from {filepath}")
    
    try:
        # Check if file is gzipped
        open_func = gzip.open if filepath.endswith('.gz') else open
        
        # Parse GTF file
        columns = ['seqname', 'source', 'feature', 'start', 'end', 
                   'score', 'strand', 'frame', 'attribute']
        rows = []
        
        with open_func(filepath, 'rt') as f:
            for line in f:
                if not line.startswith('#'):
                    rows.append(line.strip().split('\t'))
        
        # Create DataFrame
        df = pd.DataFrame(rows, columns=columns)
        
        # Parse attributes
        df['gene_id'] = df['attribute'].apply(lambda x: _extract_gtf_attribute(x, 'gene_id'))
        df['gene_name'] = df['attribute'].apply(lambda x: _extract_gtf_attribute(x, 'gene_name'))
        
        logger.info(f"Loaded {len(df)} entries from GTF file")
        return df
    
    except FileNotFoundError:
        logger.error(f"GTF file not found: {filepath}")
        raise
    except Exception as e:
        logger.error(f"Error loading GTF file: {e}")
        raise ValueError(f"Error processing GTF file: {str(e)}")


# --- Helper functions ---

def _load_gene_list(filepath: str) -> List[str]:
    """Internal function to load a list of genes from a file."""
    try:
        with open(filepath, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
        logger.info(f"Loaded {len(genes)} entries from {filepath}")
        return genes
    except FileNotFoundError:
        logger.error(f"File not found: {filepath}")
        raise


def _is_in_region(vcf_line: str, region: str) -> bool:
    """Check if a VCF line is in the specified genomic region."""
    if not region or not vcf_line:
        return False
    
    parts = vcf_line.split('\t')
    if not parts:
        return False
    
    chrom = parts[0]
    pos = int(parts[1])
    
    # Parse region string (format: 'chr:start-end')
    try:
        region_chrom, pos_range = region.split(':')
        region_start, region_end = map(int, pos_range.split('-'))
        
        return chrom == region_chrom and region_start <= pos <= region_end
    except:
        return False


def _filter_variants_by_coords(variants_df: pd.DataFrame, gene_coords: pd.DataFrame, 
                              gene_list: List[str]) -> pd.DataFrame:
    """Filter variants by gene coordinates."""
    # Filter gene coordinates to only include genes in our list
    gene_coords = gene_coords[gene_coords['gene_name'].isin(gene_list)]
    
    filtered_variants = []
    
    # For each gene, get variants that fall within its coordinates
    for _, gene_row in gene_coords.iterrows():
        chrom = gene_row['chrom']
        start = gene_row['start']
        end = gene_row['end']
        
        # Filter variants
        gene_variants = variants_df[
            (variants_df['#CHROM'] == chrom) & 
            (variants_df['POS'].astype(int) >= start) & 
            (variants_df['POS'].astype(int) <= end)
        ]
        
        filtered_variants.append(gene_variants)
    
    # Combine all filtered variants
    if filtered_variants:
        return pd.concat(filtered_variants).drop_duplicates()
    else:
        return pd.DataFrame(columns=variants_df.columns)


def _annotate_variant_frequency(variants_df: pd.DataFrame) -> pd.DataFrame:
    """Add common/rare annotations based on allele frequency."""
    variants_df = variants_df.copy()
    
    # Parse INFO field to extract allele frequency
    variants_df['AF'] = variants_df['INFO'].apply(_extract_allele_frequency)
    
    # Define common vs. rare based on frequency
    variants_df['is_common'] = variants_df['AF'] >= 0.05
    variants_df['is_rare'] = variants_df['AF'] < 0.01
    
    return variants_df


def _extract_allele_frequency(info_field: str) -> float:
    """Extract allele frequency from VCF INFO field."""
    try:
        # This assumes AF is in the INFO field in standard format
        if 'AF=' in info_field:
            af_part = info_field.split('AF=')[1].split(';')[0]
            return float(af_part)
        return 0.0  # Default if not found
    except:
        return 0.0


def _load_encode_regulatory_elements(filepath: str, cell_type: str) -> pd.DataFrame:
    """Load ENCODE regulatory elements for a specific cell type."""
    # This is a placeholder function that would be implemented based on the actual ENCODE data format
    logger.info(f"Loading ENCODE regulatory elements for {cell_type}")
    
    try:
        # Assuming a BED-like format with cell type information
        encode_df = pd.read_csv(filepath, sep='\t')
        
        # Filter for the specified cell type
        filtered_df = encode_df[encode_df['cell_type'] == cell_type]
        
        logger.info(f"Loaded {len(filtered_df)} regulatory elements for {cell_type}")
        return filtered_df
    
    except Exception as e:
        logger.error(f"Error loading ENCODE data: {e}")
        raise


def _filter_variants_by_regulatory_elements(variants_df: pd.DataFrame, 
                                           encode_df: pd.DataFrame) -> pd.DataFrame:
    """Filter variants that overlap with regulatory elements."""
    # This is a placeholder function that would be implemented based on the actual data structure
    filtered_variants = []
    
    # For each regulatory element, find overlapping variants
    for _, reg_element in encode_df.iterrows():
        chrom = reg_element['chrom']
        start = reg_element['start']
        end = reg_element['end']
        
        # Filter variants
        overlapping_variants = variants_df[
            (variants_df['#CHROM'] == chrom) & 
            (variants_df['POS'].astype(int) >= start) & 
            (variants_df['POS'].astype(int) <= end)
        ]
        
        filtered_variants.append(overlapping_variants)
    
    # Combine all filtered variants
    if filtered_variants:
        return pd.concat(filtered_variants).drop_duplicates()
    else:
        return pd.DataFrame(columns=variants_df.columns)


def _extract_gtf_attribute(attr_str: str, key: str) -> str:
    """Extract a specific attribute from a GTF attribute string."""
    try:
        if f'{key} "' in attr_str:
            value = attr_str.split(f'{key} "')[1].split('";')[0]
            return value
        return ""
    except:
        return ""