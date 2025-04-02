#!/usr/bin/env python3
"""
Variant Analysis Module for Evo2 Pipeline

This module provides utilities for analyzing genomic variants, including:
- Processing VCF files from 1000 Genomes and custom sources
- Analyzing variants in any gene or genomic region
- Scoring variant effects using various computational methods
- Integrating ENCODE data for cell-type specific analysis

Author: Shabab Khan
Date: March 26, 2025
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Union
import logging
from pathlib import Path
import cyvcf2
import allel
import pysam
from Bio import SeqIO
from pybedtools import BedTool

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("VariantAnalysis")

class VariantProcessor:
    """Class for processing and analyzing genomic variants."""
    
    def __init__(self, reference_genome: str = "GRCh38"):
        self.reference_genome = reference_genome
        self.gene_annotations = {}  # Will be populated with gene coordinates
        logger.info(f"Initialized VariantProcessor with {reference_genome}")

    def load_vcf(self, vcf_path: str) -> pd.DataFrame:
        """Load and parse VCF file into a DataFrame."""
        vcf = cyvcf2.VCF(vcf_path)
        variants = []
        
        for variant in vcf:
            variants.append({
                'chrom': variant.CHROM,
                'pos': variant.POS,
                'id': variant.ID,
                'ref': variant.REF,
                'alt': ','.join(variant.ALT),
                'qual': variant.QUAL,
                'filter': variant.FILTER,
                # Placeholder for gene annotation
                'gene': self.annotate_gene(variant.CHROM, variant.POS)
            })
            
        return pd.DataFrame(variants)

    def annotate_gene(self, chrom: str, pos: int) -> str:
        """Annotate a variant with gene information based on its position.
        
        Args:
            chrom: Chromosome of the variant
            pos: Position of the variant
            
        Returns:
            Gene name as string
        """
        for gene, coords in self.gene_annotations.items():
            if (chrom == coords['chrom'] and 
                coords['start'] <= pos <= coords['end']):
                return gene
        return "Unknown"

    def annotate_variants(self, variants_df: pd.DataFrame, 
                         encode_data: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """Annotate variants with functional information and ENCODE data.
        
        Args:
            variants_df: DataFrame of variants
            encode_data: Optional DataFrame with ENCODE annotations
            
        Returns:
            Annotated variant DataFrame
        """
        # Add basic annotations
        annotated_df = variants_df.copy()
        
        # Add ENCODE data if available
        if encode_data is not None:
            annotated_df = pd.merge(
                annotated_df, 
                encode_data,
                how='left',
                on=['chrom', 'pos']
            )
            
        return annotated_df

    def load_chromosome_vcf(self, chrom: str, vcf_dir: str = "data/1000genomes/") -> pd.DataFrame:
        """Load and parse VCF file for a specific chromosome.
        
        Args:
            chrom: Chromosome identifier (e.g., '19', 'X', 'Y')
            vcf_dir: Directory containing chromosome-specific VCF files
            
        Returns:
            DataFrame of variants for the specified chromosome
        """
        vcf_path = Path(vcf_dir) / f"chr{chrom}_variants.vcf"
        if not vcf_path.exists():
            vcf_path_gz = Path(vcf_dir) / f"chr{chrom}_variants.vcf.gz"
            if vcf_path_gz.exists():
                vcf_path = vcf_path_gz
            else:
                logger.error(f"No VCF file found for chromosome {chrom} in {vcf_dir}")
                return pd.DataFrame()
                
        logger.info(f"Loading variants from chromosome {chrom}")
        return self.load_vcf(str(vcf_path))

    def load_multiple_chromosomes(self, chroms: List[str], vcf_dir: str = "data/1000genomes/") -> pd.DataFrame:
        """Load and combine variants from multiple chromosomes.
        
        Args:
            chroms: List of chromosome identifiers (e.g., ['19', 'X', 'Y'])
            vcf_dir: Directory containing chromosome-specific VCF files
            
        Returns:
            Combined DataFrame of variants from all specified chromosomes
        """
        all_variants = []
        
        for chrom in chroms:
            chrom_variants = self.load_chromosome_vcf(chrom, vcf_dir)
            if not chrom_variants.empty:
                all_variants.append(chrom_variants)
                
        if not all_variants:
            logger.warning("No variant data loaded from any chromosome")
            return pd.DataFrame()
            
        logger.info(f"Combined variants from {len(all_variants)} chromosomes")
        return pd.concat(all_variants, ignore_index=True)

    def score_variant_impact(self, variant_info: Dict, use_api: bool = True) -> float:
        """Score the potential functional impact of a variant.
        
        Args:
            variant_info: Dictionary containing variant information
            use_api: Whether to use the NVIDIA Evo2 API for scoring
            
        Returns:
            Impact score between 0 and 1
        """
        # Use NVIDIA Evo2 API for scoring if requested
        if use_api:
            from evo2_api import Evo2API
            api = Evo2API()
            
            # Format variant for API
            api_variant = {
                "chrom": variant_info.get("chrom", ""),
                "pos": variant_info.get("pos", 0),
                "ref": variant_info.get("ref", ""),
                "alt": variant_info.get("alt", "").split(",")[0],  # Take first alt if multiple
                "id": variant_info.get("id", "")
            }
            
            # Call API and get score
            result_df = api.score_variants([api_variant])
            
            if not result_df.empty and "evo2_score" in result_df.columns:
                return float(result_df.iloc[0]["evo2_score"])
        
        # Fallback to local scoring if API is not used or fails
        impact_score = 0.0

        # Score based on variant type
        variant_type = variant_info.get("variant_type", "").upper()
        if variant_type == "SNP":
            impact_score += 0.2
        elif variant_type == "INDEL":
            impact_score += 0.3
        elif variant_type == "SV":
            impact_score += 0.4

        # Add scores based on genomic context
        if variant_info.get("in_exon"):
            impact_score += 0.3
        elif variant_info.get("in_promoter"):
            impact_score += 0.2
        elif variant_info.get("in_enhancer"):
            impact_score += 0.15

        # Normalize score to be between 0 and 1
        return min(impact_score, 1.0)

    def analyze_gene_variants(self, variants_df: pd.DataFrame,
                            gene_id: str,
                            region_info: Optional[Dict] = None,
                            use_api: bool = True) -> pd.DataFrame:
        """Analyze variants in a specific gene or genomic region.
        
        Args:
            variants_df: DataFrame of variants
            gene_id: ID of the gene to analyze
            region_info: Optional dictionary with region coordinates
            use_api: Whether to use the NVIDIA Evo2 API for analysis
            
        Returns:
            DataFrame with variant analysis
        """
        # Use NVIDIA Evo2 API for full gene analysis if requested
        if use_api and gene_id != "Unknown":
            try:
                from evo2_api import Evo2API
                api = Evo2API()
                
                # Get comprehensive gene analysis from API
                logger.info(f"Using NVIDIA Evo2 API for analysis of gene: {gene_id}")
                gene_result = api.analyze_gene(gene_id)
                
                if not isinstance(gene_result, dict) or "error" in gene_result:
                    logger.warning(f"API gene analysis failed, falling back to local analysis: {gene_result.get('error', 'Unknown error')}")
                else:
                    # Extract variants from API response
                    if "variants" in gene_result and isinstance(gene_result["variants"], list):
                        api_variants_df = pd.DataFrame(gene_result["variants"])
                        if not api_variants_df.empty:
                            logger.info(f"Retrieved {len(api_variants_df)} variants from API for gene {gene_id}")
                            # Add analysis columns
                            api_variants_df["analyzed_gene"] = gene_id
                            return api_variants_df
            except Exception as e:
                logger.error(f"Error in API gene analysis: {str(e)}")
        
        # Fall back to local analysis if API is not used or fails
        # Filter for gene variants
        if region_info:
            gene_variants = variants_df[
                (variants_df['chrom'] == region_info['chrom']) &
                (variants_df['pos'] >= region_info['start']) &
                (variants_df['pos'] <= region_info['end'])
            ].copy()
        else:
            gene_variants = variants_df[
                variants_df['gene'] == gene_id
            ].copy()

        # Add analysis columns
        gene_variants["analyzed_gene"] = gene_id
        gene_variants["impact_score"] = gene_variants.apply(
            lambda row: self.score_variant_impact(row.to_dict(), use_api=use_api),
            axis=1
        )

        return gene_variants

    def analyze_gwas_variants(self, variants_df: pd.DataFrame,
                            gwas_catalog: pd.DataFrame) -> pd.DataFrame:
        """Analyze variants based on GWAS catalog data.
        
        Args:
            variants_df: DataFrame of variants
            gwas_catalog: DataFrame with GWAS annotations
            
        Returns:
            DataFrame with GWAS variant analysis
        """
        # Merge variants with GWAS data
        gwas_variants = pd.merge(
            variants_df,
            gwas_catalog,
            how='left',
            left_on='variant_id',
            right_on='SNPS'
        )
        
        return gwas_variants


def main():
    """Main function to demonstrate variant analysis functionality."""
    processor = VariantProcessor()
    
    # Example usage
    variants_df = processor.load_vcf("data/1000genomes/chr19_variants.vcf")
    annotated_df = processor.annotate_variants(variants_df)
    
    # Generate report
    processor.generate_variant_report(annotated_df, "variant_report.txt")


if __name__ == "__main__":
    main()