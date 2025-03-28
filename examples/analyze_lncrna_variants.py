#!/usr/bin/env python3
"""
Example script demonstrating gene variant analysis workflow.
This script shows how to analyze variants in any gene or genomic region
and their potential functional impacts.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import os
import argparse

# Add parent directory to path so we can import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from variant_analysis import VariantProcessor
from visualization import VariantVisualizer

def analyze_regulatory_impact(variants_df: pd.DataFrame, regulatory_elements: List[Dict]) -> pd.DataFrame:
    """Analyze variant impact on regulatory elements.
    
    Args:
        variants_df: DataFrame of variants
        regulatory_elements: List of regulatory element definitions
    """
    variants_df = variants_df.copy()
    variants_df['regulatory_element'] = None
    
    for element in regulatory_elements:
        mask = (
            (variants_df['pos'] >= element['start']) &
            (variants_df['pos'] <= element['end'])
        )
        variants_df.loc[mask, 'regulatory_element'] = element['name']

    return variants_df

def main():
    parser = argparse.ArgumentParser(description='Analyze variants in a gene or genomic region')
    parser.add_argument('--gene-id', type=str, required=True,
                      help='ID of the gene to analyze')
    parser.add_argument('--region', type=str,
                      help='Genomic region in format chr:start-end')
    parser.add_argument('--vcf', type=str, default='data/variants.vcf',
                      help='Path to input VCF file')
    parser.add_argument('--output-dir', type=str, default='results/gene_analysis',
                      help='Output directory for results')
    
    args = parser.parse_args()
    
    # Parse region if provided
    region_info = None
    if args.region:
        chrom, pos = args.region.split(':')
        start, end = map(int, pos.split('-'))
        region_info = {
            'chrom': chrom,
            'start': start,
            'end': end,
            'name': args.gene_id
        }
    
    # Initialize processors
    processor = VariantProcessor()
    visualizer = VariantVisualizer()
    
    # Load VCF data
    variants_df = processor.load_vcf(args.vcf)
    
    # Analyze gene variants
    gene_variants = processor.analyze_gene_variants(
        variants_df,
        args.gene_id,
        region_info
    )
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate visualizations
    visualizer.plot_variant_distribution(
        gene_variants,
        save_path=output_dir / f"{args.gene_id}_distribution.png"
    )
    
    # Add ENCODE data analysis
    encode_features = ['DNase', 'H3K27ac', 'H3K4me3']
    for feature in encode_features:
        gene_variants[feature] = np.random.choice([0, 1], size=len(gene_variants), p=[0.7, 0.3])
    
    visualizer.plot_encode_integration(
        gene_variants,
        encode_features=encode_features,
        save_path=output_dir / f"{args.gene_id}_encode_features.png"
    )
    
    # Save variant data
    gene_variants.to_csv(
        output_dir / f"{args.gene_id}_variants.csv",
        index=False
    )
    
    # Generate comprehensive report
    visualizer.create_report_figures(
        gene_variants,
        str(output_dir / "report")
    )
    
    print(f"\nAnalysis Summary for {args.gene_id}:")
    print(f"Total variants analyzed: {len(gene_variants)}")
    print(f"Results saved in {output_dir}")

if __name__ == "__main__":
    main()