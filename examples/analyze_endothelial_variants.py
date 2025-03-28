#!/usr/bin/env python3
"""
Example script demonstrating gene set analysis workflow.
This script shows how to analyze variants across a set of related genes
and their functional impacts.
"""

import pandas as pd
from pathlib import Path
import argparse
import sys
import os

# Add parent directory to path so we can import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from variant_analysis import VariantProcessor
from visualization import VariantVisualizer

def main():
    parser = argparse.ArgumentParser(description='Analyze variants in a set of genes')
    parser.add_argument('--genes', type=str, nargs='+', required=True,
                      help='List of genes to analyze')
    parser.add_argument('--vcf', type=str, default='data/variants.vcf',
                      help='Path to input VCF file')
    parser.add_argument('--output-dir', type=str, default='results/gene_set_analysis',
                      help='Output directory for results')
    
    args = parser.parse_args()
    
    # Initialize processors
    processor = VariantProcessor()
    visualizer = VariantVisualizer()
    
    # Load VCF data
    variants_df = processor.load_vcf(args.vcf)
    
    # Filter variants for genes of interest
    gene_variants = variants_df[variants_df['gene'].isin(args.genes)].copy()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate visualizations
    visualizer.plot_variant_distribution(
        gene_variants,
        save_path=output_dir / "variant_distribution.png"
    )
    
    visualizer.plot_impact_heatmap(
        gene_variants,
        genes=args.genes,
        save_path=output_dir / "impact_heatmap.png"
    )
    
    # Add ENCODE data analysis
    encode_features = ['DNase', 'H3K27ac', 'H3K4me3']
    visualizer.plot_encode_integration(
        gene_variants,
        encode_features=encode_features,
        save_path=output_dir / "encode_features.png"
    )
    
    # Generate report
    visualizer.create_report_figures(
        gene_variants,
        str(output_dir / "report")
    )
    
    print(f"\nAnalysis Summary:")
    print(f"Total variants analyzed: {len(gene_variants)}")
    print(f"Genes analyzed: {', '.join(args.genes)}")
    for gene in args.genes:
        count = len(gene_variants[gene_variants['gene'] == gene])
        print(f"  {gene}: {count} variants")
    print(f"\nResults saved in {output_dir}")

if __name__ == "__main__":
    main()