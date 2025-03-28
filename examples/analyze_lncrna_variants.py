#!/usr/bin/env python3
"""
Example script demonstrating lncRNA variant analysis workflow.
This script focuses on analyzing variants in GATA2-AS1 and their
potential impact on endothelial function.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import os

# Add parent directory to path so we can import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from variant_analysis import VariantProcessor
from visualization import VariantVisualizer

# Define relevant genomic regions
GATA2_AS1_REGION = {
    'chrom': 'chr3',
    'start': 128198266,  # GRCh38 coordinates
    'end': 128207434,
    'name': 'GATA2-AS1'
}

# Define nearby regulatory elements
REGULATORY_ELEMENTS = [
    {'name': 'GATA2_enhancer', 'start': 128198000, 'end': 128199000},
    {'name': 'GATA2_promoter', 'start': 128207000, 'end': 128208000}
]

def analyze_regulatory_impact(variants_df: pd.DataFrame) -> pd.DataFrame:
    """Analyze variant impact on regulatory elements."""
    # Create a copy of the DataFrame to avoid modifying the original
    variants_df = variants_df.copy()
    
    # Initialize the regulatory_element column with None/NaN
    variants_df['regulatory_element'] = None
    
    for element in REGULATORY_ELEMENTS:
        mask = (
            (variants_df['pos'] >= element['start']) &
            (variants_df['pos'] <= element['end'])
        )
        # Using mask indexing instead of .loc with a mask
        variants_df.loc[mask, 'regulatory_element'] = element['name']

    return variants_df

def main():
    # Initialize processors
    processor = VariantProcessor()
    visualizer = VariantVisualizer()
    
    # Load VCF data
    vcf_path = "data/variants.vcf"
    variants_df = processor.load_vcf(vcf_path)
    
    # Filter for GATA2-AS1 region
    gata2as1_variants = variants_df[
        (variants_df['chrom'] == GATA2_AS1_REGION['chrom']) &
        (variants_df['pos'] >= GATA2_AS1_REGION['start']) &
        (variants_df['pos'] <= GATA2_AS1_REGION['end'])
    ].copy()
    
    # Analyze lncRNA variants
    gata2as1_variants = processor.analyze_lncrna_variants(
        gata2as1_variants,
        lncrna_id='GATA2-AS1'
    )
    
    # Add regulatory impact analysis
    gata2as1_variants = analyze_regulatory_impact(gata2as1_variants)
    
    # Add impact scores for visualization (required by plot_lncrna_analysis)
    gata2as1_variants['impact_score'] = np.random.uniform(0.1, 0.9, size=len(gata2as1_variants))
    
    # Create output directory
    output_dir = Path("results/lncrna_analysis")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate visualizations
    visualizer.plot_lncrna_analysis(
        gata2as1_variants,
        lncrna_id='GATA2-AS1',
        save_path=output_dir / "gata2as1_analysis.html"
    )
    
    # Plot variant distribution
    visualizer.plot_variant_distribution(
        gata2as1_variants,
        save_path=output_dir / "gata2as1_distribution.png"
    )
    
    # Analyze ENCODE data integration
    encode_features = ['DNase', 'H3K27ac', 'H3K4me3']
    
    # Add mock ENCODE data for visualization
    for feature in encode_features:
        gata2as1_variants[feature] = np.random.choice([0, 1], size=len(gata2as1_variants), p=[0.7, 0.3])
    
    visualizer.plot_encode_integration(
        gata2as1_variants,
        encode_features=encode_features,
        save_path=output_dir / "gata2as1_encode_features.png"
    )
    
    # Save variant data
    gata2as1_variants.to_csv(
        output_dir / "gata2as1_variants.csv",
        index=False
    )
    
    # Generate comprehensive report
    visualizer.create_report_figures(
        gata2as1_variants,
        str(output_dir / "report")
    )
    
    # Print summary statistics
    print("\nGATA2-AS1 Analysis Summary:")
    print(f"Total variants analyzed: {len(gata2as1_variants)}")
    print(f"Variants in regulatory elements: "
          f"{gata2as1_variants['regulatory_element'].notna().sum()}")
    print(f"\nResults saved in {output_dir}")

if __name__ == "__main__":
    main()