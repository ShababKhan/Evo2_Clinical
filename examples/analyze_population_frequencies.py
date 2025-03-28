#!/usr/bin/env python3
"""
Example script demonstrating population frequency analysis workflow.
This script analyzes variant frequencies across different populations
and generates comparative visualizations.
"""

import sys
from pathlib import Path

# Add the parent directory of 'evo2_pipeline' to the Python module search path
sys.path.append(str(Path(__file__).resolve().parent.parent))

import pandas as pd
import numpy as np
from variant_analysis import VariantProcessor
from visualization import VariantVisualizer

# Define population groups
POPULATION_GROUPS = {
    'EUR': 'European',
    'EAS': 'East Asian',
    'SAS': 'South Asian',
    'AFR': 'African',
    'AMR': 'Admixed American'
}

def load_1000g_data(vcf_path: str) -> pd.DataFrame:
    """Load and preprocess 1000 Genomes data."""
    processor = VariantProcessor()
    variants_df = processor.load_vcf(vcf_path)
    
    # Add simulated population frequencies for demonstration
    for pop in POPULATION_GROUPS.keys():
        variants_df[pop] = np.random.uniform(0, 1, size=len(variants_df))
    
    return variants_df

def main():
    # Initialize processors
    processor = VariantProcessor()
    visualizer = VariantVisualizer()
    
    # Load and process data
    vcf_path = "data/1000g_variants.vcf"
    variants_df = load_1000g_data(vcf_path)
    
    # Create output directory
    output_dir = Path("results/population_analysis")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Analyze population frequencies
    for pop in POPULATION_GROUPS.keys():
        pop_variants = processor.analyze_population_frequencies(
            variants_df,
            population=pop
        )
        
        # Save population-specific data
        pop_variants.to_csv(
            output_dir / f"{pop}_frequencies.csv",
            index=False
        )
    
    # Generate visualizations
    visualizer.plot_population_frequencies(
        variants_df,
        populations=list(POPULATION_GROUPS.keys()),
        save_path=output_dir / "population_frequencies.png"
    )
    
    # Additional analysis for specific variant types
    variant_types = ['missense', 'synonymous', 'regulatory']
    for v_type in variant_types:
        type_variants = variants_df[variants_df['variant_type'] == v_type]
        if len(type_variants) > 0:
            visualizer.plot_population_frequencies(
                type_variants,
                populations=list(POPULATION_GROUPS.keys()),
                save_path=output_dir / f"{v_type}_frequencies.png"
            )
    
    # Generate comprehensive report
    visualizer.create_report_figures(
        variants_df,
        str(output_dir / "report")
    )
    
    print(f"Population analysis complete. Results saved in {output_dir}")

if __name__ == "__main__":
    main()