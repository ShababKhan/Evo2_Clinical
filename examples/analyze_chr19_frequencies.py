#!/usr/bin/env python3
"""
Example script demonstrating population frequency analysis workflow using chromosome 19 data.
This script analyzes variant frequencies across different populations from chromosome 19
and generates comparative visualizations.
"""
import sys
from pathlib import Path
# Add the parent directory to the Python module search path
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

def load_chr19_data(vcf_path: str) -> pd.DataFrame:
    """Load and preprocess chromosome 19 data."""
    processor = VariantProcessor()
    print(f"Loading VCF data from {vcf_path}...")
    variants_df = processor.load_vcf(vcf_path)
    print(f"Loaded {len(variants_df)} variants")
    
    return variants_df

def main():
    # Initialize processors
    processor = VariantProcessor()
    visualizer = VariantVisualizer()
    
    # Load and process data - using the chr19.vcf.gz file
    vcf_path = "data/chr19.vcf.gz"
    variants_df = load_chr19_data(vcf_path)
    
    # For demonstration purposes, limit to a smaller subset if the dataset is large
    if len(variants_df) > 1000:
        print(f"Limiting analysis to first 1000 variants for demonstration")
        variants_df = variants_df.head(1000)
    
    # Add simulated population frequencies for demonstration
    for pop in POPULATION_GROUPS.keys():
        variants_df[pop] = np.random.uniform(0, 1, size=len(variants_df))
    
    # Create output directory
    output_dir = Path("results/chr19_analysis")
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
        print(f"Saved {pop} frequency data")
    
    # Generate visualizations
    visualizer.plot_population_frequencies(
        variants_df,
        populations=list(POPULATION_GROUPS.keys()),
        save_path=output_dir / "population_frequencies.png"
    )
    print(f"Generated population frequency plot")
    
    # Generate comprehensive report
    visualizer.create_report_figures(
        variants_df,
        str(output_dir / "report")
    )
    print(f"Created report figures")
    
    print(f"Chromosome 19 population analysis complete. Results saved in {output_dir}")

if __name__ == "__main__":
    main()