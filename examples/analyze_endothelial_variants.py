#!/usr/bin/env python3
"""
Example script demonstrating endothelial variant analysis workflow.
This script analyzes variants in key endothelial genes and generates
visualizations of their distribution and impact.
"""

import pandas as pd
from pathlib import Path
from variant_analysis import VariantProcessor
from visualization import VariantVisualizer

# Define key endothelial genes
ENDOTHELIAL_GENES = [
    'GATA2',    # Key endothelial transcription factor
    'EPAS1',    # HIF-2α, important for angiogenesis
    'KDR',      # VEGFR2, critical for vascular development
    'TEK',      # TIE2, angiopoietin receptor
    'NOS3',     # eNOS, endothelial nitric oxide synthase
]

def main():
    # Initialize processors
    processor = VariantProcessor()
    visualizer = VariantVisualizer()
    
    # Load VCF data (using example VCF path)
    vcf_path = "data/variants.vcf"
    variants_df = processor.load_vcf(vcf_path)
    
    # Analyze endothelial variants
    endothelial_variants = processor.analyze_endothelial_variants(
        variants_df,
        ENDOTHELIAL_GENES
    )
    
    # Create output directory
    output_dir = Path("results/endothelial_analysis")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate visualizations
    visualizer.plot_variant_distribution(
        endothelial_variants,
        save_path=output_dir / "endothelial_variant_distribution.png"
    )
    
    visualizer.plot_impact_heatmap(
        endothelial_variants,
        genes=ENDOTHELIAL_GENES,
        save_path=output_dir / "endothelial_impact_heatmap.png"
    )
    
    # Add ENCODE data analysis
    encode_features = ['DNase', 'H3K27ac', 'H3K4me3']
    visualizer.plot_encode_integration(
        endothelial_variants,
        encode_features=encode_features,
        save_path=output_dir / "endothelial_encode_features.png"
    )
    
    # Generate report
    visualizer.create_report_figures(
        endothelial_variants,
        str(output_dir / "report")
    )
    
    print(f"Analysis complete. Results saved in {output_dir}")

if __name__ == "__main__":
    main()