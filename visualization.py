#!/usr/bin/env python3
"""
Visualization Module for Evo2 Pipeline

This module provides visualization utilities for the endothelial research program,
specifically focusing on:
- Variant distribution plots
- Gene-level impact visualizations
- Population frequency comparisons
- ENCODE data integration plots
- Endothelial-specific feature visualization

Author: Shabab Khan
Date: March 26, 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from typing import Dict, List, Tuple, Optional, Union
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("VisualizationModule")

class VariantVisualizer:
    """Class for creating visualizations of variant analysis results."""
    
    def __init__(self):
        """Initialize the visualizer with default plotting settings."""
        sns.set_theme()  # Use seaborn's default theme
        self.default_colors = sns.color_palette("husl", 8)
        logger.info("Initialized VariantVisualizer with seaborn theme")
        
    def plot_variant_distribution(self, variants_df: pd.DataFrame,
                                chromosome: Optional[str] = None,
                                save_path: Optional[str] = None) -> None:
        """Plot the distribution of variants across the genome.
        
        Args:
            variants_df: DataFrame of variants
            chromosome: Optional chromosome to focus on
            save_path: Optional path to save the plot
        """
        plt.figure(figsize=(12, 6))
        
        if chromosome:
            data = variants_df[variants_df['chrom'] == chromosome]
            title = f"Variant Distribution on Chromosome {chromosome}"
        else:
            data = variants_df
            title = "Genome-wide Variant Distribution"
        
        if len(data) > 0:
            # Only attempt to plot if we have data
            sns.histplot(data=data, x='pos', hue='chrom', multiple="stack", kde=len(data) > 1)
        else:
            # If no data, just show an empty plot with a message
            plt.text(0.5, 0.5, "No variants found in this region",
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=plt.gca().transAxes)
        
        plt.title(title)
        plt.xlabel("Position")
        plt.ylabel("Count")
        
        if save_path:
            plt.savefig(save_path)
            logger.info(f"Saved variant distribution plot to {save_path}")
        plt.close()

    def plot_impact_heatmap(self, variants_df: pd.DataFrame,
                           genes: List[str],
                           save_path: Optional[str] = None) -> None:
        """Create a heatmap of variant impact scores across genes.
        
        Args:
            variants_df: DataFrame of variants
            genes: List of genes to include
            save_path: Optional path to save the plot
        """
        impact_matrix = variants_df.pivot_table(
            index='gene',
            columns='variant_type',
            values='impact_score',
            aggfunc='mean'
        ).fillna(0)
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(impact_matrix, annot=True, cmap='coolwarm', center=0.5, fmt=".2f")
        plt.title("Variant Impact Scores by Gene and Type")
        
        if save_path:
            plt.savefig(save_path)
            logger.info(f"Saved impact heatmap to {save_path}")
        plt.close()

    def plot_population_frequencies(self, variants_df: pd.DataFrame,
                                  populations: List[str],
                                  save_path: Optional[str] = None) -> None:
        """Plot variant frequencies across different populations.
        
        Args:
            variants_df: DataFrame of variants
            populations: List of population IDs
            save_path: Optional path to save the plot
        """
        plt.figure(figsize=(12, 6))
        
        freq_data = variants_df[populations].melt(var_name="Population", value_name="Frequency")
        sns.boxplot(data=freq_data, x='Population', y='Frequency', palette="Set2")
        
        plt.title("Variant Frequencies Across Populations")
        plt.xlabel("Population")
        plt.ylabel("Frequency")
        plt.xticks(rotation=45)
        
        if save_path:
            plt.savefig(save_path)
            logger.info(f"Saved population frequency plot to {save_path}")
        plt.close()

    def plot_encode_integration(self, variants_df: pd.DataFrame,
                              encode_features: List[str],
                              save_path: Optional[str] = None) -> None:
        """Visualize variant overlap with ENCODE features.
        
        Args:
            variants_df: DataFrame of variants
            encode_features: List of ENCODE feature types
            save_path: Optional path to save the plot
        """
        plt.figure(figsize=(10, 6))
        
        feature_counts = variants_df[encode_features].sum()
        sns.barplot(x=feature_counts.index, y=feature_counts.values, palette="viridis")
        
        plt.title("Variants Overlapping ENCODE Features")
        plt.xlabel("Feature Type")
        plt.ylabel("Number of Variants")
        plt.xticks(rotation=45)
        
        if save_path:
            plt.savefig(save_path)
            logger.info(f"Saved ENCODE integration plot to {save_path}")
        plt.close()

    def plot_lncrna_analysis(self, variants_df: pd.DataFrame,
                            lncrna_id: str = "GATA2-AS1",
                            save_path: Optional[str] = None) -> None:
        """Create visualization for lncRNA variant analysis.
        
        Args:
            variants_df: DataFrame of variants
            lncrna_id: ID of the lncRNA (default: GATA2-AS1)
            save_path: Optional path to save the plot
        """
        fig = go.Figure()
        
        # Add variant positions
        fig.add_trace(go.Scatter(
            x=variants_df['pos'],
            y=variants_df['impact_score'],
            mode='markers',
            name='Variants',
            marker=dict(
                size=10,
                color=variants_df['impact_score'],
                colorscale='Plasma',
                showscale=True
            )
        ))
        
        fig.update_layout(
            title=f"Variant Analysis for {lncrna_id}",
            xaxis_title="Position",
            yaxis_title="Impact Score",
            template="plotly_white"
        )
        
        if save_path:
            fig.write_html(save_path)
            logger.info(f"Saved lncRNA analysis plot to {save_path}")
            
    def create_report_figures(self, variants_df: pd.DataFrame,
                            output_dir: str) -> None:
        """Generate all figures for a comprehensive report.
        
        Args:
            variants_df: DataFrame of variants
            output_dir: Directory to save figures
        """
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True, parents=True)
        
        # Generate all plots
        self.plot_variant_distribution(
            variants_df,
            save_path=output_path / "variant_distribution.png"
        )
        
        # Only create impact heatmap if variant_type column exists
        if 'variant_type' in variants_df.columns and len(variants_df) > 0:
            self.plot_impact_heatmap(
                variants_df,
                genes=variants_df['gene'].unique(),
                save_path=output_path / "impact_heatmap.png"
            )
        
        # Only create ENCODE integration plot if ENCODE features exist
        encode_features = ['DNase', 'H3K27ac', 'H3K4me3']
        if all(feature in variants_df.columns for feature in encode_features) and len(variants_df) > 0:
            self.plot_encode_integration(
                variants_df,
                encode_features=encode_features,
                save_path=output_path / "encode_integration.png"
            )
        
        # Create lncRNA analysis plot if impact_score exists
        if 'impact_score' in variants_df.columns and len(variants_df) > 0:
            self.plot_lncrna_analysis(
                variants_df,
                save_path=output_path / "gata2as1_analysis.html"
            )
        
        logger.info(f"Generated report figures in {output_dir}")


def main():
    """Main function to demonstrate visualization functionality."""
    visualizer = VariantVisualizer()
    
    # Example usage with dummy data
    df = pd.DataFrame({
        'chrom': ['chr1'] * 100,
        'pos': np.random.randint(1, 1000000, 100),
        'impact_score': np.random.random(100),
        'gene': np.random.choice(['GATA2', 'EPAS1', 'KDR'], 100)
    })
    
    visualizer.create_report_figures(df, "example_output")


if __name__ == "__main__":
    main()