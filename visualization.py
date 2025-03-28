#!/usr/bin/env python3
"""
Visualization Module for Evo2 Pipeline

This module provides visualization utilities for variant analysis results,
including distribution plots, impact heatmaps, and ENCODE feature integration.

Author: Shabab Khan
Date: March 26, 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Optional, Dict
import logging
from pathlib import Path
import plotly.graph_objects as go
import plotly.express as px

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("Visualization")

class VariantVisualizer:
    """Class for creating visualizations of variant analysis results."""
    
    def __init__(self):
        """Initialize the visualizer with default style settings."""
        plt.style.use('seaborn')
        self.default_figsize = (10, 6)
        self.color_palette = "viridis"
        
    def plot_variant_distribution(self, variants_df: pd.DataFrame,
                                save_path: Optional[str] = None) -> None:
        """Plot the distribution of variants across genomic positions.
        
        Args:
            variants_df: DataFrame of variants
            save_path: Optional path to save the plot
        """
        plt.figure(figsize=self.default_figsize)
        
        sns.histplot(
            data=variants_df,
            x='pos',
            hue='gene',
            multiple="stack",
            palette=self.color_palette
        )
        
        plt.title("Variant Distribution Across Genomic Positions")
        plt.xlabel("Genomic Position")
        plt.ylabel("Number of Variants")
        plt.xticks(rotation=45)
        
        if save_path:
            plt.savefig(save_path, bbox_inches='tight')
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
        plt.figure(figsize=self.default_figsize)
        
        # Prepare data for heatmap
        impact_matrix = variants_df.pivot_table(
            index='gene',
            columns='variant_type',
            values='impact_score',
            aggfunc='mean'
        ).fillna(0)
        
        # Plot heatmap
        sns.heatmap(
            impact_matrix,
            annot=True,
            cmap=self.color_palette,
            fmt='.2f'
        )
        
        plt.title("Variant Impact Scores by Gene and Type")
        plt.ylabel("Gene")
        plt.xlabel("Variant Type")
        
        if save_path:
            plt.savefig(save_path, bbox_inches='tight')
            logger.info(f"Saved impact heatmap to {save_path}")
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
        plt.figure(figsize=self.default_figsize)
        
        feature_counts = variants_df[encode_features].sum()
        sns.barplot(x=feature_counts.index, y=feature_counts.values, palette=self.color_palette)
        
        plt.title("Variants Overlapping ENCODE Features")
        plt.xlabel("Feature Type")
        plt.ylabel("Number of Variants")
        plt.xticks(rotation=45)
        
        if save_path:
            plt.savefig(save_path, bbox_inches='tight')
            logger.info(f"Saved ENCODE integration plot to {save_path}")
        plt.close()
        
    def create_interactive_plot(self, variants_df: pd.DataFrame,
                              x: str, y: str,
                              color: Optional[str] = None,
                              title: str = "Interactive Variant Plot",
                              save_path: Optional[str] = None) -> None:
        """Create an interactive scatter plot of variant data.
        
        Args:
            variants_df: DataFrame of variants
            x: Column name for x-axis
            y: Column name for y-axis
            color: Optional column name for color coding
            title: Plot title
            save_path: Optional path to save the plot
        """
        fig = px.scatter(
            variants_df,
            x=x,
            y=y,
            color=color,
            title=title,
            hover_data=['gene', 'variant_type', 'impact_score']
        )
        
        if save_path:
            fig.write_html(save_path)
            logger.info(f"Saved interactive plot to {save_path}")
            
    def create_report_figures(self, variants_df: pd.DataFrame,
                            output_dir: str) -> None:
        """Generate a complete set of figures for analysis report.
        
        Args:
            variants_df: DataFrame of variants
            output_dir: Directory to save report figures
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Basic distribution plot
        self.plot_variant_distribution(
            variants_df,
            save_path=str(output_path / "variant_distribution.png")
        )
        
        # Impact score distribution
        if 'impact_score' in variants_df.columns:
            plt.figure(figsize=self.default_figsize)
            sns.histplot(data=variants_df, x='impact_score', bins=20)
            plt.title("Distribution of Variant Impact Scores")
            plt.savefig(output_path / "impact_score_distribution.png")
            plt.close()
            
        # Gene-wise variant counts
        if 'gene' in variants_df.columns:
            plt.figure(figsize=self.default_figsize)
            gene_counts = variants_df['gene'].value_counts()
            sns.barplot(x=gene_counts.index, y=gene_counts.values)
            plt.title("Variants per Gene")
            plt.xticks(rotation=45)
            plt.savefig(output_path / "variants_per_gene.png")
            plt.close()
            
        logger.info(f"Generated report figures in {output_dir}")


def main():
    """Main function to demonstrate visualization functionality."""
    # Example usage
    visualizer = VariantVisualizer()
    
    # Create sample data
    data = {
        'pos': np.random.randint(1000, 2000, 100),
        'gene': np.random.choice(['GENE1', 'GENE2', 'GENE3'], 100),
        'impact_score': np.random.uniform(0, 1, 100),
        'variant_type': np.random.choice(['SNP', 'INDEL', 'SV'], 100)
    }
    df = pd.DataFrame(data)
    
    # Generate example plots
    visualizer.plot_variant_distribution(df, "example_distribution.png")
    visualizer.plot_impact_heatmap(df, ['GENE1', 'GENE2', 'GENE3'], "example_heatmap.png")


if __name__ == "__main__":
    main()