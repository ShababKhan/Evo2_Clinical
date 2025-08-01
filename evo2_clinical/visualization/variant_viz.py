"""
Variant Visualization Module.

This module provides visualization components for genetic variants,
their effects, and functional impacts.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional, Union, Tuple
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import io
import base64


class VariantVisualizer:
    """Class for visualizing genetic variants and their effects."""
    
    def __init__(self, theme: str = "light"):
        """
        Initialize the variant visualizer.
        
        Args:
            theme (str): Color theme for visualizations ("light" or "dark")
        """
        self.theme = theme
        self._setup_styles()
    
    def _setup_styles(self):
        """Set up visualization styles based on theme."""
        if self.theme == "dark":
            self.colors = {
                "background": "#1E1E1E",
                "text": "#FFFFFF",
                "primary": "#BB86FC",
                "secondary": "#03DAC6",
                "accent": "#CF6679",
                "neutral": "#BBBBBB",
            }
            plt.style.use("dark_background")
        else:
            self.colors = {
                "background": "#FFFFFF",
                "text": "#000000",
                "primary": "#6200EE",
                "secondary": "#03DAC6",
                "accent": "#B00020",
                "neutral": "#555555",
            }
            plt.style.use("default")
    
    def plot_variant_distribution(
        self, 
        variants_df: pd.DataFrame, 
        groupby_column: str = "chrom",
        title: str = "Variant Distribution by Chromosome",
        interactive: bool = True
    ) -> Union[go.Figure, plt.Figure]:
        """
        Plot the distribution of variants.
        
        Args:
            variants_df (pd.DataFrame): DataFrame with variant information
            groupby_column (str): Column to group variants by (default: "chrom")
            title (str): Plot title
            interactive (bool): Whether to use interactive Plotly (True) 
                              or static Matplotlib (False)
        
        Returns:
            Union[go.Figure, plt.Figure]: Plotly or Matplotlib figure
        """
        if variants_df.empty:
            return None
            
        # Count variants by groupby column
        if groupby_column in variants_df.columns:
            counts = variants_df[groupby_column].value_counts().sort_index()
        else:
            raise ValueError(f"Column '{groupby_column}' not found in variants DataFrame")
        
        # Create interactive or static plot
        if interactive:
            fig = px.bar(
                x=counts.index,
                y=counts.values,
                labels={"x": groupby_column, "y": "Count"},
                title=title,
                color=counts.values,
                color_continuous_scale="viridis"
            )
            fig.update_layout(
                xaxis_title=groupby_column,
                yaxis_title="Number of Variants",
                template="plotly_white" if self.theme == "light" else "plotly_dark"
            )
            return fig
        else:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.bar(counts.index.astype(str), counts.values, color=self.colors["primary"])
            ax.set_xlabel(groupby_column)
            ax.set_ylabel("Number of Variants")
            ax.set_title(title)
            plt.xticks(rotation=45)
            plt.tight_layout()
            return fig

    def plot_variant_scores(
        self, 
        variants_df: pd.DataFrame, 
        score_column: str,
        threshold: Optional[float] = None,
        title: str = "Variant Functional Impact Scores",
        interactive: bool = True
    ) -> Union[go.Figure, plt.Figure]:
        """
        Plot variant scores with optional threshold.
        
        Args:
            variants_df (pd.DataFrame): DataFrame with variant information
            score_column (str): Column containing scores to plot
            threshold (float, optional): Highlight scores above this threshold
            title (str): Plot title
            interactive (bool): Whether to use interactive Plotly or static Matplotlib
        
        Returns:
            Union[go.Figure, plt.Figure]: Plotly or Matplotlib figure
        """
        if variants_df.empty or score_column not in variants_df.columns:
            return None
        
        # Sort by score for better visualization
        sorted_df = variants_df.sort_values(score_column)
        
        # Create index labels - try to use variant id if available
        if "id" in sorted_df.columns:
            labels = sorted_df["id"]
        elif "variant_id" in sorted_df.columns:
            labels = sorted_df["variant_id"]
        else:
            # Create labels from chromosome and position
            if "chrom" in sorted_df.columns and "pos" in sorted_df.columns:
                labels = sorted_df.apply(lambda row: f"{row['chrom']}:{row['pos']}", axis=1)
            else:
                labels = [f"Variant {i}" for i in range(len(sorted_df))]
        
        scores = sorted_df[score_column]
        
        # Create interactive or static plot
        if interactive:
            # Color based on threshold
            if threshold is not None:
                colors = [self.colors["accent"] if s >= threshold else self.colors["primary"] 
                         for s in scores]
                
                fig = px.bar(
                    x=labels,
                    y=scores,
                    color=scores,
                    labels={"x": "Variant", "y": score_column},
                    title=title,
                    color_continuous_scale="RdBu_r"
                )
                
                # Add threshold line
                fig.add_shape(
                    type="line",
                    x0=-0.5,
                    x1=len(scores)-0.5,
                    y0=threshold,
                    y1=threshold,
                    line=dict(color="red", width=2, dash="dash")
                )
            else:
                fig = px.bar(
                    x=labels,
                    y=scores,
                    labels={"x": "Variant", "y": score_column},
                    title=title,
                    color=scores,
                    color_continuous_scale="viridis"
                )
                
            fig.update_layout(
                xaxis_title="Variant",
                yaxis_title=score_column,
                xaxis_tickangle=45,
                template="plotly_white" if self.theme == "light" else "plotly_dark"
            )
            return fig
        else:
            fig, ax = plt.subplots(figsize=(12, 6))
            
            if threshold is not None:
                colors = [self.colors["accent"] if s >= threshold else self.colors["primary"] 
                         for s in scores]
                ax.bar(range(len(scores)), scores, color=colors)
                ax.axhline(y=threshold, color='red', linestyle='--', alpha=0.7)
            else:
                ax.bar(range(len(scores)), scores, color=self.colors["primary"])
            
            ax.set_xlabel("Variant")
            ax.set_ylabel(score_column)
            ax.set_title(title)
            
            # Add variant labels but limit for readability
            max_labels = min(20, len(labels))
            step = max(1, len(labels) // max_labels)
            plt.xticks(
                range(0, len(scores), step), 
                labels[::step], 
                rotation=45, 
                ha='right'
            )
            
            plt.tight_layout()
            return fig
    
    def plot_variant_lollipop(
        self, 
        variants_df: pd.DataFrame,
        gene_column: str = "gene",
        score_column: str = "impact_score",
        position_column: str = "pos",
        gene: Optional[str] = None,
        interactive: bool = True
    ) -> Union[go.Figure, plt.Figure]:
        """
        Create a lollipop plot of variants along a gene/chromosome.
        
        Args:
            variants_df (pd.DataFrame): DataFrame with variant information
            gene_column (str): Column with gene names
            score_column (str): Column with variant impact scores
            position_column (str): Column with genomic positions
            gene (str, optional): Filter to specific gene
            interactive (bool): Whether to use interactive Plotly or static Matplotlib
            
        Returns:
            Union[go.Figure, plt.Figure]: Plotly or Matplotlib figure
        """
        if variants_df.empty:
            return None
            
        # Filter by gene if specified
        if gene and gene_column in variants_df.columns:
            plot_df = variants_df[variants_df[gene_column] == gene].copy()
            if plot_df.empty:
                return None
        else:
            plot_df = variants_df.copy()
        
        # Need position and score columns
        if position_column not in plot_df.columns or score_column not in plot_df.columns:
            return None
            
        # Sort by position
        plot_df = plot_df.sort_values(position_column)
        
        # Create variant labels
        if "id" in plot_df.columns:
            labels = plot_df["id"]
        else:
            if "ref" in plot_df.columns and "alt" in plot_df.columns:
                labels = plot_df.apply(
                    lambda row: f"{row['ref']}>{row['alt']}", axis=1
                )
            else:
                labels = [f"Var{i}" for i in range(len(plot_df))]
        
        positions = plot_df[position_column]
        scores = plot_df[score_column]
        
        # Create title
        if gene:
            title = f"Variants in {gene}"
        else:
            title = "Variant Positions and Impact Scores"
        
        # Create interactive or static plot
        if interactive:
            fig = go.Figure()
            
            # Add stems
            fig.add_trace(
                go.Scatter(
                    x=positions,
                    y=scores,
                    mode='lines',
                    line=dict(color='gray', width=1),
                    hoverinfo='skip',
                    showlegend=False
                )
            )
            
            # Add markers
            fig.add_trace(
                go.Scatter(
                    x=positions,
                    y=scores,
                    mode='markers',
                    marker=dict(
                        size=10,
                        color=scores,
                        colorscale='Viridis',
                        showscale=True,
                        colorbar=dict(title=score_column)
                    ),
                    text=labels,
                    hovertemplate="<b>%{text}</b><br>" + 
                                 f"Position: %{{x}}<br>" + 
                                 f"{score_column}: %{{y:.3f}}"
                )
            )
            
            fig.update_layout(
                title=title,
                xaxis_title="Genomic Position",
                yaxis_title=score_column,
                template="plotly_white" if self.theme == "light" else "plotly_dark"
            )
            
            return fig
        else:
            fig, ax = plt.subplots(figsize=(12, 6))
            
            # Plot stems
            for x, y in zip(positions, scores):
                ax.plot([x, x], [0, y], color='gray', alpha=0.5)
            
            # Plot markers
            scatter = ax.scatter(
                positions, 
                scores, 
                c=scores, 
                cmap='viridis', 
                s=80, 
                zorder=2
            )
            
            # Add colorbar
            plt.colorbar(scatter, label=score_column)
            
            ax.set_title(title)
            ax.set_xlabel("Genomic Position")
            ax.set_ylabel(score_column)
            
            plt.tight_layout()
            return fig
    
    def plot_variant_effect_comparison(
        self, 
        variants_df: pd.DataFrame,
        effect_columns: List[str],
        variant_id_column: str = "id",
        title: str = "Comparison of Variant Effect Scores",
        interactive: bool = True
    ) -> Union[go.Figure, plt.Figure]:
        """
        Create a comparative visualization of different effect scores for variants.
        
        Args:
            variants_df (pd.DataFrame): DataFrame with variant information
            effect_columns (List[str]): List of columns with different effect scores
            variant_id_column (str): Column with variant identifiers
            title (str): Plot title
            interactive (bool): Whether to use interactive Plotly or static Matplotlib
            
        Returns:
            Union[go.Figure, plt.Figure]: Plotly or Matplotlib figure
        """
        if variants_df.empty:
            return None
            
        # Check if all required columns exist
        missing_cols = [col for col in effect_columns if col not in variants_df.columns]
        if missing_cols:
            print(f"Warning: Missing columns: {missing_cols}")
            effect_columns = [col for col in effect_columns if col in variants_df.columns]
            if not effect_columns:
                return None
                
        # Get variant IDs
        if variant_id_column in variants_df.columns:
            variant_ids = variants_df[variant_id_column]
        else:
            # Generate IDs from chromosome and position if available
            if "chrom" in variants_df.columns and "pos" in variants_df.columns:
                variant_ids = variants_df.apply(
                    lambda row: f"{row['chrom']}:{row['pos']}", axis=1
                )
            else:
                variant_ids = [f"Variant {i}" for i in range(len(variants_df))]
        
        # Prepare data for plotting
        if interactive:
            # Create a radar chart for each variant
            fig = go.Figure()
            
            for i, variant_id in enumerate(variant_ids):
                values = [variants_df.iloc[i][col] for col in effect_columns]
                # Close the polygon
                values.append(values[0])
                effect_cols_plot = effect_columns + [effect_columns[0]]
                
                fig.add_trace(go.Scatterpolar(
                    r=values,
                    theta=effect_cols_plot,
                    fill='toself',
                    name=str(variant_id)
                ))
            
            fig.update_layout(
                title=title,
                polar=dict(
                    radialaxis=dict(
                        visible=True,
                        range=[0, 1]
                    )
                ),
                showlegend=True,
                template="plotly_white" if self.theme == "light" else "plotly_dark"
            )
            
            return fig
        else:
            # Create a heatmap for static visualization
            plt.figure(figsize=(12, 8))
            
            # Create a matrix of scores
            score_matrix = variants_df[effect_columns].values
            
            # Create heatmap
            ax = sns.heatmap(
                score_matrix, 
                annot=True, 
                fmt=".2f",
                xticklabels=effect_columns,
                yticklabels=variant_ids,
                cmap="viridis"
            )
            
            plt.title(title)
            plt.tight_layout()
            
            return plt.gcf()
    
    def get_variant_summary_table(
        self, 
        variants_df: pd.DataFrame,
        include_columns: Optional[List[str]] = None,
        style: bool = True
    ) -> pd.DataFrame:
        """
        Generate a styled summary table of variants.
        
        Args:
            variants_df (pd.DataFrame): DataFrame with variant information
            include_columns (List[str], optional): Specific columns to include
            style (bool): Whether to apply styling (for notebook display)
            
        Returns:
            pd.DataFrame: Summary table with optional styling
        """
        if variants_df.empty:
            return pd.DataFrame()
        
        # Select columns to include
        if include_columns:
            # Only include columns that exist
            cols_to_include = [col for col in include_columns if col in variants_df.columns]
            if not cols_to_include:
                # If none of the specified columns exist, use default subset
                summary_df = self._get_default_variant_columns(variants_df)
            else:
                summary_df = variants_df[cols_to_include].copy()
        else:
            # Use default subset of columns
            summary_df = self._get_default_variant_columns(variants_df)
        
        if not style:
            return summary_df
        
        # Apply styling for better visualization
        # This will only have an effect when displayed in notebooks
        styled_df = summary_df.style
        
        # Add styling for score columns
        score_columns = [col for col in summary_df.columns if "score" in col.lower()]
        for col in score_columns:
            styled_df = styled_df.background_gradient(cmap="viridis", subset=[col])
        
        return styled_df
    
    def plot_variant_heatmap(
        self, 
        variants_df: pd.DataFrame,
        x_column: str,
        y_column: str,
        value_column: str,
        title: str = "Variant Effect Heatmap",
        interactive: bool = True
    ) -> Union[go.Figure, plt.Figure]:
        """
        Create a heatmap visualization of variant effects.
        
        Args:
            variants_df (pd.DataFrame): DataFrame with variant information
            x_column (str): Column to use for x-axis categories
            y_column (str): Column to use for y-axis categories
            value_column (str): Column containing values for color intensity
            title (str): Plot title
            interactive (bool): Whether to use interactive Plotly or static Matplotlib
            
        Returns:
            Union[go.Figure, plt.Figure]: Plotly or Matplotlib figure
        """
        if variants_df.empty:
            return None
            
        # Check if required columns exist
        for col in [x_column, y_column, value_column]:
            if col not in variants_df.columns:
                return None
        
        # Create pivot table for heatmap
        try:
            pivot_df = variants_df.pivot_table(
                values=value_column, 
                index=y_column, 
                columns=x_column, 
                aggfunc='mean'
            )
        except:
            # If pivot fails, try a different approach
            grouped = variants_df.groupby([y_column, x_column])[value_column].mean().reset_index()
            pivot_df = grouped.pivot(index=y_column, columns=x_column, values=value_column)
        
        # Fill NaN values with 0
        pivot_df = pivot_df.fillna(0)
        
        if interactive:
            fig = px.imshow(
                pivot_df,
                labels=dict(
                    x=x_column,
                    y=y_column,
                    color=value_column
                ),
                x=pivot_df.columns,
                y=pivot_df.index,
                color_continuous_scale="viridis",
                title=title
            )
            
            fig.update_layout(
                template="plotly_white" if self.theme == "light" else "plotly_dark"
            )
            
            return fig
        else:
            plt.figure(figsize=(12, 8))
            
            sns.heatmap(
                pivot_df, 
                cmap="viridis",
                annot=True, 
                fmt=".2f"
            )
            
            plt.title(title)
            plt.tight_layout()
            
            return plt.gcf()
    
    def _get_default_variant_columns(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Get a default subset of columns for variant summary."""
        # Define priority order for columns to include
        priority_columns = [
            "id", "variant_id", "chrom", "pos", "ref", "alt", "gene", "impact_score",
            "functional_score", "effect", "consequence"
        ]
        
        # Get columns that exist in the DataFrame
        available_cols = [col for col in priority_columns if col in variants_df.columns]
        
        # If we have less than 3 columns, include more columns from the DataFrame
        if len(available_cols) < 3:
            available_cols = list(variants_df.columns[:min(10, len(variants_df.columns))])
            
        return variants_df[available_cols].copy()
    
    def fig_to_base64(self, fig: Union[go.Figure, plt.Figure]) -> str:
        """
        Convert a figure to base64 string for embedding in HTML/Markdown.
        
        Args:
            fig: Matplotlib or Plotly figure object
            
        Returns:
            str: Base64 encoded image string
        """
        if fig is None:
            return ""
            
        buffer = io.BytesIO()
        
        # Handle different figure types
        if isinstance(fig, go.Figure):
            fig.write_image(buffer, format="png")
        else:
            # Matplotlib figure
            fig.savefig(buffer, format="png", bbox_inches="tight")
            plt.close(fig)
        
        # Get the image as base64 string
        buffer.seek(0)
        image_base64 = base64.b64encode(buffer.read()).decode("utf-8")
        return f"data:image/png;base64,{image_base64}"