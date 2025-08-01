"""
Pathway Visualization Module.

This module provides visualization components for biological pathways,
focusing on EMT pathways, gene networks, and signaling cascades.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional, Union, Tuple, Any
import plotly.express as px
import plotly.graph_objects as go
import networkx as nx
import io
import base64
from matplotlib.colors import LinearSegmentedColormap


class PathwayVisualizer:
    """Class for visualizing biological pathways and gene networks."""
    
    def __init__(self, theme: str = "light"):
        """
        Initialize the pathway visualizer.
        
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
                "nodes": {
                    "gene": "#BB86FC",
                    "variant": "#03DAC6",
                    "pathway": "#CF6679",
                    "transcription_factor": "#FFB74D",
                    "drug": "#64FFDA",
                    "phenotype": "#FF7043"
                },
                "edges": {
                    "activation": "#4FC3F7",
                    "inhibition": "#EF5350",
                    "unknown": "#BBBBBB"
                }
            }
            plt.style.use("dark_background")
        else:
            self.colors = {
                "background": "#FFFFFF",
                "text": "#000000",
                "nodes": {
                    "gene": "#7B1FA2",
                    "variant": "#00796B",
                    "pathway": "#C62828",
                    "transcription_factor": "#EF6C00",
                    "drug": "#00695C",
                    "phenotype": "#D84315"
                },
                "edges": {
                    "activation": "#0288D1",
                    "inhibition": "#D32F2F",
                    "unknown": "#757575"
                }
            }
            plt.style.use("default")
    
    def plot_pathway_network(
        self,
        nodes: pd.DataFrame,
        edges: pd.DataFrame,
        node_id_col: str = "id",
        source_col: str = "source",
        target_col: str = "target",
        node_type_col: str = "type",
        node_score_col: Optional[str] = "score",
        edge_type_col: Optional[str] = "interaction",
        title: str = "Pathway Network",
        interactive: bool = True,
        highlight_nodes: Optional[List[str]] = None
    ) -> Union[go.Figure, plt.Figure]:
        """
        Plot a network visualization of a biological pathway.
        
        Args:
            nodes (pd.DataFrame): DataFrame with node information
            edges (pd.DataFrame): DataFrame with edge information
            node_id_col (str): Column in nodes with node IDs
            source_col (str): Column in edges with source node IDs
            target_col (str): Column in edges with target node IDs
            node_type_col (str): Column in nodes with node types
            node_score_col (str, optional): Column in nodes with scores for sizing
            edge_type_col (str, optional): Column in edges with edge/interaction types
            title (str): Plot title
            interactive (bool): Whether to use interactive Plotly or static Matplotlib
            highlight_nodes (List[str], optional): List of node IDs to highlight
            
        Returns:
            Union[go.Figure, plt.Figure]: Plotly or Matplotlib figure
        """
        if nodes.empty or edges.empty:
            return None
        
        # Create networkx graph
        G = nx.DiGraph()
        
        # Add nodes with attributes
        for _, node in nodes.iterrows():
            node_id = str(node[node_id_col])
            
            # Set node attributes
            attrs = {"id": node_id}
            
            # Set node type
            if node_type_col in nodes.columns:
                node_type = node[node_type_col]
                attrs["type"] = node_type
            else:
                node_type = "gene"  # default
            
            # Set node score/size if available
            if node_score_col and node_score_col in nodes.columns:
                attrs["score"] = node[node_score_col]
            
            # Add other attributes
            for col in nodes.columns:
                if col not in [node_id_col]:
                    attrs[col] = node[col]
            
            G.add_node(node_id, **attrs)
        
        # Add edges with attributes
        for _, edge in edges.iterrows():
            source = str(edge[source_col])
            target = str(edge[target_col])
            
            # Skip if source or target doesn't exist
            if source not in G or target not in G:
                continue
            
            # Set edge attributes
            attrs = {}
            
            # Set edge type if available
            if edge_type_col and edge_type_col in edges.columns:
                attrs["interaction"] = edge[edge_type_col]
            
            # Add other attributes
            for col in edges.columns:
                if col not in [source_col, target_col]:
                    attrs[col] = edge[col]
            
            G.add_edge(source, target, **attrs)
        
        # Use layout algorithms for node positioning
        try:
            pos = nx.kamada_kawai_layout(G)
        except:
            pos = nx.spring_layout(G, seed=42)
        
        # Set node color and size
        node_colors = []
        node_sizes = []
        for node in G.nodes():
            node_type = G.nodes[node].get("type", "gene")
            node_color = self.colors["nodes"].get(node_type, self.colors["nodes"]["gene"])
            node_colors.append(node_color)
            
            # Set size based on score or default
            if "score" in G.nodes[node]:
                node_sizes.append(20 + 80 * G.nodes[node]["score"])
            else:
                node_sizes.append(30)
        
        # Set edge colors and widths
        edge_colors = []
        edge_widths = []
        for u, v, data in G.edges(data=True):
            interaction = data.get("interaction", "unknown")
            
            if interaction.lower() in ["activation", "activates", "activating"]:
                edge_color = self.colors["edges"]["activation"]
            elif interaction.lower() in ["inhibition", "inhibits", "inhibiting"]:
                edge_color = self.colors["edges"]["inhibition"]
            else:
                edge_color = self.colors["edges"]["unknown"]
            
            edge_colors.append(edge_color)
            
            # Set edge width based on weight or default
            if "weight" in data:
                edge_widths.append(1 + 3 * data["weight"])
            else:
                edge_widths.append(1.5)
        
        # Create visualization
        if interactive:
            # Create Plotly network visualization
            edge_traces = []
            
            # Add edges
            for (u, v, data), color, width in zip(G.edges(data=True), edge_colors, edge_widths):
                x0, y0 = pos[u]
                x1, y1 = pos[v]
                
                # Edge trace
                edge_trace = go.Scatter(
                    x=[x0, x1, None],
                    y=[y0, y1, None],
                    line=dict(width=width, color=color),
                    hoverinfo='text',
                    text=data.get("interaction", ""),
                    mode='lines',
                    showlegend=False
                )
                edge_traces.append(edge_trace)
                
                # Add arrow for directed graph
                # Calculate the position for the arrow (80% along the edge)
                arrow_x = x0 + 0.8*(x1-x0)
                arrow_y = y0 + 0.8*(y1-y0)
                
                # Calculate angle for arrow
                angle = np.arctan2(y1-y0, x1-x0)
                
                # Add arrow marker
                arrow_trace = go.Scatter(
                    x=[arrow_x],
                    y=[arrow_y],
                    mode='markers',
                    marker=dict(
                        symbol='triangle-right',
                        size=8,
                        color=color,
                        angle=np.degrees(angle)
                    ),
                    hoverinfo='skip',
                    showlegend=False
                )
                edge_traces.append(arrow_trace)
            
            # Node trace
            node_trace = go.Scatter(
                x=[pos[node][0] for node in G.nodes()],
                y=[pos[node][1] for node in G.nodes()],
                mode='markers',
                hoverinfo='text',
                marker=dict(
                    size=node_sizes,
                    color=node_colors,
                    line=dict(width=1, color='black')
                ),
                text=[f"ID: {node}<br>Type: {G.nodes[node].get('type', 'gene')}<br>"
                      f"Score: {G.nodes[node].get('score', 'N/A')}" 
                      for node in G.nodes()],
                name='Nodes'
            )
            
            # Create figure
            fig = go.Figure(data=edge_traces + [node_trace])
            
            # Add node labels
            node_labels = go.Scatter(
                x=[pos[node][0] for node in G.nodes()],
                y=[pos[node][1] for node in G.nodes()],
                mode='text',
                text=[node for node in G.nodes()],
                textposition='bottom center',
                hoverinfo='none',
                showlegend=False
            )
            fig.add_trace(node_labels)
            
            # Update layout
            fig.update_layout(
                title=title,
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20, l=5, r=5, t=40),
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                template="plotly_white" if self.theme == "light" else "plotly_dark"
            )
            
            return fig
        else:
            # Create Matplotlib network visualization
            plt.figure(figsize=(12, 10))
            
            # Draw edges
            nx.draw_networkx_edges(
                G, pos,
                width=edge_widths,
                edge_color=edge_colors,
                arrows=True,
                arrowsize=15,
                alpha=0.7
            )
            
            # Draw nodes
            nx.draw_networkx_nodes(
                G, pos,
                node_size=node_sizes,
                node_color=node_colors,
                edgecolors='black',
                linewidths=1,
                alpha=0.9
            )
            
            # Draw labels
            nx.draw_networkx_labels(
                G, pos,
                font_size=10,
                font_weight='bold',
                font_color='black' if self.theme == 'light' else 'white'
            )
            
            plt.title(title)
            plt.axis('off')
            plt.tight_layout()
            
            return plt.gcf()
    
    def plot_emt_pathway(
        self, 
        pathway_scores: Dict[str, float],
        gene_scores: Optional[Dict[str, float]] = None,
        title: str = "EMT Pathway Analysis",
        interactive: bool = True
    ) -> Union[go.Figure, plt.Figure]:
        """
        Plot EMT pathway scores with gene-level data if available.
        
        Args:
            pathway_scores (Dict[str, float]): Dictionary of pathway names and scores
            gene_scores (Dict[str, float], optional): Dictionary of gene names and scores
            title (str): Plot title
            interactive (bool): Whether to use interactive Plotly or static Matplotlib
            
        Returns:
            Union[go.Figure, plt.Figure]: Plotly or Matplotlib figure
        """
        if not pathway_scores:
            return None
        
        # Prepare pathway data
        pathways = list(pathway_scores.keys())
        scores = list(pathway_scores.values())
        
        # Sort by absolute score for better visualization
        sorted_indices = np.argsort(np.abs(scores))[::-1]
        pathways = [pathways[i] for i in sorted_indices]
        scores = [scores[i] for i in sorted_indices]
        
        # Create visualization
        if interactive:
            # Create a bar chart for pathways
            fig = go.Figure()
            
            # Add pathway bars
            colors = [self.colors["edges"]["activation"] if s > 0 else self.colors["edges"]["inhibition"] 
                    for s in scores]
            
            fig.add_trace(go.Bar(
                x=pathways,
                y=scores,
                marker_color=colors,
                name="Pathway Scores",
                hovertemplate="%{x}: %{y:.3f}"
            ))
            
            # Add horizontal line at y=0
            fig.add_shape(
                type="line",
                x0=-0.5,
                x1=len(pathways)-0.5,
                y0=0,
                y1=0,
                line=dict(color="gray", width=1)
            )
            
            # Update layout
            fig.update_layout(
                title=title,
                xaxis_title="Pathway",
                yaxis_title="Activation/Inhibition Score",
                xaxis_tickangle=45,
                template="plotly_white" if self.theme == "light" else "plotly_dark"
            )
            
            # Add gene scores in a separate subplot if available
            if gene_scores:
                # Sort genes by absolute score
                gene_items = sorted(gene_scores.items(), key=lambda x: abs(x[1]), reverse=True)
                gene_names = [item[0] for item in gene_items]
                gene_vals = [item[1] for item in gene_items]
                
                # Limit to top 20 genes
                if len(gene_names) > 20:
                    gene_names = gene_names[:20]
                    gene_vals = gene_vals[:20]
                
                # Create subplot with gene scores
                fig = make_subplots(
                    rows=2, 
                    cols=1, 
                    subplot_titles=("Pathway Scores", "Top Contributing Genes"),
                    row_heights=[0.6, 0.4]
                )
                
                # Add pathway bars to first subplot
                fig.add_trace(
                    go.Bar(
                        x=pathways,
                        y=scores,
                        marker_color=colors,
                        name="Pathway Scores",
                        hovertemplate="%{x}: %{y:.3f}"
                    ),
                    row=1, col=1
                )
                
                # Add gene bars to second subplot
                gene_colors = [self.colors["edges"]["activation"] if s > 0 else self.colors["edges"]["inhibition"] 
                             for s in gene_vals]
                
                fig.add_trace(
                    go.Bar(
                        x=gene_names,
                        y=gene_vals,
                        marker_color=gene_colors,
                        name="Gene Scores",
                        hovertemplate="%{x}: %{y:.3f}"
                    ),
                    row=2, col=1
                )
                
                # Add horizontal lines at y=0
                fig.add_shape(
                    type="line",
                    x0=-0.5,
                    x1=len(pathways)-0.5,
                    y0=0,
                    y1=0,
                    line=dict(color="gray", width=1),
                    row=1, col=1
                )
                
                fig.add_shape(
                    type="line",
                    x0=-0.5,
                    x1=len(gene_names)-0.5,
                    y0=0,
                    y1=0,
                    line=dict(color="gray", width=1),
                    row=2, col=1
                )
                
                # Update layout
                fig.update_layout(
                    title=title,
                    xaxis_tickangle=45,
                    xaxis2_tickangle=45,
                    template="plotly_white" if self.theme == "light" else "plotly_dark"
                )
                
                fig.update_xaxes(title_text="Pathway", row=1, col=1)
                fig.update_xaxes(title_text="Gene", row=2, col=1)
                fig.update_yaxes(title_text="Activation/Inhibition Score", row=1, col=1)
                fig.update_yaxes(title_text="Contribution Score", row=2, col=1)
            
            return fig
        else:
            # Create static visualization with Matplotlib
            if gene_scores:
                # Create figure with two subplots
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), height_ratios=[3, 2])
                
                # Plot pathway scores in first subplot
                colors = [self.colors["edges"]["activation"] if s > 0 else self.colors["edges"]["inhibition"] 
                        for s in scores]
                ax1.bar(pathways, scores, color=colors)
                ax1.axhline(y=0, color='gray', linestyle='-', alpha=0.7)
                ax1.set_title("Pathway Scores")
                ax1.set_xlabel("Pathway")
                ax1.set_ylabel("Activation/Inhibition Score")
                plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
                
                # Sort genes by absolute score
                gene_items = sorted(gene_scores.items(), key=lambda x: abs(x[1]), reverse=True)
                gene_names = [item[0] for item in gene_items]
                gene_vals = [item[1] for item in gene_items]
                
                # Limit to top 20 genes
                if len(gene_names) > 20:
                    gene_names = gene_names[:20]
                    gene_vals = gene_vals[:20]
                
                # Plot gene scores in second subplot
                gene_colors = [self.colors["edges"]["activation"] if s > 0 else self.colors["edges"]["inhibition"] 
                             for s in gene_vals]
                ax2.bar(gene_names, gene_vals, color=gene_colors)
                ax2.axhline(y=0, color='gray', linestyle='-', alpha=0.7)
                ax2.set_title("Top Contributing Genes")
                ax2.set_xlabel("Gene")
                ax2.set_ylabel("Contribution Score")
                plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
                
                plt.suptitle(title, fontsize=16)
                plt.tight_layout()
            else:
                # Create figure with single plot
                plt.figure(figsize=(12, 6))
                
                colors = [self.colors["edges"]["activation"] if s > 0 else self.colors["edges"]["inhibition"] 
                        for s in scores]
                plt.bar(pathways, scores, color=colors)
                plt.axhline(y=0, color='gray', linestyle='-', alpha=0.7)
                plt.title(title)
                plt.xlabel("Pathway")
                plt.ylabel("Activation/Inhibition Score")
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()
            
            return plt.gcf()
    
    def plot_pathway_enrichment(
        self, 
        pathway_df: pd.DataFrame,
        pathway_col: str = "pathway",
        pvalue_col: str = "pvalue",
        enrichment_col: Optional[str] = "enrichment",
        gene_count_col: Optional[str] = "gene_count",
        show_genes: bool = True,
        genes_col: Optional[str] = "genes",
        max_pathways: int = 15,
        title: str = "Pathway Enrichment Analysis",
        interactive: bool = True
    ) -> Union[go.Figure, plt.Figure]:
        """
        Plot pathway enrichment results.
        
        Args:
            pathway_df (pd.DataFrame): DataFrame with pathway enrichment results
            pathway_col (str): Column with pathway names
            pvalue_col (str): Column with p-values
            enrichment_col (str, optional): Column with enrichment scores
            gene_count_col (str, optional): Column with gene counts
            show_genes (bool): Whether to show genes in hover text
            genes_col (str, optional): Column with gene lists
            max_pathways (int): Maximum number of pathways to show
            title (str): Plot title
            interactive (bool): Whether to use interactive Plotly or static Matplotlib
            
        Returns:
            Union[go.Figure, plt.Figure]: Plotly or Matplotlib figure
        """
        if pathway_df.empty:
            return None
        
        # Check required columns
        if pathway_col not in pathway_df.columns or pvalue_col not in pathway_df.columns:
            return None
        
        # Sort by p-value
        sorted_df = pathway_df.sort_values(pvalue_col).reset_index(drop=True)
        
        # Limit to max pathways
        if len(sorted_df) > max_pathways:
            sorted_df = sorted_df.iloc[:max_pathways]
        
        # Create plot
        if interactive:
            # Create Plotly figure
            fig = go.Figure()
            
            # Calculate -log10(pvalue) for better visualization
            neg_log_p = -np.log10(sorted_df[pvalue_col])
            
            # Create hover text
            hover_text = []
            for i, row in sorted_df.iterrows():
                text = f"Pathway: {row[pathway_col]}<br>p-value: {row[pvalue_col]:.2e}"
                
                if enrichment_col in sorted_df.columns:
                    text += f"<br>Enrichment: {row[enrichment_col]:.2f}"
                    
                if gene_count_col in sorted_df.columns:
                    text += f"<br>Gene count: {row[gene_count_col]}"
                    
                if show_genes and genes_col in sorted_df.columns:
                    # Convert to list if it's a string
                    genes = row[genes_col]
                    if isinstance(genes, str):
                        if '[' in genes and ']' in genes:  # Check if it's a string representation of list
                            try:
                                genes = eval(genes)
                            except:
                                genes = genes.strip('[]').split(',')
                                
                    if isinstance(genes, list):
                        # Limit number of genes shown
                        if len(genes) > 10:
                            genes_text = ', '.join(genes[:10]) + ", ..."
                        else:
                            genes_text = ', '.join(genes)
                        text += f"<br>Genes: {genes_text}"
                
                hover_text.append(text)
            
            # Determine marker size
            if gene_count_col in sorted_df.columns:
                marker_size = sorted_df[gene_count_col] * 2 + 10
            elif enrichment_col in sorted_df.columns:
                marker_size = sorted_df[enrichment_col] * 5 + 10
            else:
                marker_size = 20
            
            # Add scatter plot
            fig.add_trace(go.Scatter(
                x=neg_log_p,
                y=sorted_df[pathway_col],
                mode='markers',
                marker=dict(
                    size=marker_size,
                    color=neg_log_p,
                    colorscale='Viridis',
                    colorbar=dict(title="-log10(p-value)"),
                    showscale=True
                ),
                text=hover_text,
                hoverinfo='text'
            ))
            
            # Update layout
            fig.update_layout(
                title=title,
                xaxis_title="-log10(p-value)",
                yaxis_title="Pathway",
                template="plotly_white" if self.theme == "light" else "plotly_dark"
            )
            
            return fig
        else:
            # Create Matplotlib figure
            plt.figure(figsize=(10, max(6, max_pathways * 0.4)))
            
            # Calculate -log10(pvalue) for better visualization
            neg_log_p = -np.log10(sorted_df[pvalue_col])
            
            # Determine marker size
            if gene_count_col in sorted_df.columns:
                marker_size = sorted_df[gene_count_col] * 2 + 10
            elif enrichment_col in sorted_df.columns:
                marker_size = sorted_df[enrichment_col] * 5 + 10
            else:
                marker_size = 20
            
            # Create colormap
            cmap = plt.cm.viridis
            
            # Plot pathways
            plt.scatter(
                neg_log_p,
                range(len(sorted_df)),
                s=marker_size,
                c=neg_log_p,
                cmap=cmap,
                alpha=0.8
            )
            
            # Add pathway labels
            plt.yticks(range(len(sorted_df)), sorted_df[pathway_col])
            
            # Add colorbar
            plt.colorbar(label="-log10(p-value)")
            
            # Add labels and title
            plt.xlabel("-log10(p-value)")
            plt.title(title)
            
            plt.grid(axis='x', linestyle='--', alpha=0.7)
            plt.tight_layout()
            
            return plt.gcf()
    
    def plot_gene_interaction_heatmap(
        self,
        interaction_matrix: Union[pd.DataFrame, np.ndarray],
        gene_names: Optional[List[str]] = None,
        title: str = "Gene Interaction Heatmap",
        symmetric: bool = True,
        interactive: bool = True
    ) -> Union[go.Figure, plt.Figure]:
        """
        Plot a heatmap of gene-gene interactions.
        
        Args:
            interaction_matrix (Union[pd.DataFrame, np.ndarray]): Interaction matrix
            gene_names (List[str], optional): List of gene names
            title (str): Plot title
            symmetric (bool): Whether the matrix is symmetric
            interactive (bool): Whether to use interactive Plotly or static Matplotlib
            
        Returns:
            Union[go.Figure, plt.Figure]: Plotly or Matplotlib figure
        """
        # Convert to DataFrame if numpy array
        if isinstance(interaction_matrix, np.ndarray):
            if gene_names is None:
                gene_names = [f"Gene {i+1}" for i in range(interaction_matrix.shape[0])]
            
            interaction_df = pd.DataFrame(
                interaction_matrix,
                index=gene_names,
                columns=gene_names
            )
        else:
            interaction_df = interaction_matrix
        
        # For interactive plots
        if interactive:
            # Create a custom colorscale for interactions
            # Blue for negative, red for positive, white for zero
            colorscale = [
                [0, "blue"],               # Most negative values
                [0.45, "lightblue"],       # Slightly negative
                [0.5, "white"],            # Zero or near-zero
                [0.55, "pink"],            # Slightly positive
                [1, "red"]                 # Most positive values
            ]
            
            # Create heatmap
            fig = px.imshow(
                interaction_df,
                color_continuous_scale=colorscale,
                labels=dict(x="Gene", y="Gene", color="Interaction"),
                title=title
            )
            
            fig.update_layout(
                width=800,
                height=800,
                template="plotly_white" if self.theme == "light" else "plotly_dark"
            )
            
            return fig
        else:
            # Create Matplotlib figure
            plt.figure(figsize=(10, 10))
            
            # Create a custom colormap
            colors = ["blue", "lightblue", "white", "pink", "red"]
            n_bins = 100
            cmap_name = "interaction_cmap"
            cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)
            
            # Plot heatmap
            sns.heatmap(
                interaction_df,
                cmap=cmap,
                center=0,
                annot=interaction_df.shape[0] < 15,  # Only annotate if small enough
                fmt=".2f" if interaction_df.shape[0] < 15 else "",
                linewidths=0.5,
                square=True
            )
            
            plt.title(title)
            plt.tight_layout()
            
            return plt.gcf()
    
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