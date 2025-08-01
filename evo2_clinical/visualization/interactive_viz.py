"""
Interactive visualization components for Evo2_Clinical.

This module provides interactive visualization tools and dashboard components
for exploring variant data, lncRNA data, and analysis results.
"""

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np


def create_interactive_scatterplot(data_df, x_col, y_col, color_col=None, size_col=None,
                                hover_cols=None, title=None):
    """
    Create an interactive scatter plot with hover information and customizable appearance.
    
    Parameters:
    -----------
    data_df : pandas.DataFrame
        DataFrame containing the data to visualize
    x_col : str
        Column to use for x-axis
    y_col : str
        Column to use for y-axis
    color_col : str, optional
        Column to use for point colors
    size_col : str, optional
        Column to use for point sizes
    hover_cols : list, optional
        List of columns to include in hover information
    title : str, optional
        Plot title
        
    Returns:
    --------
    plotly.graph_objects.Figure
        Interactive plotly figure
    """
    if title is None:
        title = f"{y_col} vs {x_col}"
        
    # Prepare hover data
    if hover_cols is None:
        hover_cols = []
        
    # Create scatter plot
    fig = px.scatter(
        data_df,
        x=x_col,
        y=y_col,
        color=color_col,
        size=size_col,
        hover_data=hover_cols,
        title=title
    )
    
    # Update layout for better interactivity
    fig.update_layout(
        hovermode='closest',
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    
    return fig


def create_genomic_browser(variants_df, chr_col='#CHROM', pos_col='POS', 
                        ref_col='REF', alt_col='ALT', info_col='INFO', 
                        score_col=None, region=None):
    """
    Create an interactive genomic browser visualization.
    
    Parameters:
    -----------
    variants_df : pandas.DataFrame
        DataFrame containing variant data
    chr_col : str, optional
        Column containing chromosome information
    pos_col : str, optional
        Column containing position information
    ref_col : str, optional
        Column containing reference allele
    alt_col : str, optional
        Column containing alternate allele
    info_col : str, optional
        Column containing variant information
    score_col : str, optional
        Column containing a score to use for variant color
    region : tuple, optional
        Tuple of (chromosome, start, end) to restrict view
        
    Returns:
    --------
    plotly.graph_objects.Figure
        Interactive genomic browser visualization
    """
    # Filter by region if specified
    if region:
        chromosome, start_pos, end_pos = region
        filtered_df = variants_df[(variants_df[chr_col] == chromosome) & 
                                (variants_df[pos_col] >= start_pos) & 
                                (variants_df[pos_col] <= end_pos)]
    else:
        filtered_df = variants_df.copy()
    
    # Check if we have data
    if filtered_df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No variants in the selected region",
            showarrow=False,
            xref="paper",
            yref="paper",
            x=0.5,
            y=0.5
        )
        return fig
    
    # Sort by position
    filtered_df = filtered_df.sort_values(pos_col)
    
    # Prepare tooltip information
    hover_text = []
    for i, row in filtered_df.iterrows():
        hover_info = f"Pos: {row[pos_col]}<br>"
        hover_info += f"Variant: {row[ref_col]}>{row[alt_col]}<br>"
        if info_col in row and row[info_col]:
            hover_info += f"Info: {row[info_col]}<br>"
        if score_col and score_col in row:
            hover_info += f"Score: {row[score_col]:.4f}"
        hover_text.append(hover_info)
    
    # Create the figure
    if score_col and score_col in filtered_df.columns:
        fig = px.scatter(
            filtered_df,
            x=pos_col,
            y=[0.5] * len(filtered_df),  # Place all points on a single line
            color=score_col,
            size=[10] * len(filtered_df),  # Fixed size
            hover_name=[f"{r[ref_col]}>{r[alt_col]}" for _, r in filtered_df.iterrows()],
            hover_data={pos_col: True},
            color_continuous_scale='Viridis',
            title=f"Genomic View: {filtered_df[chr_col].iloc[0]} - {filtered_df[pos_col].min():,} to {filtered_df[pos_col].max():,}"
        )
    else:
        fig = px.scatter(
            filtered_df,
            x=pos_col,
            y=[0.5] * len(filtered_df),  # Place all points on a single line
            size=[10] * len(filtered_df),  # Fixed size
            hover_name=[f"{r[ref_col]}>{r[alt_col]}" for _, r in filtered_df.iterrows()],
            hover_data={pos_col: True},
            title=f"Genomic View: {filtered_df[chr_col].iloc[0]} - {filtered_df[pos_col].min():,} to {filtered_df[pos_col].max():,}"
        )
    
    # Customize layout
    fig.update_layout(
        yaxis=dict(
            visible=False,
            showticklabels=False,
            fixedrange=True
        ),
        xaxis=dict(
            title="Genomic Position",
            showgrid=True
        ),
        hoverlabel=dict(
            bgcolor="white",
            font_size=12
        )
    )
    
    return fig


def create_interactive_heatmap(data_matrix, x_labels=None, y_labels=None, 
                            title="Interactive Heatmap", colorscale='Viridis'):
    """
    Create an interactive heatmap visualization.
    
    Parameters:
    -----------
    data_matrix : numpy.ndarray or pandas.DataFrame
        2D data matrix to visualize
    x_labels : list, optional
        Labels for the x-axis
    y_labels : list, optional
        Labels for the y-axis
    title : str, optional
        Plot title
    colorscale : str, optional
        Colorscale to use
        
    Returns:
    --------
    plotly.graph_objects.Figure
        Interactive heatmap visualization
    """
    # Convert to numpy array if DataFrame
    if isinstance(data_matrix, pd.DataFrame):
        if x_labels is None:
            x_labels = data_matrix.columns.tolist()
        if y_labels is None:
            y_labels = data_matrix.index.tolist()
        data_matrix = data_matrix.values
    
    # Create heatmap figure
    fig = go.Figure(data=go.Heatmap(
        z=data_matrix,
        x=x_labels,
        y=y_labels,
        colorscale=colorscale,
        colorbar=dict(title="Value")
    ))
    
    # Update layout
    fig.update_layout(
        title=title,
        xaxis=dict(title="", tickangle=45),
        yaxis=dict(title="")
    )
    
    return fig


def create_correlation_explorer(data_df, numeric_cols=None, default_x=None, default_y=None):
    """
    Create an interactive correlation explorer with multiple visualization types.
    
    Parameters:
    -----------
    data_df : pandas.DataFrame
        DataFrame containing the data to explore
    numeric_cols : list, optional
        List of numeric columns to include (if None, all numeric columns are used)
    default_x : str, optional
        Default column for x-axis
    default_y : str, optional
        Default column for y-axis
        
    Returns:
    --------
    dict
        Dictionary containing different visualization figures
    """
    # Get numeric columns if not provided
    if numeric_cols is None:
        numeric_cols = data_df.select_dtypes(include=[np.number]).columns.tolist()
    
    if len(numeric_cols) < 2:
        return {"error": "Not enough numeric columns for correlation analysis"}
    
    # Set defaults if not provided
    if default_x is None:
        default_x = numeric_cols[0]
    if default_y is None:
        default_y = numeric_cols[1] if len(numeric_cols) > 1 else numeric_cols[0]
    
    result = {}
    
    # Create correlation matrix
    corr_matrix = data_df[numeric_cols].corr()
    
    # Correlation heatmap
    heatmap_fig = create_interactive_heatmap(
        corr_matrix,
        title="Correlation Matrix"
    )
    result['heatmap'] = heatmap_fig
    
    # Scatter plot for selected variables
    scatter_fig = create_interactive_scatterplot(
        data_df,
        x_col=default_x,
        y_col=default_y,
        title=f"Scatter Plot: {default_y} vs {default_x}"
    )
    result['scatter'] = scatter_fig
    
    # Pair plot for multiple columns
    if len(numeric_cols) > 2:
        # Create a sample of the data if it's too large
        sample_size = min(1000, len(data_df))
        sample_df = data_df.sample(sample_size) if len(data_df) > sample_size else data_df
        
        # Create pair plot matrix (up to 5 columns for clarity)
        display_cols = numeric_cols[:5]
        n_cols = len(display_cols)
        
        # Create subplot grid
        fig = make_subplots(rows=n_cols, cols=n_cols,
                          shared_xaxes=False, shared_yaxes=False,
                          vertical_spacing=0.05, horizontal_spacing=0.05)
        
        # Add scatter plots
        for i, col_i in enumerate(display_cols):
            for j, col_j in enumerate(display_cols):
                if i == j:
                    # Histogram on diagonal
                    fig.add_trace(
                        go.Histogram(x=sample_df[col_i], name=col_i),
                        row=i+1, col=j+1
                    )
                else:
                    # Scatter plot on off-diagonal
                    fig.add_trace(
                        go.Scatter(
                            x=sample_df[col_j],
                            y=sample_df[col_i],
                            mode='markers',
                            marker=dict(size=3),
                            showlegend=False
                        ),
                        row=i+1, col=j+1
                    )
        
        # Update layout
        fig.update_layout(
            title="Pair Plot Matrix",
            height=250 * n_cols,
            width=250 * n_cols,
            showlegend=False
        )
        
        result['pair_plot'] = fig
    
    return result


def create_interactive_lncrna_explorer(lncrna_df, gene_variants_df=None, expression_data=None):
    """
    Create an interactive explorer for lncRNA data with multiple linked visualizations.
    
    Parameters:
    -----------
    lncrna_df : pandas.DataFrame
        DataFrame containing lncRNA data
    gene_variants_df : pandas.DataFrame, optional
        DataFrame containing variant data for genes
    expression_data : pandas.DataFrame, optional
        DataFrame containing expression data
        
    Returns:
    --------
    dict
        Dictionary containing different visualization figures
    """
    result = {}
    
    # Create basic bar chart for lncRNA metrics
    if 'length' in lncrna_df.columns:
        length_fig = px.bar(
            lncrna_df.sort_values('length', ascending=False),
            x='lncrna',
            y='length',
            title='lncRNA Length Distribution',
            labels={'lncrna': 'lncRNA', 'length': 'Length (bp)'}
        )
        result['length'] = length_fig
    
    # Create scatter plot if functionality score exists
    if 'functionality_score' in lncrna_df.columns:
        if 'length' in lncrna_df.columns:
            func_fig = px.scatter(
                lncrna_df,
                x='length',
                y='functionality_score',
                hover_name='lncrna',
                title='lncRNA Functionality Score vs Length',
                labels={'functionality_score': 'Functionality Score', 'length': 'Length (bp)'}
            )
            result['functionality'] = func_fig
    
    # Create radar chart if multiple score dimensions exist
    score_columns = [col for col in lncrna_df.columns if col.endswith('_score')]
    if len(score_columns) >= 3 and 'lncrna' in lncrna_df.columns:
        # Create radar chart data
        fig = go.Figure()
        
        for i, row in lncrna_df.iterrows():
            fig.add_trace(go.Scatterpolar(
                r=[row[col] for col in score_columns],
                theta=score_columns,
                fill='toself',
                name=row['lncrna']
            ))
        
        fig.update_layout(
            polar=dict(
                radialaxis=dict(
                    visible=True,
                    range=[0, 1]
                )),
            title='Multi-dimensional lncRNA Scores',
            showlegend=True
        )
        
        result['radar'] = fig
    
    # Create variant count visualization if gene variant data available
    if gene_variants_df is not None and 'lncrna' in lncrna_df.columns:
        # Extract gene information
        if 'gene' in gene_variants_df.columns:
            # Count variants per gene
            gene_counts = gene_variants_df['gene'].value_counts().reset_index()
            gene_counts.columns = ['gene', 'variant_count']
            
            # Merge with lncRNA data if possible
            if 'associated_gene' in lncrna_df.columns:
                merged_data = pd.merge(
                    lncrna_df, 
                    gene_counts, 
                    left_on='associated_gene', 
                    right_on='gene', 
                    how='inner'
                )
                
                if not merged_data.empty:
                    fig = px.bar(
                        merged_data,
                        x='lncrna',
                        y='variant_count',
                        title='Variant Counts in Genes Associated with lncRNAs',
                        labels={'lncrna': 'lncRNA', 'variant_count': 'Number of Variants'}
                    )
                    result['variant_counts'] = fig
    
    # Create expression heatmap if expression data available
    if expression_data is not None and 'lncrna' in lncrna_df.columns:
        lncrna_names = lncrna_df['lncrna'].tolist()
        matched_expression = expression_data[expression_data['gene_name'].isin(lncrna_names)]
        
        if not matched_expression.empty:
            # Pivot data for heatmap
            sample_cols = [col for col in matched_expression.columns if col.startswith('sample_')]
            if sample_cols:
                pivot_df = matched_expression.pivot(index='gene_name', columns=sample_cols)
                
                # Create interactive heatmap
                heatmap_fig = create_interactive_heatmap(
                    pivot_df,
                    title='lncRNA Expression Heatmap'
                )
                result['expression'] = heatmap_fig
    
    return result


def create_interactive_pathway_explorer(gene_list, interaction_data, expression_data=None):
    """
    Create an interactive pathway explorer visualization.
    
    Parameters:
    -----------
    gene_list : list
        List of genes to include
    interaction_data : pandas.DataFrame
        DataFrame containing gene interactions with 'source', 'target', and 'score' columns
    expression_data : pandas.DataFrame, optional
        DataFrame containing expression data with 'gene_name' and expression columns
        
    Returns:
    --------
    plotly.graph_objects.Figure
        Interactive pathway visualization
    """
    # Filter interactions to only include genes in our list
    filtered_interactions = interaction_data[
        (interaction_data['source'].isin(gene_list)) & 
        (interaction_data['target'].isin(gene_list))
    ]
    
    # Get unique genes for nodes
    unique_genes = set(filtered_interactions['source']).union(set(filtered_interactions['target']))
    
    # Check if we have expression data for node coloring
    expression_values = {}
    if expression_data is not None:
        gene_expression = expression_data[expression_data['gene_name'].isin(unique_genes)]
        if not gene_expression.empty and 'expression_value' in gene_expression.columns:
            expression_values = dict(zip(gene_expression['gene_name'], gene_expression['expression_value']))
    
    # Create node positions (using Fruchterman-Reingold-like algorithm)
    pos = {}
    np.random.seed(42)  # For reproducibility
    
    # Initial random positions
    for gene in unique_genes:
        pos[gene] = [np.random.uniform(-1, 1), np.random.uniform(-1, 1)]
    
    # Simple force-directed layout (limited iterations for demonstration)
    iterations = 50
    k = 0.1  # Optimal distance
    
    for _ in range(iterations):
        # Calculate repulsive forces
        for gene_i in unique_genes:
            disp = [0, 0]
            for gene_j in unique_genes:
                if gene_i != gene_j:
                    dx = pos[gene_i][0] - pos[gene_j][0]
                    dy = pos[gene_i][1] - pos[gene_j][1]
                    dist = max(0.01, np.sqrt(dx*dx + dy*dy))  # Avoid division by zero
                    
                    # Repulsive force
                    disp[0] += k*k/dist * dx/dist
                    disp[1] += k*k/dist * dy/dist
            
            # Update position based on repulsive force
            mag = max(0.01, np.sqrt(disp[0]*disp[0] + disp[1]*disp[1]))
            pos[gene_i][0] += disp[0]/mag * min(mag, 0.1)
            pos[gene_i][1] += disp[1]/mag * min(mag, 0.1)
    
    # Calculate attractive forces from edges
    for _, row in filtered_interactions.iterrows():
        source = row['source']
        target = row['target']
        score = row.get('score', 1.0)
        
        # Get positions
        dx = pos[target][0] - pos[source][0]
        dy = pos[target][1] - pos[source][1]
        dist = max(0.01, np.sqrt(dx*dx + dy*dy))
        
        # Attractive force (proportional to score)
        factor = dist*dist/k * min(score, 1.0)
        
        # Move points toward each other
        pos[source][0] += dx * factor * 0.05
        pos[source][1] += dy * factor * 0.05
        pos[target][0] -= dx * factor * 0.05
        pos[target][1] -= dy * factor * 0.05
    
    # Normalize positions to fit in the plot
    x_values = [p[0] for p in pos.values()]
    y_values = [p[1] for p in pos.values()]
    x_min, x_max = min(x_values), max(x_values)
    y_min, y_max = min(y_values), max(y_values)
    
    for gene in pos:
        pos[gene][0] = (pos[gene][0] - x_min) / (x_max - x_min) * 0.9 + 0.05
        pos[gene][1] = (pos[gene][1] - y_min) / (y_max - y_min) * 0.9 + 0.05
    
    # Create node trace
    node_x = [pos[gene][0] for gene in unique_genes]
    node_y = [pos[gene][1] for gene in unique_genes]
    
    # Set node colors based on expression if available
    if expression_values:
        node_color = [expression_values.get(gene, 0) for gene in unique_genes]
        colorscale = 'Viridis'
    else:
        node_color = ['lightblue'] * len(unique_genes)
        colorscale = None
    
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        hoverinfo='text',
        marker=dict(
            size=15,
            color=node_color,
            colorscale=colorscale,
            line=dict(width=2)
        ),
        text=list(unique_genes),
        textposition='top center'
    )
    
    # Create edge traces
    edge_traces = []
    
    # Dictionary to store interactions for creating edges
    edge_dict = {}
    
    for _, row in filtered_interactions.iterrows():
        source = row['source']
        target = row['target']
        
        # Get positions
        source_pos = pos[source]
        target_pos = pos[target]
        
        # Get score if available
        score = row.get('score', 1.0)
        
        # Create edge trace
        edge_trace = go.Scatter(
            x=[source_pos[0], target_pos[0], None],
            y=[source_pos[1], target_pos[1], None],
            mode='lines',
            line=dict(width=2, color=f'rgba(100, 100, 100, {min(score, 1.0)})'),
            hoverinfo='text',
            hovertext=f"{source} - {target} (Score: {score:.2f})"
        )
        edge_traces.append(edge_trace)
    
    # Create figure
    fig = go.Figure(data=edge_traces + [node_trace],
                   layout=go.Layout(
                       title='Interactive Pathway Network',
                       showlegend=False,
                       hovermode='closest',
                       margin=dict(b=20, l=5, r=5, t=40),
                       xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       width=800,
                       height=800
                   ))
    
    return fig


def create_variant_timeline_visualization(variants_df, time_col, value_col, group_col=None):
    """
    Create an interactive timeline visualization of variant data.
    
    Parameters:
    -----------
    variants_df : pandas.DataFrame
        DataFrame containing variant data with time information
    time_col : str
        Column containing time/date information
    value_col : str
        Column containing values to plot
    group_col : str, optional
        Column to use for grouping variants
        
    Returns:
    --------
    plotly.graph_objects.Figure
        Interactive timeline visualization
    """
    # Create base figure
    if group_col and group_col in variants_df.columns:
        fig = px.line(
            variants_df,
            x=time_col,
            y=value_col,
            color=group_col,
            title=f'Variant Timeline: {value_col} Over {time_col}',
            labels={time_col: 'Time', value_col: 'Value'}
        )
    else:
        fig = px.line(
            variants_df,
            x=time_col,
            y=value_col,
            title=f'Variant Timeline: {value_col} Over {time_col}',
            labels={time_col: 'Time', value_col: 'Value'}
        )
    
    # Add range slider
    fig.update_layout(
        xaxis=dict(
            rangeslider=dict(
                visible=True
            )
        )
    )
    
    return fig