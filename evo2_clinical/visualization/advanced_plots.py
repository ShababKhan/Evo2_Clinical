"""
Advanced plotting functionality for the Evo2_Clinical package.

This module provides sophisticated visualization capabilities beyond the basic
visualizers, including interactive plots, multi-dimensional analyses, and 
special purpose visualizations for epigenetic and GWAS data.
"""

import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_variant_distribution(variants_df, column, title=None, color_column=None, nbins=20):
    """
    Create a histogram showing the distribution of values in the specified column.
    
    Parameters:
    -----------
    variants_df : pandas.DataFrame
        DataFrame containing variant data
    column : str
        The column to plot
    title : str, optional
        Plot title
    color_column : str, optional
        Column to use for color grouping
    nbins : int, optional
        Number of bins for histogram
        
    Returns:
    --------
    plotly.graph_objects.Figure
        The plotly figure object
    """
    if title is None:
        title = f"Distribution of {column}"
    
    fig = px.histogram(variants_df, x=column, color=color_column, 
                      nbins=nbins, title=title)
    return fig


def plot_variant_heatmap(variants_df, x_col, y_col, z_col=None, title=None):
    """
    Create a heatmap of variant data.
    
    Parameters:
    -----------
    variants_df : pandas.DataFrame
        DataFrame containing variant data
    x_col : str
        Column for x-axis
    y_col : str
        Column for y-axis
    z_col : str, optional
        Column for z-axis (intensity). If None, counts occurrences.
    title : str, optional
        Plot title
        
    Returns:
    --------
    plotly.graph_objects.Figure
        The plotly figure object
    """
    if title is None:
        title = f"{z_col if z_col else 'Count'} by {x_col} and {y_col}"
        
    if z_col is None:
        # Create contingency table
        contingency = pd.crosstab(variants_df[y_col], variants_df[x_col])
        fig = px.imshow(contingency, labels=dict(x=x_col, y=y_col, color="Count"),
                      title=title)
    else:
        # Aggregate data
        pivot = variants_df.pivot_table(index=y_col, columns=x_col, values=z_col, aggfunc='mean')
        fig = px.imshow(pivot, labels=dict(x=x_col, y=y_col, color=z_col),
                      title=title)
        
    return fig


def create_pathway_network(gene_list, interactions_df, score_column=None):
    """
    Create a network visualization of gene pathways.
    
    Parameters:
    -----------
    gene_list : list
        List of genes to include
    interactions_df : pandas.DataFrame
        DataFrame containing gene interaction data with columns 'source', 'target', and optional score
    score_column : str, optional
        Column containing interaction scores
        
    Returns:
    --------
    plotly.graph_objects.Figure
        The plotly figure object containing the network
    """
    # Filter interactions to only include genes in our list
    filtered_interactions = interactions_df[
        (interactions_df['source'].isin(gene_list)) & 
        (interactions_df['target'].isin(gene_list))
    ]
    
    # Get unique genes for nodes
    unique_genes = set(filtered_interactions['source']).union(set(filtered_interactions['target']))
    
    # Create node positions (using circular layout)
    n_genes = len(unique_genes)
    angles = np.linspace(0, 2*np.pi, n_genes, endpoint=False)
    
    node_x = [np.cos(angle) for angle in angles]
    node_y = [np.sin(angle) for angle in angles]
    
    # Create node trace
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        hoverinfo='text',
        marker=dict(
            size=15,
            color='lightblue',
            line=dict(width=2)
        ),
        text=list(unique_genes),
        textposition='top center'
    )
    
    # Create edge traces
    edge_traces = []
    
    # Mapping of genes to their positions
    gene_to_pos = {gene: (node_x[i], node_y[i]) for i, gene in enumerate(unique_genes)}
    
    for _, row in filtered_interactions.iterrows():
        source = row['source']
        target = row['target']
        
        # Get positions
        source_pos = gene_to_pos[source]
        target_pos = gene_to_pos[target]
        
        # Get score if available
        score = row[score_column] if score_column and score_column in row else 1.0
        
        # Create edge trace
        edge_trace = go.Scatter(
            x=[source_pos[0], target_pos[0], None],
            y=[source_pos[1], target_pos[1], None],
            mode='lines',
            line=dict(width=2, color=f'rgba(100, 100, 100, {min(score, 1.0)})'),
            hoverinfo='none'
        )
        edge_traces.append(edge_trace)
    
    # Create figure
    fig = go.Figure(data=edge_traces + [node_trace],
                   layout=go.Layout(
                       title='Pathway Interaction Network',
                       showlegend=False,
                       hovermode='closest',
                       margin=dict(b=20, l=5, r=5, t=40),
                       xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
                   ))
    
    return fig


def plot_lncrna_multidimensional(lncrna_df, dimensions, lncrna_col='lncrna'):
    """
    Create a radar chart for multi-dimensional analysis of lncRNAs.
    
    Parameters:
    -----------
    lncrna_df : pandas.DataFrame
        DataFrame containing lncRNA data
    dimensions : list
        List of column names to include in radar chart
    lncrna_col : str, optional
        Column containing lncRNA names
        
    Returns:
    --------
    plotly.graph_objects.Figure
        The plotly figure object
    """
    fig = go.Figure()
    
    for i, row in lncrna_df.iterrows():
        fig.add_trace(go.Scatterpolar(
            r=[row[dim] for dim in dimensions],
            theta=dimensions,
            fill='toself',
            name=row[lncrna_col]
        ))
    
    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 1]
            )),
        title='Multidimensional Analysis of lncRNAs',
        showlegend=True
    )
    
    return fig


def plot_epigenetic_correlation(variants_df, epigenetic_df, variant_id_col, epigenetic_id_col, 
                             score_col, epigenetic_marker_col, title=None):
    """
    Create a visualization of correlations between variants and epigenetic markers.
    
    Parameters:
    -----------
    variants_df : pandas.DataFrame
        DataFrame containing variant data
    epigenetic_df : pandas.DataFrame
        DataFrame containing epigenetic data
    variant_id_col : str
        Column in variants_df to use for joining
    epigenetic_id_col : str
        Column in epigenetic_df to use for joining
    score_col : str
        Column containing correlation scores
    epigenetic_marker_col : str
        Column containing epigenetic marker names
    title : str, optional
        Plot title
        
    Returns:
    --------
    plotly.graph_objects.Figure
        The plotly figure object
    """
    if title is None:
        title = f"Variant-Epigenetic Marker Correlation"
    
    # Join datasets
    merged_df = pd.merge(variants_df, epigenetic_df, 
                        left_on=variant_id_col, 
                        right_on=epigenetic_id_col)
    
    # Pivot table for heatmap
    pivot_df = merged_df.pivot_table(index=variant_id_col, 
                                   columns=epigenetic_marker_col, 
                                   values=score_col)
    
    # Create heatmap
    fig = px.imshow(pivot_df, 
                  title=title,
                  labels=dict(x="Epigenetic Marker", y="Variant", color="Correlation"))
    
    return fig


def create_manhattan_plot(gwas_df, chr_col, pos_col, pval_col, gene_col=None, 
                      highlight_genes=None, title="Manhattan Plot"):
    """
    Create a Manhattan plot for GWAS results.
    
    Parameters:
    -----------
    gwas_df : pandas.DataFrame
        DataFrame containing GWAS results
    chr_col : str
        Column containing chromosome information
    pos_col : str
        Column containing position information
    pval_col : str
        Column containing p-values
    gene_col : str, optional
        Column containing gene names
    highlight_genes : list, optional
        List of genes to highlight
    title : str, optional
        Plot title
        
    Returns:
    --------
    plotly.graph_objects.Figure
        The plotly figure object
    """
    # Create a copy of the data
    data = gwas_df.copy()
    
    # Convert chromosome to numeric for plotting
    if data[chr_col].dtype == 'object':
        data[chr_col] = data[chr_col].str.replace('chr', '').replace('X', '23').replace('Y', '24')
        data[chr_col] = pd.to_numeric(data[chr_col], errors='coerce')
    
    # Calculate -log10 of p-values
    data['-log10p'] = -np.log10(data[pval_col])
    
    # Sort by chromosome and position
    data = data.sort_values([chr_col, pos_col])
    
    # Get chromosome sizes
    chromosome_sizes = data.groupby(chr_col)[pos_col].max()
    
    # Calculate cumulative position
    data['cumulative_pos'] = 0
    data['color'] = 'gray'
    
    # Prepare alternating colors for chromosomes
    chromosome_colors = {i: 'darkblue' if i % 2 == 0 else 'lightblue' for i in range(1, 25)}
    
    # Calculate cumulative position and set colors
    cumulative_pos = 0
    for chrom in range(1, 25):
        if chrom in data[chr_col].values:
            data.loc[data[chr_col] == chrom, 'cumulative_pos'] = data.loc[data[chr_col] == chrom, pos_col] + cumulative_pos
            data.loc[data[chr_col] == chrom, 'color'] = chromosome_colors[chrom]
            cumulative_pos += chromosome_sizes.get(chrom, 0) + 1000000  # Add a gap between chromosomes
    
    # Highlight specified genes if applicable
    if highlight_genes and gene_col in data.columns:
        data.loc[data[gene_col].isin(highlight_genes), 'color'] = 'red'
    
    # Create Manhattan plot
    fig = px.scatter(
        data,
        x='cumulative_pos',
        y='-log10p',
        color='color',
        title=title,
        labels={'cumulative_pos': 'Chromosome Position', '-log10p': '-log10(p-value)'},
        color_discrete_map='identity'  # Use the color values directly
    )
    
    # Add horizontal line for significance threshold
    fig.add_shape(
        type='line',
        x0=data['cumulative_pos'].min(),
        x1=data['cumulative_pos'].max(),
        y0=5,  # -log10(1e-5)
        y1=5,
        line=dict(color='red', width=2, dash='dash')
    )
    
    # Hide legend and update layout
    fig.update_layout(showlegend=False)
    
    # Add chromosome labels
    chromosome_ticks = []
    for chrom in sorted(data[chr_col].unique()):
        subset = data[data[chr_col] == chrom]
        if not subset.empty:
            chromosome_ticks.append((subset['cumulative_pos'].min() + subset['cumulative_pos'].max()) / 2)
    
    fig.update_layout(
        xaxis=dict(
            tickvals=chromosome_ticks,
            ticktext=[str(int(c)) if c < 23 else 'X' if c == 23 else 'Y' for c in sorted(data[chr_col].unique())],
            title='Chromosome'
        )
    )
    
    return fig


def variant_effect_visualization(variants_df, effect_col, annotation_col=None):
    """
    Visualize variant effects by creating a pie chart and/or bar chart.
    
    Parameters:
    -----------
    variants_df : pandas.DataFrame
        DataFrame containing variant data
    effect_col : str
        Column containing variant effect information
    annotation_col : str, optional
        Additional column for grouping variants
        
    Returns:
    --------
    dict
        Dictionary containing plotly figure objects
    """
    result_figs = {}
    
    # Create pie chart of effects
    effect_counts = variants_df[effect_col].value_counts()
    pie_fig = px.pie(
        values=effect_counts.values,
        names=effect_counts.index,
        title='Distribution of Variant Effects'
    )
    result_figs['pie'] = pie_fig
    
    # Create grouped bar chart if annotation provided
    if annotation_col and annotation_col in variants_df.columns:
        grouped_counts = variants_df.groupby([effect_col, annotation_col]).size().reset_index(name='count')
        bar_fig = px.bar(
            grouped_counts,
            x=effect_col,
            y='count',
            color=annotation_col,
            title=f'Variant Effects by {annotation_col}',
            barmode='group'
        )
        result_figs['bar'] = bar_fig
    
    return result_figs


def opentargets_visualization(opentargets_data, score_column, source_column=None):
    """
    Create visualizations for OpenTargets data.
    
    Parameters:
    -----------
    opentargets_data : pandas.DataFrame
        DataFrame containing OpenTargets data
    score_column : str
        Column containing association scores
    source_column : str, optional
        Column containing data sources
        
    Returns:
    --------
    dict
        Dictionary containing plotly figure objects
    """
    result_figs = {}
    
    # Histogram of association scores
    hist_fig = px.histogram(
        opentargets_data,
        x=score_column,
        nbins=20,
        title='Distribution of Target-Disease Association Scores'
    )
    result_figs['histogram'] = hist_fig
    
    # Bar chart by source if available
    if source_column and source_column in opentargets_data.columns:
        source_scores = opentargets_data.groupby(source_column)[score_column].mean().reset_index()
        source_scores = source_scores.sort_values(score_column, ascending=False)
        
        bar_fig = px.bar(
            source_scores,
            x=source_column,
            y=score_column,
            title='Average Association Score by Source',
            labels={score_column: 'Mean Association Score'}
        )
        result_figs['bar'] = bar_fig
    
    return result_figs