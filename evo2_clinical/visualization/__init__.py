"""
Visualization module for Evo2_Clinical.

This module provides visualization components for genetic variants, pathways,
treatment outcomes, and other analysis results.
"""

from .variant_viz import VariantVisualizer
from .pathway_viz import PathwayVisualizer
from .treatment_viz import TreatmentOutcomeVisualizer
from .dashboard import create_dashboard
from .advanced_plots import (
    plot_variant_distribution, 
    plot_variant_heatmap, 
    create_pathway_network,
    plot_lncrna_multidimensional,
    plot_epigenetic_correlation,
    create_manhattan_plot,
    variant_effect_visualization,
    opentargets_visualization
)
from .interactive_viz import (
    create_interactive_scatterplot,
    create_genomic_browser,
    create_interactive_heatmap,
    create_correlation_explorer,
    create_interactive_lncrna_explorer,
    create_interactive_pathway_explorer,
    create_variant_timeline_visualization
)

# Explicitly import modules to make them available
from . import advanced_plots
from . import interactive_viz

__all__ = [
    # Basic visualizers
    'VariantVisualizer',
    'PathwayVisualizer',
    'TreatmentOutcomeVisualizer',
    'create_dashboard',
    
    # Advanced plots
    'plot_variant_distribution', 
    'plot_variant_heatmap',
    'create_pathway_network',
    'plot_lncrna_multidimensional',
    'plot_epigenetic_correlation',
    'create_manhattan_plot',
    'variant_effect_visualization',
    'opentargets_visualization',
    
    # Interactive visualizations
    'create_interactive_scatterplot',
    'create_genomic_browser',
    'create_interactive_heatmap',
    'create_correlation_explorer',
    'create_interactive_lncrna_explorer',
    'create_interactive_pathway_explorer',
    'create_variant_timeline_visualization',
    
    # Module imports for backward compatibility
    'advanced_plots',
    'interactive_viz'
]