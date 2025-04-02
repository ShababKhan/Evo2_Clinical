#!/usr/bin/env python3
"""
Evo2 Clinical Streamlit App

Interactive web interface for genomic variant analysis using Streamlit and PyGWalker.
Provides visualization and analysis tools for exploring genetic variants across genes.
"""

import streamlit as st
import pygwalker as pyg
import pandas as pd
from pathlib import Path
from typing import List, Optional
import plotly.express as px
from evo2_pipeline import Evo2Pipeline
from variant_analysis import VariantProcessor
from visualization import VariantVisualizer

# Configure page settings
st.set_page_config(
    page_title="Evo2 Clinical Analysis",
    page_icon="🧬",
    layout="wide"
)

# Initialize session state
if 'pipeline' not in st.session_state:
    st.session_state.pipeline = Evo2Pipeline({
        "data_dir": "data",
        "output_dir": "results"
    })

def load_variant_data(vcf_file: Optional[str] = None, chromosomes: List[str] = None) -> pd.DataFrame:
    """Load variant data from VCF file or specific chromosomes.
    
    Args:
        vcf_file: Optional VCF file uploaded by user
        chromosomes: List of chromosomes to analyze
        
    Returns:
        DataFrame containing variant data
    """
    processor = VariantProcessor()
    
    # If user uploaded a file, use that
    if vcf_file:
        return processor.load_vcf(vcf_file)
    
    # Otherwise load from chromosome-specific files
    elif chromosomes:
        return processor.load_multiple_chromosomes(chromosomes)
        
    # Return empty DataFrame if no data source specified
    return pd.DataFrame()

def main():
    st.title("🧬 Evo2 Clinical Analysis")
    
    # API Settings
    st.sidebar.header("NVIDIA Evo2 API Settings")
    use_nvidia_api = st.sidebar.checkbox(
        "Use NVIDIA Evo2 API", 
        value=True,
        help="Use NVIDIA Evo2 cloud service for variant scoring"
    )
    
    # Only show API key input if API is being used
    if use_nvidia_api:
        api_key = st.sidebar.text_input(
            "NVIDIA API Key", 
            type="password",
            help="Enter your NVIDIA Evo2 API key (leave empty to use NVIDIA_EVO2_API_KEY environment variable)"
        )
        
        if api_key:
            # Store API key in session state
            st.session_state.api_key = api_key
        
    # Analysis Settings    
    st.sidebar.header("Analysis Settings")

    # File upload section
    vcf_file = st.sidebar.file_uploader("Upload VCF File", type=["vcf", "vcf.gz"])
    
    # Chromosome selection
    available_chromosomes = [str(i) for i in range(1, 23)] + ["X", "Y"]
    selected_chromosomes = st.sidebar.multiselect(
        "Select Chromosomes",
        available_chromosomes,
        default=["19"],
        help="Select one or more chromosomes to analyze"
    )
    
    # Gene selection
    genes = st.sidebar.text_input(
        "Enter genes to analyze (comma-separated)",
        help="e.g., BRCA1, BRCA2, TP53"
    ).split(",")
    genes = [g.strip() for g in genes if g.strip()]

    # Cell type selection
    cell_types = st.sidebar.multiselect(
        "Select Cell Types",
        ["BREAST", "BLOOD", "LUNG", "LIVER", "KIDNEY", "BRAIN"],
        help="Select cell types for ENCODE data filtering"
    )

    # Disease/trait selection
    traits = st.sidebar.text_input(
        "Enter diseases/traits (comma-separated)",
        help="e.g., breast cancer, diabetes"
    ).split(",")
    traits = [t.strip() for t in traits if t.strip()]

    # Analysis type selection
    analysis_type = st.sidebar.selectbox(
        "Select Analysis Type",
        ["Variant Distribution", "Impact Scores", "Population Frequencies", "ENCODE Features"]
    )

    # Load data based on user input
    data_loaded = False
    variants_df = pd.DataFrame()
    
    if vcf_file:
        st.info("Loading and processing variant data from uploaded file...")
        variants_df = load_variant_data(vcf_file=vcf_file)
        data_loaded = not variants_df.empty
    elif selected_chromosomes:
        st.info(f"Loading and processing variant data from chromosomes: {', '.join(selected_chromosomes)}...")
        variants_df = load_variant_data(chromosomes=selected_chromosomes)
        data_loaded = not variants_df.empty
        
    if data_loaded:
        # Create tabs for different views
        tab1, tab2 = st.tabs(["PyGWalker Analysis", "Preset Visualizations"])
        
        with tab1:
            st.subheader("Interactive Data Analysis")
            pyg.walk(variants_df, env='Streamlit')
        
        with tab2:
            st.subheader("Preset Visualizations")
            
            if analysis_type == "Variant Distribution":
                fig = px.histogram(
                    variants_df,
                    x='pos',
                    color='gene',
                    title="Variant Distribution Across Genomic Positions"
                )
                st.plotly_chart(fig)
            
            elif analysis_type == "Impact Scores":
                if 'impact_score' in variants_df.columns:
                    fig = px.box(
                        variants_df,
                        x='gene',
                        y='impact_score',
                        title="Impact Score Distribution by Gene"
                    )
                    st.plotly_chart(fig)
            
            elif analysis_type == "Population Frequencies":
                if any(col.startswith('AF_') for col in variants_df.columns):
                    freq_cols = [col for col in variants_df.columns if col.startswith('AF_')]
                    for col in freq_cols:
                        fig = px.violin(
                            variants_df,
                            x='gene',
                            y=col,
                            title=f"Allele Frequency Distribution - {col}"
                        )
                        st.plotly_chart(fig)
            
            elif analysis_type == "ENCODE Features":
                encode_features = ['DNase', 'H3K27ac', 'H3K4me3']
                if all(feat in variants_df.columns for feat in encode_features):
                    fig = px.bar(
                        variants_df[encode_features].sum().reset_index(),
                        x='index',
                        y=0,
                        title="ENCODE Feature Counts"
                    )
                    st.plotly_chart(fig)
        
        # Add download button for processed data
        st.download_button(
            "Download Processed Data",
            variants_df.to_csv(index=False).encode('utf-8'),
            "processed_variants.csv",
            "text/csv",
            key='download-csv'
        )
    else:
        st.error("No variant data found in the uploaded file or selected chromosomes.")

    # Add documentation section
    with st.sidebar.expander("Documentation"):
        st.markdown("""
        ### How to Use
        1. Upload a VCF file containing variant data
        2. Enter genes of interest
        3. Select cell types and traits
        4. Choose analysis type
        5. Explore data using PyGWalker or preset visualizations
        
        ### Analysis Types
        - **Variant Distribution**: View variant positions across genes
        - **Impact Scores**: Analyze predicted functional impacts
        - **Population Frequencies**: Compare allele frequencies
        - **ENCODE Features**: Explore regulatory elements
        """)

if __name__ == "__main__":
    main()