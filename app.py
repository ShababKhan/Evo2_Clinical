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

def load_variant_data(vcf_file: Optional[str] = None) -> pd.DataFrame:
    """Load variant data from VCF file or use example data."""
    if vcf_file:
        processor = VariantProcessor()
        return processor.load_vcf(vcf_file)
    return pd.DataFrame()  # Return empty DataFrame if no file provided

def main():
    st.title("🧬 Evo2 Clinical Analysis")
    st.sidebar.header("Analysis Settings")

    # File upload section
    vcf_file = st.sidebar.file_uploader("Upload VCF File", type=["vcf", "vcf.gz"])
    
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

    # Load data when user provides input
    if vcf_file:
        st.info("Loading and processing variant data...")
        variants_df = load_variant_data(vcf_file)
        
        if not variants_df.empty:
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
            st.error("No variant data found in the uploaded file.")
    else:
        st.info("Please upload a VCF file to begin analysis.")

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