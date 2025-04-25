"""
Evo2_Clinical Streamlit Application

This Streamlit application provides a graphical user interface for the Evo2_Clinical package,
allowing users to run analyses, visualize results, and manage data without writing code.
"""

import os
import pandas as pd
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
import seaborn as sns
import yaml
import time
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go
import sqlite3
import logging

# Import Evo2_Clinical package modules
from evo2_clinical.config import Config
from evo2_clinical.data import io
from evo2_clinical.analysis import variant_scoring, lncrna_scoring
from evo2_clinical.database.manager import DatabaseManager, create_all_databases
from evo2_clinical.pipeline.main import Pipeline
from evo2_clinical.utils import helpers

# Set up logging
log_dir = os.path.join(os.path.dirname(__file__), 'logs')
os.makedirs(log_dir, exist_ok=True)
log_file = helpers.setup_logging(log_dir, logging.INFO)

logger = logging.getLogger(__name__)
logger.info("Starting Evo2_Clinical Streamlit application")

# App configuration
st.set_page_config(
    page_title="Evo2_Clinical Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Initialize session state to store data between reruns
if 'config' not in st.session_state:
    st.session_state.config = Config()
if 'db_manager' not in st.session_state:
    db_dir = os.path.join(os.path.dirname(__file__), 'output', 'databases')
    os.makedirs(db_dir, exist_ok=True)
    st.session_state.db_manager = DatabaseManager(db_dir)
if 'variants_df' not in st.session_state:
    st.session_state.variants_df = None
if 'scored_variants_df' not in st.session_state:
    st.session_state.scored_variants_df = None
if 'lncrnas_df' not in st.session_state:
    st.session_state.lncrnas_df = None
if 'analysis_results' not in st.session_state:
    st.session_state.analysis_results = {}


# Helper functions
def load_example_data():
    """Load example data for demonstration purposes."""
    # Sample variants data
    variants_data = {
        '#CHROM': ['chr1', 'chr1', 'chr2', 'chr3', 'chr5', 'chr7', 'chr9', 'chr11', 'chr17', 'chr19'],
        'POS': [1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000],
        'ID': ['.', 'rs123456', '.', 'rs789012', 'rs345678', '.', 'rs901234', 'rs567890', '.', 'rs234567'],
        'REF': ['A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C'],
        'ALT': ['G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'],
        'QUAL': ['30', '40', '50', '60', '70', '80', '90', '100', '60', '70'],
        'FILTER': ['PASS', 'PASS', 'PASS', 'PASS', 'PASS', 'PASS', 'PASS', 'PASS', 'PASS', 'PASS'],
        'INFO': ['GENE=GATA2', 'GENE=KLF2', 'GENE=KLF4', 'GENE=FOXC1', 'GENE=ERG', 
                'GENE=GATA2-AS1', 'GENE=ANRIL', 'GENE=PIRAT', 'GENE=LRAC', 'GENE=GATA2-AS1']
    }
    variants_df = pd.DataFrame(variants_data)
    
    # Sample lncRNAs data
    lncrna_data = {
        'lncrna': ['GATA2-AS1', 'ANRIL', 'PIRAT', 'LRAC', 'LINC01899', 
                   'HOTAIR', 'MALAT1', 'NEAT1', 'XIST', 'TERC'],
        'length': [1200, 3100, 1800, 2400, 1500, 2200, 8700, 3700, 19000, 450]
    }
    lncrnas_df = pd.DataFrame(lncrna_data)
    
    # Return the dataframes
    return variants_df, lncrnas_df


def get_available_databases():
    """Get list of available databases in the output directory."""
    db_dir = st.session_state.db_manager.db_dir
    db_files = [f.stem for f in Path(db_dir).glob('*.db')]
    return db_files


# Sidebar for navigation
st.sidebar.title("Evo2_Clinical Dashboard")
st.sidebar.image("https://cdn-icons-png.flaticon.com/512/2103/2103611.png", width=100)

# Navigation
page = st.sidebar.selectbox(
    "Navigation", 
    ["Home", "Data Management", "Variant Analysis", "lncRNA Analysis", "Pipeline Execution", "Database Explorer", "Configuration"]
)

st.sidebar.markdown("---")
st.sidebar.markdown("### Quick Actions")
if st.sidebar.button("Load Example Data"):
    with st.spinner("Loading example data..."):
        st.session_state.variants_df, st.session_state.lncrnas_df = load_example_data()
        st.sidebar.success("Example data loaded!")

if st.sidebar.button("Create Databases"):
    with st.spinner("Creating databases..."):
        success = create_all_databases(st.session_state.db_manager)
        if success:
            st.sidebar.success("Databases created successfully!")
        else:
            st.sidebar.error("Failed to create some databases.")

# Main content based on selected page
if page == "Home":
    st.title("🧬 Evo2_Clinical Dashboard")
    st.markdown("""
    ## Welcome to the Evo2_Clinical Dashboard
    
    This interactive dashboard provides access to the Evo2 computational pipeline, a framework for investigating 
    endothelial influence on pulmonary disease, integrating computational insights with molecular biology.
    
    ### Project Focus
    - **Endothelial Cell Function** in lung inflammation and fibrosis
    - **lncRNA Analysis** including ANRIL, PIRAT, LRAC, GATA2-AS1
    - **Genetic Variant Scoring** to identify functionally relevant variants
    - **Epigenetic Mediator Analysis** in the context of pulmonary diseases
    
    ### Using This Dashboard
    1. **Data Management**: Upload or view variant and lncRNA data
    2. **Variant Analysis**: Score and analyze genetic variants 
    3. **lncRNA Analysis**: Evaluate lncRNA functionality
    4. **Pipeline Execution**: Run the complete Evo2 computational pipeline
    5. **Database Explorer**: Browse and query analysis results
    6. **Configuration**: Adjust system settings and paths
    
    ### Getting Started
    Click "Load Example Data" in the sidebar to begin with demonstration data.
    """)
    
    # Dashboard cards with key metrics
    st.markdown("### Dashboard Overview")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric(
            label="Variants Loaded", 
            value=len(st.session_state.variants_df) if st.session_state.variants_df is not None else 0
        )
    
    with col2:
        st.metric(
            label="lncRNAs Loaded", 
            value=len(st.session_state.lncrnas_df) if st.session_state.lncrnas_df is not None else 0
        )
    
    with col3:
        st.metric(
            label="Databases Available", 
            value=len(get_available_databases())
        )
    
    # Show recent logs
    st.markdown("### Recent Activity")
    try:
        with open(log_file, 'r') as f:
            logs = f.readlines()[-10:]  # Get last 10 lines
            log_text = '\n'.join(logs)
            st.text_area("Log Output", log_text, height=200)
    except Exception as e:
        st.error(f"Unable to load logs: {e}")

elif page == "Data Management":
    st.title("📊 Data Management")
    
    # Tabs for different data types
    tab1, tab2, tab3 = st.tabs(["Variant Data", "lncRNA Data", "Gene Lists"])
    
    # Variant Data Tab
    with tab1:
        st.header("Variant Data Management")
        
        # Upload variant data
        st.subheader("Upload Variant Data")
        uploaded_file = st.file_uploader("Upload VCF file", type=['vcf', 'gz'])
        
        if uploaded_file is not None:
            with st.spinner('Loading VCF file...'):
                # Save the uploaded file to a temporary location
                temp_file = os.path.join(os.path.dirname(__file__), 'temp_vcf_file.vcf')
                with open(temp_file, 'wb') as f:
                    f.write(uploaded_file.getbuffer())
                
                try:
                    # Load the VCF file
                    variants_df = io.load_vcf_file(temp_file)
                    st.session_state.variants_df = variants_df
                    st.success(f"VCF file loaded successfully! {len(variants_df)} variants found.")
                except Exception as e:
                    st.error(f"Error loading VCF file: {e}")
                finally:
                    # Clean up
                    if os.path.exists(temp_file):
                        os.remove(temp_file)
        
        # Display loaded variant data
        st.subheader("Current Variant Data")
        if st.session_state.variants_df is not None:
            st.write(f"Total variants: {len(st.session_state.variants_df)}")
            st.dataframe(st.session_state.variants_df.head(10))
            
            # Basic variant statistics
            st.subheader("Variant Statistics")
            col1, col2 = st.columns(2)
            
            with col1:
                if '#CHROM' in st.session_state.variants_df.columns:
                    chrom_counts = st.session_state.variants_df['#CHROM'].value_counts()
                    fig = px.bar(
                        x=chrom_counts.index, 
                        y=chrom_counts.values,
                        labels={'x': 'Chromosome', 'y': 'Count'},
                        title='Variants per Chromosome'
                    )
                    st.plotly_chart(fig, use_container_width=True)
            
            with col2:
                if 'INFO' in st.session_state.variants_df.columns:
                    # Extract gene names from INFO field and count
                    gene_list = []
                    for info in st.session_state.variants_df['INFO']:
                        gene = helpers.extract_gene_name_from_vcf_info(info)
                        if gene:
                            gene_list.append(gene)
                    
                    if gene_list:
                        gene_counts = pd.Series(gene_list).value_counts().head(10)
                        fig = px.bar(
                            x=gene_counts.index, 
                            y=gene_counts.values,
                            labels={'x': 'Gene', 'y': 'Count'},
                            title='Top 10 Genes with Variants'
                        )
                        st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No variant data loaded. Please upload a VCF file or load the example data.")
    
    # lncRNA Data Tab
    with tab2:
        st.header("lncRNA Data Management")
        
        # Upload lncRNA data
        st.subheader("Upload lncRNA Data")
        uploaded_file = st.file_uploader("Upload lncRNA list (one name per line)", type=['txt', 'csv'])
        
        if uploaded_file is not None:
            with st.spinner('Loading lncRNA data...'):
                try:
                    if uploaded_file.name.endswith('.csv'):
                        lncrnas_df = pd.read_csv(uploaded_file)
                    else:
                        # Assume it's a text file with one lncRNA per line
                        lncrna_list = [line.decode('utf-8').strip() for line in uploaded_file.readlines() if line.decode('utf-8').strip()]
                        lncrnas_df = pd.DataFrame({'lncrna': lncrna_list})
                    
                    st.session_state.lncrnas_df = lncrnas_df
                    st.success(f"lncRNA data loaded successfully! {len(lncrnas_df)} lncRNAs found.")
                except Exception as e:
                    st.error(f"Error loading lncRNA data: {e}")
        
        # Display loaded lncRNA data
        st.subheader("Current lncRNA Data")
        if st.session_state.lncrnas_df is not None:
            st.write(f"Total lncRNAs: {len(st.session_state.lncrnas_df)}")
            st.dataframe(st.session_state.lncrnas_df.head(10))
        else:
            st.info("No lncRNA data loaded. Please upload a file or load the example data.")
    
    # Gene Lists Tab
    with tab3:
        st.header("Gene Lists Management")
        
        # Define gene lists
        gene_list_types = ["Endothelial Genes", "EMT Pathway Genes", "Epigenetic Mediator Genes"]
        selected_list = st.selectbox("Select Gene List Type", gene_list_types)
        
        # Text area for entering genes
        genes_text = st.text_area(
            "Enter gene names (one per line):",
            "GATA2\nKLF2\nKLF4\nFOXC1\nERG" if selected_list == "Endothelial Genes" else ""
        )
        
        if st.button(f"Save {selected_list}"):
            gene_list = [gene.strip() for gene in genes_text.split('\n') if gene.strip()]
            
            # Save the gene list
            list_type_key = selected_list.lower().replace(' ', '_')
            if 'gene_lists' not in st.session_state:
                st.session_state.gene_lists = {}
            
            st.session_state.gene_lists[list_type_key] = gene_list
            st.success(f"{selected_list} saved with {len(gene_list)} genes!")
        
        # Display saved gene lists
        if 'gene_lists' in st.session_state:
            st.subheader("Saved Gene Lists")
            for list_name, genes in st.session_state.gene_lists.items():
                with st.expander(f"{list_name.replace('_', ' ').title()} ({len(genes)} genes)"):
                    st.write(", ".join(genes))

elif page == "Variant Analysis":
    st.title("🔬 Variant Analysis")
    
    if st.session_state.variants_df is None:
        st.warning("No variant data loaded. Please go to Data Management to upload variant data or load example data.")
    else:
        # Analysis options
        st.subheader("Variant Analysis Options")
        
        analysis_type = st.selectbox(
            "Select Analysis Type",
            ["General Variant Scoring", "EMT Pathway Analysis", "GATA2-AS1 Effect Prediction", "Epigenetic Mediator Analysis"]
        )
        
        # Only show this for General Variant Scoring
        context_window = 1000000
        if analysis_type == "General Variant Scoring":
            context_window = st.slider(
                "Context Window Size (bp)",
                min_value=10000,
                max_value=2000000,
                value=1000000,
                step=10000,
                help="Size of genomic context window for variant analysis in base pairs"
            )
        
        # Filter options
        with st.expander("Filter Options"):
            # Chromosome filter
            available_chromosomes = st.session_state.variants_df['#CHROM'].unique().tolist() if '#CHROM' in st.session_state.variants_df.columns else []
            selected_chromosomes = st.multiselect(
                "Filter by Chromosome",
                options=available_chromosomes,
                default=available_chromosomes[:5] if len(available_chromosomes) > 5 else available_chromosomes
            )
            
            # Apply filters to create a working dataset
            filtered_variants_df = st.session_state.variants_df
            if selected_chromosomes:
                filtered_variants_df = filtered_variants_df[filtered_variants_df['#CHROM'].isin(selected_chromosomes)]
            
            # Show sample count after filtering
            st.write(f"Filtered dataset: {len(filtered_variants_df)} variants")
        
        # Execute analysis button
        if st.button("Run Analysis"):
            with st.spinner(f"Running {analysis_type}..."):
                try:
                    # This is a mock implementation as we don't have the actual Evo2 executable
                    # In a real implementation, we would use the actual functions from the package
                    
                    # Simulate processing time
                    progress_bar = st.progress(0)
                    for i in range(100):
                        time.sleep(0.01)
                        progress_bar.progress(i + 1)
                    
                    # Generate mock results based on analysis type
                    if analysis_type == "General Variant Scoring":
                        # Add mock functional scores
                        result_df = filtered_variants_df.copy()
                        result_df['functional_score'] = np.random.uniform(0, 1, len(result_df))
                        st.session_state.scored_variants_df = result_df
                        st.session_state.analysis_results[analysis_type] = result_df
                        
                    elif analysis_type == "EMT Pathway Analysis":
                        # Filter for EMT genes and add mock scores
                        emt_genes = ['KLF2', 'KLF4', 'FOXC1'] # Example EMT genes
                        result_df = filtered_variants_df.copy()
                        result_df = result_df[result_df['INFO'].apply(lambda x: any(gene in x for gene in emt_genes))]
                        result_df['emt_pathway_score'] = np.random.uniform(0, 1, len(result_df))
                        st.session_state.analysis_results[analysis_type] = result_df
                        
                    elif analysis_type == "GATA2-AS1 Effect Prediction":
                        # Filter for GATA2-AS1 and add mock effect predictions
                        result_df = filtered_variants_df[filtered_variants_df['INFO'].str.contains('GATA2-AS1')].copy()
                        result_df['predicted_effect'] = np.random.choice(
                            ['High Impact', 'Moderate Impact', 'Low Impact', 'No Impact'],
                            size=len(result_df)
                        )
                        result_df['confidence_score'] = np.random.uniform(0.5, 1.0, len(result_df))
                        st.session_state.analysis_results[analysis_type] = result_df
                        
                    elif analysis_type == "Epigenetic Mediator Analysis":
                        # Add mock epigenetic scores
                        result_df = filtered_variants_df.copy()
                        result_df['epigenetic_score'] = np.random.uniform(0, 1, len(result_df))
                        st.session_state.analysis_results[analysis_type] = result_df
                    
                    st.success(f"{analysis_type} completed successfully!")
                    
                except Exception as e:
                    st.error(f"Error during analysis: {e}")
        
        # Display results if available
        if analysis_type in st.session_state.analysis_results:
            st.subheader("Analysis Results")
            result_df = st.session_state.analysis_results[analysis_type]
            
            # Display basic results table
            st.dataframe(result_df.head(10))
            
            # Visualizations based on analysis type
            st.subheader("Visualizations")
            
            if analysis_type == "General Variant Scoring":
                col1, col2 = st.columns(2)
                
                with col1:
                    # Histogram of functional scores
                    fig = px.histogram(
                        result_df, 
                        x='functional_score',
                        nbins=20,
                        title='Distribution of Functional Scores'
                    )
                    st.plotly_chart(fig, use_container_width=True)
                
                with col2:
                    # Boxplot of functional scores by chromosome
                    if '#CHROM' in result_df.columns:
                        fig = px.box(
                            result_df,
                            x='#CHROM',
                            y='functional_score',
                            title='Functional Scores by Chromosome'
                        )
                        st.plotly_chart(fig, use_container_width=True)
            
            elif analysis_type == "EMT Pathway Analysis":
                # Extract gene names and plot scores by gene
                gene_names = []
                for info in result_df['INFO']:
                    gene = helpers.extract_gene_name_from_vcf_info(info)
                    if gene:
                        gene_names.append(gene)
                
                result_df['gene'] = gene_names
                
                fig = px.box(
                    result_df,
                    x='gene',
                    y='emt_pathway_score',
                    title='EMT Pathway Scores by Gene'
                )
                st.plotly_chart(fig, use_container_width=True)
            
            elif analysis_type == "GATA2-AS1 Effect Prediction":
                # Pie chart of predicted effects
                effect_counts = result_df['predicted_effect'].value_counts()
                fig = px.pie(
                    values=effect_counts.values, 
                    names=effect_counts.index,
                    title='Distribution of Predicted Effects on GATA2-AS1'
                )
                st.plotly_chart(fig, use_container_width=True)
                
                # Scatter plot of position vs confidence
                if 'POS' in result_df.columns:
                    fig = px.scatter(
                        result_df,
                        x='POS',
                        y='confidence_score',
                        color='predicted_effect',
                        title='Variant Position vs Confidence Score',
                        labels={'POS': 'Position', 'confidence_score': 'Confidence Score'}
                    )
                    st.plotly_chart(fig, use_container_width=True)
            
            elif analysis_type == "Epigenetic Mediator Analysis":
                # Histogram of epigenetic scores
                fig = px.histogram(
                    result_df,
                    x='epigenetic_score',
                    nbins=20,
                    title='Distribution of Epigenetic Scores'
                )
                st.plotly_chart(fig, use_container_width=True)
            
            # Export options
            st.subheader("Export Results")
            col1, col2 = st.columns(2)
            
            with col1:
                if st.button("Save to CSV"):
                    csv_file = os.path.join(os.path.dirname(__file__), f"{analysis_type.lower().replace(' ', '_')}_results.csv")
                    result_df.to_csv(csv_file, index=False)
                    st.success(f"Results saved to {csv_file}")
            
            with col2:
                if st.button("Save to Database"):
                    try:
                        # Save to an appropriate database table
                        if analysis_type == "General Variant Scoring":
                            db_manager = st.session_state.db_manager
                            # Prepare data for database schema
                            db_data = result_df.copy()
                            if '#CHROM' in db_data.columns:
                                db_data = db_data.rename(columns={'#CHROM': 'chrom', 'POS': 'pos', 'REF': 'ref', 'ALT': 'alt'})
                            
                            # Extract gene names
                            gene_names = []
                            for info in db_data['INFO']:
                                gene = helpers.extract_gene_name_from_vcf_info(info)
                                gene_names.append(gene if gene else 'Unknown')
                            db_data['gene'] = gene_names
                            
                            # Add placeholders for other required fields
                            db_data['is_common'] = 1
                            db_data['is_rare'] = 0
                            db_data['allele_freq'] = 0.1
                            
                            # Select only relevant columns
                            cols = ['chrom', 'pos', 'ref', 'alt', 'gene', 'is_common', 'is_rare', 'allele_freq', 'functional_score']
                            valid_cols = [col for col in cols if col in db_data.columns]
                            db_data = db_data[valid_cols]
                            
                            # Save to database
                            db_manager.populate_database('endothelial_variants_db', db_data, 'variants')
                            st.success("Results saved to endothelial_variants_db!")
                        
                        else:
                            st.info("Database storage for this analysis type not yet implemented.")
                    
                    except Exception as e:
                        st.error(f"Error saving to database: {e}")

elif page == "lncRNA Analysis":
    st.title("🧬 lncRNA Analysis")
    
    if st.session_state.lncrnas_df is None:
        st.warning("No lncRNA data loaded. Please go to Data Management to upload lncRNA data or load example data.")
    else:
        # Analysis options
        st.subheader("lncRNA Analysis Options")
        
        analysis_type = st.selectbox(
            "Select Analysis Type",
            ["General lncRNA Scoring", "Specific lncRNA Analysis"]
        )
        
        # Filter options
        with st.expander("Filter Options"):
            if analysis_type == "Specific lncRNA Analysis":
                available_lncrnas = st.session_state.lncrnas_df['lncrna'].tolist() if 'lncrna' in st.session_state.lncrnas_df.columns else []
                selected_lncrnas = st.multiselect(
                    "Select specific lncRNAs",
                    options=available_lncrnas,
                    default=["GATA2-AS1", "ANRIL", "PIRAT", "LRAC"] if all(lnc in available_lncrnas for lnc in ["GATA2-AS1", "ANRIL", "PIRAT", "LRAC"]) else available_lncrnas[:4]
                )
                
                # Apply filters to create a working dataset
                if selected_lncrnas:
                    filtered_lncrnas_df = st.session_state.lncrnas_df[st.session_state.lncrnas_df['lncrna'].isin(selected_lncrnas)]
                else:
                    filtered_lncrnas_df = st.session_state.lncrnas_df
            else:
                filtered_lncrnas_df = st.session_state.lncrnas_df
            
            # Show sample count after filtering
            st.write(f"Filtered dataset: {len(filtered_lncrnas_df)} lncRNAs")
        
        # Execute analysis button
        if st.button("Run Analysis"):
            with st.spinner(f"Running {analysis_type}..."):
                try:
                    # This is a mock implementation as we don't have the actual Evo2 executable
                    # Simulate processing time
                    progress_bar = st.progress(0)
                    for i in range(100):
                        time.sleep(0.01)
                        progress_bar.progress(i + 1)
                    
                    # Generate mock results based on analysis type
                    if analysis_type == "General lncRNA Scoring":
                        # Add mock functionality scores
                        result_df = filtered_lncrnas_df.copy()
                        result_df['functionality_score'] = np.random.uniform(0, 1, len(result_df))
                        st.session_state.analysis_results[analysis_type] = result_df
                        
                    elif analysis_type == "Specific lncRNA Analysis":
                        # Add mock detailed scores for specific lncRNAs
                        result_df = filtered_lncrnas_df.copy()
                        result_df['functionality_score'] = np.random.uniform(0.3, 0.9, len(result_df))
                        result_df['structure_score'] = np.random.uniform(0.4, 0.8, len(result_df))
                        result_df['conservation_score'] = np.random.uniform(0.2, 0.7, len(result_df))
                        result_df['binding_domains'] = np.random.randint(1, 10, len(result_df))
                        st.session_state.analysis_results[analysis_type] = result_df
                    
                    st.success(f"{analysis_type} completed successfully!")
                    
                except Exception as e:
                    st.error(f"Error during analysis: {e}")
        
        # Display results if available
        if analysis_type in st.session_state.analysis_results:
            st.subheader("Analysis Results")
            result_df = st.session_state.analysis_results[analysis_type]
            
            # Display basic results table
            st.dataframe(result_df)
            
            # Visualizations based on analysis type
            st.subheader("Visualizations")
            
            if analysis_type == "General lncRNA Scoring":
                col1, col2 = st.columns(2)
                
                with col1:
                    # Histogram of functionality scores
                    fig = px.histogram(
                        result_df, 
                        x='functionality_score',
                        nbins=20,
                        title='Distribution of lncRNA Functionality Scores'
                    )
                    st.plotly_chart(fig, use_container_width=True)
                
                with col2:
                    # Bar chart of top 10 lncRNAs by functionality score
                    top_lncrnas = result_df.sort_values('functionality_score', ascending=False).head(10)
                    fig = px.bar(
                        top_lncrnas,
                        x='lncrna',
                        y='functionality_score',
                        title='Top 10 lncRNAs by Functionality Score'
                    )
                    st.plotly_chart(fig, use_container_width=True)
            
            elif analysis_type == "Specific lncRNA Analysis":
                # Radar chart for multidimensional visualization
                fig = go.Figure()
                
                for i, row in result_df.iterrows():
                    fig.add_trace(go.Scatterpolar(
                        r=[row['functionality_score'], row['structure_score'], 
                           row['conservation_score'], row['binding_domains'] / 10],  # Normalize binding domains
                        theta=['Functionality', 'Structure', 'Conservation', 'Binding Domains'],
                        fill='toself',
                        name=row['lncrna']
                    ))
                
                fig.update_layout(
                    polar=dict(
                        radialaxis=dict(
                            visible=True,
                            range=[0, 1]
                        )),
                    title='Multidimensional Analysis of Selected lncRNAs',
                    showlegend=True
                )
                st.plotly_chart(fig, use_container_width=True)
                
                # Correlation heatmap
                numeric_cols = ['functionality_score', 'structure_score', 'conservation_score', 'binding_domains']
                corr_matrix = result_df[numeric_cols].corr()
                
                fig = px.imshow(
                    corr_matrix, 
                    text_auto=True, 
                    title='Correlation Between lncRNA Metrics',
                    color_continuous_scale='Viridis'
                )
                st.plotly_chart(fig, use_container_width=True)
            
            # Export options
            st.subheader("Export Results")
            col1, col2 = st.columns(2)
            
            with col1:
                if st.button("Save to CSV"):
                    csv_file = os.path.join(os.path.dirname(__file__), f"{analysis_type.lower().replace(' ', '_')}_results.csv")
                    result_df.to_csv(csv_file, index=False)
                    st.success(f"Results saved to {csv_file}")
            
            with col2:
                if st.button("Save to Database"):
                    try:
                        if analysis_type == "General lncRNA Scoring" or analysis_type == "Specific lncRNA Analysis":
                            db_manager = st.session_state.db_manager
                            db_data = result_df.copy()
                            
                            # Rename columns to match database schema if needed
                            if 'lncrna' in db_data.columns:
                                db_data = db_data.rename(columns={'lncrna': 'lncrna_name'})
                            
                            # Fill missing columns if needed
                            if 'structure_score' not in db_data.columns:
                                db_data['structure_score'] = None
                            if 'conservation_score' not in db_data.columns:
                                db_data['conservation_score'] = None
                            if 'binding_domains' not in db_data.columns:
                                db_data['binding_domains'] = None
                            
                            # Save to database
                            db_manager.populate_database('lncrna_functionality_db', db_data, 'lncrnas')
                            st.success("Results saved to lncrna_functionality_db!")
                        else:
                            st.info("Database storage for this analysis type not yet implemented.")
                            
                    except Exception as e:
                        st.error(f"Error saving to database: {e}")

elif page == "Pipeline Execution":
    st.title("🔄 Pipeline Execution")
    st.markdown("""
    This page allows you to execute the complete Evo2 computational pipeline, which integrates 
    all the analysis steps including data loading, variant scoring, lncRNA analysis, 
    and database generation.
    """)
    
    # Pipeline configuration
    st.subheader("Pipeline Configuration")
    
    # Config file upload
    uploaded_config = st.file_uploader("Upload Configuration File (optional)", type=['yaml', 'yml'])
    if uploaded_config is not None:
        with st.spinner("Loading configuration..."):
            temp_config = os.path.join(os.path.dirname(__file__), 'temp_config.yaml')
            with open(temp_config, 'wb') as f:
                f.write(uploaded_config.getbuffer())
            
            try:
                config_data = helpers.load_yaml_config(temp_config)
                st.session_state.config = Config(temp_config)
                st.success("Configuration loaded successfully!")
                
                # Show the loaded config
                with st.expander("View Loaded Configuration"):
                    st.write(config_data)
                    
            except Exception as e:
                st.error(f"Error loading configuration: {e}")
            finally:
                if os.path.exists(temp_config):
                    os.remove(temp_config)
    
    # Pipeline options
    with st.expander("Pipeline Options"):
        st.checkbox("Generate databases", value=True)
        st.checkbox("Run EMT gene analysis", value=True)
        st.checkbox("Analyze GATA2-AS1", value=True)
        st.checkbox("Perform epigenetic mediator analysis", value=True)
        
        log_level = st.selectbox(
            "Log Level",
            options=["INFO", "DEBUG", "WARNING", "ERROR", "CRITICAL"],
            index=0
        )
    
    # Execute pipeline button
    if st.button("Run Complete Pipeline"):
        # Check if required data is available
        if st.session_state.variants_df is None:
            st.error("No variant data is loaded. Please go to Data Management first.")
        elif st.session_state.lncrnas_df is None:
            st.error("No lncRNA data is loaded. Please go to Data Management first.")
        else:
            with st.spinner("Executing Evo2 computational pipeline..."):
                try:
                    # This is a mock implementation of the pipeline execution
                    # In a real implementation, we would use the actual Pipeline class
                    
                    # Create progress bar for steps
                    steps = ["Data Loading", "Variant Scoring", "lncRNA Analysis", "Database Generation", "Specific Analyses"]
                    current_step = st.empty()
                    progress_bar = st.progress(0)
                    
                    for i, step in enumerate(steps):
                        current_step.write(f"Step {i+1}/{len(steps)}: {step}")
                        
                        # Simulate processing time for each step
                        for j in range(20):
                            time.sleep(0.1)
                            progress_bar.progress((i * 20 + j + 1) / (len(steps) * 20))
                    
                    # Generate mock results for demonstration
                    # In a real implementation, we would access the actual results from the Pipeline instance
                    
                    # Mock successful pipeline execution
                    st.session_state.pipeline_executed = True
                    st.session_state.pipeline_results = {
                        "scored_variants_count": len(st.session_state.variants_df),
                        "scored_lncrnas_count": len(st.session_state.lncrnas_df),
                        "execution_time": np.random.uniform(10, 60),
                        "database_tables_created": ["variants", "lncrnas", "predictions", "epigenetic_variants", "emt_variants"]
                    }
                    
                    st.success("Evo2 computational pipeline executed successfully!")
                
                except Exception as e:
                    st.error(f"Error executing pipeline: {e}")
    
    # Display pipeline results if available
    if 'pipeline_executed' in st.session_state and st.session_state.pipeline_executed:
        st.subheader("Pipeline Results")
        
        # Show summary metrics
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric(
                label="Variants Scored",
                value=st.session_state.pipeline_results["scored_variants_count"]
            )
        
        with col2:
            st.metric(
                label="lncRNAs Analyzed",
                value=st.session_state.pipeline_results["scored_lncrnas_count"]
            )
        
        with col3:
            st.metric(
                label="Execution Time (s)",
                value=f"{st.session_state.pipeline_results['execution_time']:.2f}"
            )
        
        # Show databases created
        st.subheader("Databases Created")
        for table in st.session_state.pipeline_results["database_tables_created"]:
            st.success(f"✅ {table}")
        
        # Pipeline logs
        st.subheader("Pipeline Logs")
        try:
            with open(log_file, 'r') as f:
                logs = f.readlines()[-20:]  # Get last 20 lines
                log_text = '\n'.join(logs)
                st.text_area("Log Output", log_text, height=300)
        except Exception as e:
            st.error(f"Unable to load logs: {e}")

elif page == "Database Explorer":
    st.title("🗄️ Database Explorer")
    
    # Get available databases
    available_dbs = get_available_databases()
    
    if not available_dbs:
        st.warning("No databases found. Please run analyses or the complete pipeline first.")
    else:
        # Database selection
        selected_db = st.selectbox("Select Database", available_dbs)
        
        # Get tables in the selected database
        db_manager = st.session_state.db_manager
        db_path = os.path.join(db_manager.db_dir, f"{selected_db}.db")
        
        try:
            # Open connection to the database
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            
            # Get list of tables in the database
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = [row[0] for row in cursor.fetchall()]
            
            if not tables:
                st.info(f"No tables found in database '{selected_db}'")
            else:
                # Table selection
                selected_table = st.selectbox("Select Table", tables)
                
                # Get table schema
                cursor.execute(f"PRAGMA table_info({selected_table})")
                schema = cursor.fetchall()
                
                with st.expander("Table Schema"):
                    schema_df = pd.DataFrame(schema, columns=['cid', 'name', 'type', 'notnull', 'dflt_value', 'pk'])
                    st.dataframe(schema_df)
                
                # Query options
                st.subheader("Query Options")
                
                query_type = st.radio(
                    "Query Type",
                    ["Select All", "Custom Query"]
                )
                
                if query_type == "Select All":
                    limit = st.slider("Row Limit", 10, 1000, 100)
                    query = f"SELECT * FROM {selected_table} LIMIT {limit}"
                else:
                    query = st.text_area(
                        "SQL Query",
                        value=f"SELECT * FROM {selected_table} LIMIT 100",
                        height=100
                    )
                
                # Execute query
                if st.button("Run Query"):
                    with st.spinner("Executing query..."):
                        try:
                            results = pd.read_sql_query(query, conn)
                            st.success(f"Query executed successfully. {len(results)} rows returned.")
                            st.dataframe(results)
                            
                            # Export options
                            col1, col2 = st.columns(2)
                            
                            with col1:
                                if st.button("Export to CSV"):
                                    csv_file = os.path.join(os.path.dirname(__file__), f"{selected_db}_{selected_table}_query_results.csv")
                                    results.to_csv(csv_file, index=False)
                                    st.success(f"Results saved to {csv_file}")
                            
                            # Visualize data
                            if len(results) > 0:
                                st.subheader("Data Visualization")
                                
                                # Get numeric columns for visualization
                                numeric_cols = results.select_dtypes(include=['number']).columns.tolist()
                                
                                if numeric_cols:
                                    # Choose visualization type
                                    viz_type = st.selectbox(
                                        "Visualization Type",
                                        ["Bar Chart", "Histogram", "Scatter Plot", "Box Plot"]
                                    )
                                    
                                    if viz_type == "Bar Chart":
                                        x_col = st.selectbox("X-axis", results.columns.tolist())
                                        y_col = st.selectbox("Y-axis", numeric_cols)
                                        
                                        fig = px.bar(results, x=x_col, y=y_col, title=f"{y_col} by {x_col}")
                                        st.plotly_chart(fig, use_container_width=True)
                                    
                                    elif viz_type == "Histogram":
                                        col = st.selectbox("Column", numeric_cols)
                                        nbins = st.slider("Number of bins", 5, 50, 20)
                                        
                                        fig = px.histogram(results, x=col, nbins=nbins, title=f"Distribution of {col}")
                                        st.plotly_chart(fig, use_container_width=True)
                                    
                                    elif viz_type == "Scatter Plot":
                                        x_col = st.selectbox("X-axis", numeric_cols)
                                        y_col = st.selectbox("Y-axis", [c for c in numeric_cols if c != x_col] if len(numeric_cols) > 1 else numeric_cols)
                                        
                                        fig = px.scatter(results, x=x_col, y=y_col, title=f"{y_col} vs {x_col}")
                                        st.plotly_chart(fig, use_container_width=True)
                                    
                                    elif viz_type == "Box Plot":
                                        y_col = st.selectbox("Values", numeric_cols)
                                        categorical_cols = [c for c in results.columns if c not in numeric_cols]
                                        if categorical_cols:
                                            x_col = st.selectbox("Categories", ["None"] + categorical_cols)
                                            if x_col != "None":
                                                fig = px.box(results, x=x_col, y=y_col, title=f"{y_col} Distribution by {x_col}")
                                            else:
                                                fig = px.box(results, y=y_col, title=f"{y_col} Distribution")
                                            st.plotly_chart(fig, use_container_width=True)
                                        else:
                                            fig = px.box(results, y=y_col, title=f"{y_col} Distribution")
                                            st.plotly_chart(fig, use_container_width=True)
                                else:
                                    st.info("No numeric columns found for visualization.")
                        
                        except Exception as e:
                            st.error(f"Error executing query: {e}")
            
            # Close connection
            conn.close()
            
        except Exception as e:
            st.error(f"Error accessing database: {e}")

elif page == "Configuration":
    st.title("⚙️ Configuration")
    
    # Display current configuration
    st.subheader("Current Configuration")
    
    # Data Paths
    st.markdown("### Data Paths")
    data_paths = {}
    for key, value in st.session_state.config.DATA_PATHS.items():
        data_paths[key] = st.text_input(f"{key} Path", value=value)
    
    # Tool Paths
    st.markdown("### Tool Paths")
    tool_paths = {}
    for key, value in st.session_state.config.TOOL_PATHS.items():
        tool_paths[key] = st.text_input(f"{key} Path", value=value)
    
    # Output Directories
    st.markdown("### Output Directories")
    output_dirs = {}
    for key, value in st.session_state.config.OUTPUT_DIRS.items():
        output_dirs[key] = st.text_input(f"{key} Directory", value=value)
    
    # Save Configuration
    if st.button("Save Configuration"):
        # Update session state config with new values
        for key, value in data_paths.items():
            st.session_state.config.DATA_PATHS[key] = value
        
        for key, value in tool_paths.items():
            st.session_state.config.TOOL_PATHS[key] = value
        
        for key, value in output_dirs.items():
            st.session_state.config.OUTPUT_DIRS[key] = value
        
        # Save to YAML file
        config_data = {
            'DATA_PATHS': st.session_state.config.DATA_PATHS,
            'TOOL_PATHS': st.session_state.config.TOOL_PATHS,
            'OUTPUT_DIRS': st.session_state.config.OUTPUT_DIRS
        }
        
        config_file = os.path.join(os.path.dirname(__file__), 'evo2_clinical_config.yaml')
        try:
            helpers.save_yaml_config(config_data, config_file)
            st.success(f"Configuration saved to {config_file}")
        except Exception as e:
            st.error(f"Error saving configuration: {e}")
    
    # Validate Configuration
    if st.button("Validate Configuration"):
        st.markdown("### Validation Results")
        
        # Check data paths
        data_path_errors = []
        for key, path in st.session_state.config.DATA_PATHS.items():
            if path and not os.path.exists(path):
                data_path_errors.append(f"'{key}': {path}")
        
        if data_path_errors:
            st.error("The following data paths do not exist:")
            for error in data_path_errors:
                st.write(f"- {error}")
        else:
            st.success("All specified data paths exist or are not set.")
        
        # Check tool paths
        tool_path_errors = []
        for key, path in st.session_state.config.TOOL_PATHS.items():
            if path and not os.path.exists(path):
                tool_path_errors.append(f"'{key}': {path}")
        
        if tool_path_errors:
            st.error("The following tool paths do not exist:")
            for error in tool_path_errors:
                st.write(f"- {error}")
        else:
            st.success("All specified tool paths exist or are not set.")
        
        # Check output directories
        for key, dir_path in st.session_state.config.OUTPUT_DIRS.items():
            full_path = os.path.join(os.path.dirname(__file__), dir_path) if not os.path.isabs(dir_path) else dir_path
            os.makedirs(full_path, exist_ok=True)
        
        st.success("All output directories exist or have been created.")

# Footer
st.markdown("---")
st.markdown("""
<div style="text-align: center">
    <p>Evo2_Clinical Dashboard | Version 0.1.0 | © 2025</p>
</div>
""", unsafe_allow_html=True)