"""
Evo2 Clinical Genomic Analysis Platform
Streamlit application that provides a user interface for genomic data analysis,
focusing on CTEPH and lncRNA analysis.
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import io
from Bio import SeqIO
from Bio.Seq import Seq

# Import the classes from the calculations module
from calculations import GWASAnalyzer, RNAAnalyzer, DiseaseAssociationAnalyzer

st.set_page_config(
    page_title="Evo2 Clinical Genomics Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Apply custom styling
st.markdown("""
<style>
    .reportview-container {
        background-color: #f0f2f6;
    }
    .sidebar .sidebar-content {
        background-color: #f0f2f6;
    }
    h1 {
        color: #2c3e50;
    }
    h2 {
        color: #34495e;
    }
</style>
""", unsafe_allow_html=True)

# Title and introduction
st.title("Evo2 Clinical Genomics Platform")
st.markdown("""
This platform provides tools for analyzing genomic and transcriptomic data, with a focus on 
Chronic Thromboembolic Pulmonary Hypertension (CTEPH) and long non-coding RNAs (lncRNAs).
""")

# Sidebar for navigation
st.sidebar.title("Navigation")
analysis_option = st.sidebar.radio(
    "Select Analysis Type",
    ["GWAS Analysis", "RNA Sequence Analysis", "Disease Association Analysis"]
)

# Global variables to store loaded data
gwas_data = None
rna_data = None
phenotype_data = None

# Utility Functions
def load_gwas_data(uploaded_file):
    """Load GWAS data from an uploaded file"""
    if uploaded_file.name.endswith('.tsv') or uploaded_file.name.endswith('.txt'):
        df = pd.read_csv(uploaded_file, sep='\t')
    else:
        df = pd.read_csv(uploaded_file)
    
    # Map column names to expected format if needed
    column_mapping = {
        'SNPS': 'SNP', 
        'CHR_ID': 'CHR', 
        'CHR_POS': 'POS', 
        'P-VALUE': 'P',
        'STRONGEST SNP-RISK ALLELE': 'A1',
        'DISEASE/TRAIT': 'DISEASE'  # Add mapping for disease column
    }
    
    # Rename columns if they exist
    for old_col, new_col in column_mapping.items():
        if old_col in df.columns and new_col not in df.columns:
            df.rename(columns={old_col: new_col}, inplace=True)
    
    # Check if required columns exist
    required_columns = ['SNP', 'CHR', 'POS', 'P', 'A1', 'A2']
    missing_columns = [col for col in required_columns if col not in df.columns]
    
    # If A2 is missing, create a dummy one for demonstration
    if 'A2' in missing_columns and 'A1' in df.columns:
        df['A2'] = df['A1'].apply(lambda x: 'G' if x == 'A' else 'A')
        missing_columns.remove('A2')
    
    if missing_columns:
        st.error(f"Missing required columns: {', '.join(missing_columns)}")
        return None
    
    return df

def load_fasta_data(uploaded_file):
    """Load RNA sequence data from an uploaded FASTA file"""
    # Save the uploaded file to a temporary file
    temp_path = f"temp_{uploaded_file.name}"
    with open(temp_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    
    # Return the path to the temporary file
    return temp_path

def load_phenotype_data(uploaded_file):
    """Load phenotype data from an uploaded file"""
    if uploaded_file.name.endswith('.tsv') or uploaded_file.name.endswith('.txt'):
        df = pd.read_csv(uploaded_file, sep='\t')
    else:
        df = pd.read_csv(uploaded_file)
    
    # Check if required columns exist
    required_columns = ['SUBJECT_ID', 'DISEASE_STATUS']
    missing_columns = [col for col in required_columns if col not in df.columns]
    
    if missing_columns:
        st.error(f"Missing required columns: {', '.join(missing_columns)}")
        return None
    
    return df

# GWAS Analysis Section
if analysis_option == "GWAS Analysis":
    st.header("GWAS Analysis")
    st.markdown("""
    Upload your GWAS data to identify significant genetic variants associated with a disease.
    """)
    
    # File uploader for GWAS data - allow for large files
    uploaded_gwas = st.file_uploader("Upload GWAS Data (CSV or TSV)", type=["csv", "tsv", "txt"])
    
    if uploaded_gwas is not None:
        gwas_data = load_gwas_data(uploaded_gwas)
        
        if gwas_data is not None:
            st.success(f"GWAS data loaded successfully! {gwas_data.shape[0]} variants found.")
            
            # Display sample data
            st.subheader("Sample Data")
            st.dataframe(gwas_data.head())
            
            # Combined filtering section
            st.subheader("Filter Variants")
            st.markdown("""
            Apply one or more filters to the GWAS data. All selected filters will be applied when you click the "Apply Filters" button.
            """)
            
            with st.expander("Available Filters", expanded=True):
                # 1. Disease filter
                if 'DISEASE' in gwas_data.columns:
                    st.markdown("### Filter by Disease")
                    
                    # Get unique disease values for reference
                    unique_diseases = sorted(gwas_data['DISEASE'].dropna().unique())
                    
                    # Create text input for disease search
                    disease_query = st.text_input("Enter disease name (will filter by partial matches):", "")
                    
                    # Show examples of available diseases
                    if unique_diseases and len(unique_diseases) > 0:
                        sample_size = min(5, len(unique_diseases))
                        st.info(f"Examples of diseases in the dataset: {', '.join(unique_diseases[:sample_size])}")
                    
                    # Add checkbox to enable the filter
                    use_disease_filter = st.checkbox("Apply disease filter", value=False)
                    
                    # If checkbox is checked but no query is entered, show a warning
                    if use_disease_filter and not disease_query:
                        st.warning("Please enter a disease term to filter by, or uncheck the filter.")
                else:
                    use_disease_filter = False
                    disease_query = ""
                
                # 2. P-value filter
                st.markdown("### Filter by P-value")
                p_threshold = st.number_input(
                    "P-value threshold (variants with p-values below this threshold will be selected)",
                    min_value=0.0,
                    max_value=1.0,
                    value=0.0005,
                    format="%.8f",
                    step=0.00001
                )
                use_pvalue_filter = st.checkbox("Apply p-value filter", value=False)
                
                # 3. Genomic region filter
                st.markdown("### Filter by Genomic Region")
                col1, col2, col3 = st.columns(3)
                with col1:
                    chromosome = st.text_input("Chromosome", "1")
                with col2:
                    start_pos = st.number_input("Start Position", value=1000000)
                with col3:
                    end_pos = st.number_input("End Position", value=2000000)
                use_region_filter = st.checkbox("Apply genomic region filter", value=False)
            
            # Apply all selected filters
            if st.button("Apply Filters"):
                try:
                    # Start with full dataset
                    filtered_data = gwas_data.copy()
                    filters_applied = []
                    
                    # Apply disease filter if selected
                    if use_disease_filter and disease_query:
                        filtered_data = filtered_data[filtered_data['DISEASE'].str.contains(disease_query, case=False, na=False)]
                        filters_applied.append(f"Disease contains: '{disease_query}'")
                    
                    # Apply p-value filter if selected
                    if use_pvalue_filter:
                        filtered_data = filtered_data[filtered_data['P'] < p_threshold]
                        filters_applied.append(f"P-value < {p_threshold}")
                    
                    # Apply genomic region filter if selected
                    if use_region_filter:
                        # Make a local copy to avoid modifying the filtered data
                        region_data = filtered_data.copy()
                        
                        # Ensure POS column is numeric
                        if not pd.api.types.is_numeric_dtype(region_data['POS']):
                            region_data['POS'] = pd.to_numeric(region_data['POS'], errors='coerce')
                            # Drop rows where conversion failed
                            region_data = region_data.dropna(subset=['POS'])
                        
                        # Normalize chromosome format
                        sample_chrs = region_data['CHR'].astype(str).unique()[:5]
                        has_chr_prefix = any(str(c).startswith('chr') for c in sample_chrs)
                        
                        search_chromosome = chromosome
                        if has_chr_prefix and not str(search_chromosome).startswith('chr'):
                            search_chromosome = f"chr{search_chromosome}"
                        elif not has_chr_prefix and str(search_chromosome).startswith('chr'):
                            search_chromosome = str(search_chromosome)[3:]
                        
                        # Filter by region
                        filtered_data = region_data[
                            (region_data['CHR'].astype(str) == str(search_chromosome)) & 
                            (region_data['POS'].astype(float) >= float(start_pos)) & 
                            (region_data['POS'].astype(float) <= float(end_pos))
                        ]
                        
                        filters_applied.append(f"Region: Chromosome {search_chromosome}, positions {start_pos}-{end_pos}")
                    
                    # Display results
                    if not filtered_data.empty:
                        if filters_applied:
                            st.success(f"Found {filtered_data.shape[0]} variants matching filters: {', '.join(filters_applied)}")
                        else:
                            st.success(f"No filters applied. Showing all {filtered_data.shape[0]} variants.")
                        
                        st.dataframe(filtered_data)
                        
                        # Visualization if we have filtered data
                        if 'P' in filtered_data.columns and len(filtered_data) > 0:
                            st.subheader("Visualization")
                            
                            # Manhattan plot-like visualization (modified to be short and wide)
                            fig, ax = plt.subplots(figsize=(20, 5))
                            
                            # Convert chromosome to numeric for plotting
                            def chr_to_num(chrom):
                                if str(chrom).startswith('chr'):
                                    chrom = str(chrom)[3:]
                                try:
                                    return int(chrom)
                                except ValueError:
                                    # Handle X, Y, etc.
                                    if chrom == 'X':
                                        return 23
                                    elif chrom == 'Y':
                                        return 24
                                    elif chrom == 'MT':
                                        return 25
                                    else:
                                        return 0
                            
                            # Add CHR_NUM column for plotting if not already present
                            if 'CHR_NUM' not in filtered_data.columns:
                                filtered_data['CHR_NUM'] = filtered_data['CHR'].apply(chr_to_num)
                            
                            # Check if a specific chromosome is selected through the region filter
                            if use_region_filter:
                                # Chromosome-specific plot: position vs. -log10(p-value)
                                scatter = ax.scatter(
                                    filtered_data['POS'], 
                                    -np.log10(filtered_data['P']),
                                    alpha=0.7, 
                                    s=15,
                                    c=-np.log10(filtered_data['P']),  # Color by significance
                                    cmap='viridis'
                                )
                                
                                # Increase tick frequency for high detail
                                ax.xaxis.set_major_locator(plt.MaxNLocator(20))
                                
                                ax.set_xlabel(f'Chromosome {chromosome} Position')
                                ax.set_title(f'Manhattan Plot for Chromosome {chromosome}')
                            else:
                                # Standard Manhattan plot with all chromosomes
                                scatter = ax.scatter(
                                    filtered_data['CHR_NUM'], 
                                    -np.log10(filtered_data['P']),
                                    alpha=0.7, 
                                    s=15,
                                    c=-np.log10(filtered_data['P']),  # Color by significance
                                    cmap='viridis'
                                )
                                
                                # Increase tick frequency for high detail
                                ax.xaxis.set_major_locator(plt.MaxNLocator(24))
                                
                                ax.set_xlabel('Chromosome')
                                ax.set_title('Manhattan Plot of Filtered Variants')
                            
                            # Add threshold line if p-value filter is applied
                            if use_pvalue_filter:
                                ax.axhline(-np.log10(p_threshold), color='red', linestyle='--', 
                                           label=f'p-value threshold: {p_threshold}')
                                ax.legend()
                            
                            ax.set_ylabel('-log10(p-value)')
                            ax.yaxis.set_major_locator(plt.MaxNLocator(10))
                            
                            # Add colorbar
                            cbar = plt.colorbar(scatter)
                            cbar.set_label('-log10(p-value)')
                            
                            # Adjust layout to make sure everything is visible
                            plt.tight_layout()
                            
                            st.pyplot(fig)
                        
                        # Option to download results
                        csv = filtered_data.to_csv(index=False)
                        
                        # Create a descriptive filename based on applied filters
                        if filters_applied:
                            filename = "_".join([f.split(":")[0].lower().replace(" ", "_") for f in filters_applied])
                            filename = f"filtered_variants_{filename}.csv"
                        else:
                            filename = "all_variants.csv"
                        
                        st.download_button(
                            label="Download Filtered Variants as CSV",
                            data=csv,
                            file_name=filename,
                            mime="text/csv"
                        )
                        
                        # Update the global gwas_data for subsequent analyses
                        st.info("Filtered data is now available for further analysis.")
                        gwas_data = filtered_data
                    else:
                        st.warning("No variants match the selected filters. Try adjusting your filter criteria.")
                        
                        # If genomic region filter was applied, show sample chromosomes
                        if use_region_filter:
                            sample_chrs = gwas_data['CHR'].astype(str).unique()[:5]
                            st.info(f"Your data contains chromosomes like: {', '.join(str(c) for c in sample_chrs)}")
                
                except Exception as e:
                    st.error(f"Error during filtering: {str(e)}")
                    import traceback
                    st.error(f"Detailed error: {traceback.format_exc()}")
            
            # Additional analysis options
            st.subheader("Additional Analysis")
            
            # Calculate allele frequencies
            if st.button("Calculate Allele Frequencies"):
                try:
                    gwas_analyzer = GWASAnalyzer(gwas_data)
                    allele_freq_data = gwas_analyzer.calculate_allele_frequencies()
                    
                    st.success("Allele frequencies calculated successfully!")
                    st.dataframe(allele_freq_data)
                    
                    # Visualization of MAF distribution
                    fig, ax = plt.subplots(figsize=(10, 6))
                    sns.histplot(allele_freq_data['MAF'], bins=20, kde=True, ax=ax)
                    ax.set_xlabel('Minor Allele Frequency (MAF)')
                    ax.set_ylabel('Count')
                    ax.set_title('Distribution of Minor Allele Frequencies')
                    st.pyplot(fig)
                    
                    # Option to download results
                    csv = allele_freq_data.to_csv(index=False)
                    st.download_button(
                        label="Download Allele Frequencies as CSV",
                        data=csv,
                        file_name="allele_frequencies.csv",
                        mime="text/csv"
                    )
                except Exception as e:
                    st.error(f"Error during analysis: {str(e)}")

# RNA Sequence Analysis Section
elif analysis_option == "RNA Sequence Analysis":
    st.header("RNA Sequence Analysis")
    st.markdown("""
    Upload RNA sequence data in FASTA format to analyze sequence properties and predict structures.
    """)
    
    # File uploader for RNA sequence data
    uploaded_fasta = st.file_uploader("Upload RNA Sequences (FASTA)", type=["fasta", "fa", "fna"])
    
    if uploaded_fasta is not None:
        fasta_path = load_fasta_data(uploaded_fasta)
        
        try:
            rna_analyzer = RNAAnalyzer(fasta_path)
            sequences = rna_analyzer.sequences
            
            if sequences:
                st.success(f"RNA sequence data loaded successfully! {len(sequences)} sequences found.")
                
                # Display sample data
                st.subheader("Sample Sequences")
                sample_data = []
                for i, record in enumerate(sequences[:5]):  # Show first 5 sequences
                    seq_str = str(record.seq)
                    sample_data.append({
                        "ID": record.id,
                        "Length": len(seq_str),
                        "Sequence Preview": seq_str[:50] + "..." if len(seq_str) > 50 else seq_str
                    })
                
                st.table(sample_data)
                
                # Options for analysis
                st.subheader("Analysis Options")
                
                # 1. Calculate sequence statistics
                st.markdown("### Sequence Statistics")
                
                if st.button("Calculate Sequence Statistics"):
                    try:
                        sequence_stats = rna_analyzer.calculate_sequence_statistics()
                        
                        st.success("Sequence statistics calculated successfully!")
                        st.dataframe(sequence_stats)
                        
                        # Visualizations
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            # Length distribution
                            fig1, ax1 = plt.subplots(figsize=(8, 5))
                            sns.histplot(sequence_stats['length'], bins=20, kde=True, ax=ax1)
                            ax1.set_xlabel('Sequence Length')
                            ax1.set_ylabel('Count')
                            ax1.set_title('Distribution of Sequence Lengths')
                            st.pyplot(fig1)
                        
                        with col2:
                            # GC content distribution
                            fig2, ax2 = plt.subplots(figsize=(8, 5))
                            sns.histplot(sequence_stats['GC_content'], bins=20, kde=True, ax=ax2)
                            ax2.set_xlabel('GC Content (%)')
                            ax2.set_ylabel('Count')
                            ax2.set_title('Distribution of GC Content')
                            st.pyplot(fig2)
                        
                        # Option to download results
                        csv = sequence_stats.to_csv(index=False)
                        st.download_button(
                            label="Download Sequence Statistics as CSV",
                            data=csv,
                            file_name="sequence_statistics.csv",
                            mime="text/csv"
                        )
                    except Exception as e:
                        st.error(f"Error during analysis: {str(e)}")
                
                # 2. Predict RNA structure (placeholder)
                st.markdown("### RNA Structure Prediction")
                
                selected_sequence = st.selectbox(
                    "Select a sequence for structure prediction",
                    options=[record.id for record in sequences]
                )
                
                if st.button("Predict Structure"):
                    st.info("RNA structure prediction is a placeholder in the current implementation.")
                    st.warning("In a real implementation, this would integrate with external tools like RNAfold or ViennaRNA.")
                    
                    # Placeholder visualization
                    selected_record = next((r for r in sequences if r.id == selected_sequence), None)
                    if selected_record:
                        seq_str = str(selected_record.seq)
                        
                        # Simple visualization of nucleotide composition
                        nucleotides = ['A', 'C', 'G', 'U', 'T', 'N']
                        counts = [seq_str.count(n) for n in nucleotides]
                        
                        fig, ax = plt.subplots(figsize=(8, 6))
                        ax.bar(nucleotides, counts)
                        ax.set_xlabel('Nucleotide')
                        ax.set_ylabel('Count')
                        ax.set_title(f'Nucleotide Composition of {selected_sequence}')
                        st.pyplot(fig)
                
                # 3. Predict lncRNA essentiality (placeholder)
                st.markdown("### lncRNA Essentiality Prediction")
                
                if st.button("Predict Essentiality"):
                    st.info("lncRNA essentiality prediction is a placeholder in the current implementation.")
                    st.warning("In a real implementation, this would use a sophisticated model based on sequence features, conservation, and more.")
                    
                    # Generate placeholder results
                    essentiality_results = []
                    for record in sequences[:10]:  # Limit to first 10 for display
                        essentiality_results.append({
                            "ID": record.id,
                            "Essentiality Score": rna_analyzer.predict_lncrna_essentiality(record.id),
                            "Confidence": "Low (Placeholder)"
                        })
                    
                    st.table(essentiality_results)
            else:
                st.error("No sequences found in the uploaded file.")
        except Exception as e:
            st.error(f"Error during RNA analysis: {str(e)}")
        
        # Clean up temporary file
        if os.path.exists(fasta_path):
            os.remove(fasta_path)

# Disease Association Analysis Section
elif analysis_option == "Disease Association Analysis":
    st.header("Disease Association Analysis")
    st.markdown("""
    Analyze associations between genetic variants, lncRNA features, and disease phenotypes.
    This analysis requires both GWAS data and phenotype data.
    """)
    
    # File uploaders
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("GWAS Data")
        uploaded_gwas = st.file_uploader("Upload GWAS Data (CSV or TSV)", 
                                         type=["csv", "tsv", "txt"],
                                         key="gwas_disease")
    
    with col2:
        st.subheader("Phenotype Data")
        uploaded_phenotype = st.file_uploader("Upload Phenotype Data (CSV or TSV)",
                                             type=["csv", "tsv", "txt"])
    
    # Load data when both files are uploaded
    if uploaded_gwas is not None and uploaded_phenotype is not None:
        gwas_data = load_gwas_data(uploaded_gwas)
        phenotype_data = load_phenotype_data(uploaded_phenotype)
        
        if gwas_data is not None and phenotype_data is not None:
            st.success("Both GWAS and phenotype data loaded successfully!")
            
            # Display sample data
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("GWAS Data Sample")
                st.dataframe(gwas_data.head())
            
            with col2:
                st.subheader("Phenotype Data Sample")
                st.dataframe(phenotype_data.head())
            
            # Options for analysis
            st.subheader("Analysis Options")
            
            # 1. Assess variant association
            st.markdown("### Variant-Disease Association Analysis")
            
            p_threshold = st.slider("P-value Threshold for Variant Selection", 
                                   0.0, 0.1, 0.05, 0.001,
                                   key="p_threshold_disease")
            
            if st.button("Perform Association Analysis"):
                try:
                    # Check if SUBJECT_ID exists in gwas_data
                    if 'SUBJECT_ID' not in gwas_data.columns:
                        st.error("GWAS data must contain a 'SUBJECT_ID' column for this analysis.")
                    else:
                        gwas_analyzer = GWASAnalyzer(gwas_data)
                        disease_analyzer = DiseaseAssociationAnalyzer(gwas_analyzer, None, phenotype_data)
                        
                        association_results = disease_analyzer.assess_variant_association(p_threshold)
                        
                        if not association_results.empty:
                            st.success("Association analysis completed successfully!")
                            st.dataframe(association_results)
                            
                            # Visualization
                            fig, ax = plt.subplots(figsize=(10, 6))
                            
                            # Create forest plot for odds ratios
                            y_pos = np.arange(len(association_results))
                            ax.errorbar(
                                x=association_results['odds_ratio'],
                                y=y_pos,
                                xerr=[
                                    association_results['odds_ratio'] - association_results['conf_lower'],
                                    association_results['conf_upper'] - association_results['odds_ratio']
                                ],
                                fmt='o',
                                capsize=5
                            )
                            
                            ax.axvline(x=1.0, color='red', linestyle='--', alpha=0.7)
                            ax.set_yticks(y_pos)
                            ax.set_yticklabels(association_results['SNP'])
                            ax.set_xlabel('Odds Ratio (95% CI)')
                            ax.set_title('Forest Plot of Variant-Disease Associations')
                            plt.tight_layout()
                            
                            st.pyplot(fig)
                            
                            # Option to download results
                            csv = association_results.to_csv(index=False)
                            st.download_button(
                                label="Download Association Results as CSV",
                                data=csv,
                                file_name="variant_disease_associations.csv",
                                mime="text/csv"
                            )
                        else:
                            st.warning("No significant associations found with the given threshold.")
                except Exception as e:
                    st.error(f"Error during association analysis: {str(e)}")
            
            # 2. Compare lncRNA expression
            st.markdown("### lncRNA Expression Comparison")
            
            # File uploader for expression data
            uploaded_expression = st.file_uploader("Upload Expression Data (CSV or TSV)",
                                                 type=["csv", "tsv", "txt"])
            
            if uploaded_expression is not None:
                try:
                    if uploaded_expression.name.endswith('.tsv') or uploaded_expression.name.endswith('.txt'):
                        expression_data = pd.read_csv(uploaded_expression, sep='\t')
                    else:
                        expression_data = pd.read_csv(uploaded_expression)
                    
                    st.success("Expression data loaded successfully!")
                    st.dataframe(expression_data.head())
                    
                    # Select lncRNA to compare
                    lncrna_options = [col for col in expression_data.columns 
                                     if col not in ['SUBJECT_ID', 'DISEASE_STATUS']]
                    
                    selected_lncrna = st.selectbox("Select lncRNA to Compare", options=lncrna_options)
                    
                    if st.button("Compare Expression"):
                        try:
                            # Create a temporary RNAAnalyzer (not used but required for constructor)
                            temp_fasta = os.path.join(os.getcwd(), "data", "reference", "GRCh38.chr19.fa")
                            if not os.path.exists(temp_fasta):
                                st.warning("Reference FASTA file not found. Using a placeholder.")
                                temp_fasta = "placeholder.fa"
                                with open(temp_fasta, "w") as f:
                                    f.write(">placeholder\nACGT\n")
                            
                            rna_analyzer = RNAAnalyzer(temp_fasta)
                            gwas_analyzer = GWASAnalyzer(gwas_data)
                            disease_analyzer = DiseaseAssociationAnalyzer(gwas_analyzer, rna_analyzer, phenotype_data)
                            
                            comparison_results = disease_analyzer.compare_lncrna_expression(
                                expression_data, selected_lncrna)
                            
                            st.success("Expression comparison completed successfully!")
                            st.dataframe(comparison_results)
                            
                            # Visualization - Box plot of expression by group
                            fig, ax = plt.subplots(figsize=(10, 6))
                            
                            sns.boxplot(
                                x='DISEASE_STATUS',
                                y=selected_lncrna,
                                data=expression_data,
                                ax=ax
                            )
                            
                            ax.set_xlabel('Disease Status (0=Control, 1=Case)')
                            ax.set_ylabel(f'{selected_lncrna} Expression')
                            ax.set_title(f'{selected_lncrna} Expression by Disease Status')
                            
                            # Add p-value to the plot
                            p_value = comparison_results['p_value'].values[0]
                            ax.text(
                                0.5, 0.95,
                                f'p-value: {p_value:.4f}',
                                horizontalalignment='center',
                                verticalalignment='center',
                                transform=ax.transAxes,
                                bbox=dict(facecolor='white', alpha=0.8)
                            )
                            
                            st.pyplot(fig)
                            
                            # Clean up temporary file if created
                            if os.path.exists("placeholder.fa"):
                                os.remove("placeholder.fa")
                        except Exception as e:
                            st.error(f"Error during expression comparison: {str(e)}")
                except Exception as e:
                    st.error(f"Error loading expression data: {str(e)}")
            
            # 3. Bayesian prioritization (placeholder)
            st.markdown("### Bayesian Prioritization of Variant-lncRNA Pairs")
            
            if st.button("Perform Bayesian Prioritization"):
                st.info("Bayesian prioritization is a placeholder in the current implementation.")
                st.warning("In a real implementation, this would use a more sophisticated model with real prior data.")
                
                # Generate placeholder data
                placeholder_pairs = [
                    {"variant_id": "rs123", "lncrna_id": "ANRIL"},
                    {"variant_id": "rs456", "lncrna_id": "HOTAIR"},
                    {"variant_id": "rs789", "lncrna_id": "MALAT1"}
                ]
                
                placeholder_priors = {
                    "rs123-ANRIL": 0.2,
                    "rs456-HOTAIR": 0.15,
                    "rs789-MALAT1": 0.1
                }
                
                # Create disease analyzer (without RNA analyzer)
                gwas_analyzer = GWASAnalyzer(gwas_data)
                disease_analyzer = DiseaseAssociationAnalyzer(gwas_analyzer, None, phenotype_data)
                
                # Get placeholder results
                prioritization_results = disease_analyzer.bayesian_prioritization(
                    placeholder_pairs, placeholder_priors)
                
                st.table(prioritization_results)

# Footer
st.markdown("---")
st.markdown("© 2025 Evo2 Clinical Genomics Platform")
st.markdown("""
This is a demo application. For research purposes only.
Developed based on the genomic analysis tools in the Evo2_Clinical package.
""")