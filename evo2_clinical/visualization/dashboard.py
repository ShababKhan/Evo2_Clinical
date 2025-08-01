"""
Dashboard Module for Evo2_Clinical.

This module provides the main dashboard functionality that integrates
all visualization components into a cohesive Streamlit-based interface.
"""

import os
import json
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Union, Tuple
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from evo2_clinical.visualization.variant_viz import VariantVisualizer
from evo2_clinical.visualization.pathway_viz import PathwayVisualizer
from evo2_clinical.visualization.treatment_viz import TreatmentOutcomeVisualizer
from evo2_clinical.config import config


class Dashboard:
    """Main dashboard class for Evo2_Clinical visualizations."""
    
    def __init__(self, theme: str = "light"):
        """
        Initialize the dashboard.
        
        Args:
            theme (str): Color theme for visualizations ("light" or "dark")
        """
        self.theme = theme
        self.variant_viz = VariantVisualizer(theme=theme)
        self.pathway_viz = PathwayVisualizer(theme=theme)
        self.treatment_viz = TreatmentOutcomeVisualizer(theme=theme)
        
    def load_variant_data(self, file_path: str) -> pd.DataFrame:
        """
        Load variant data from a file.
        
        Args:
            file_path (str): Path to the variant data file (CSV or JSON)
            
        Returns:
            pd.DataFrame: DataFrame with variant data
        """
        file_ext = os.path.splitext(file_path)[1].lower()
        
        try:
            if file_ext == '.csv':
                return pd.read_csv(file_path)
            elif file_ext == '.json':
                return pd.read_json(file_path)
            else:
                raise ValueError(f"Unsupported file type: {file_ext}")
        except Exception as e:
            st.error(f"Error loading variant data: {e}")
            return pd.DataFrame()
    
    def load_pathway_data(self, file_path: str) -> Tuple[Dict[str, float], Optional[Dict[str, float]]]:
        """
        Load pathway data from a file.
        
        Args:
            file_path (str): Path to the pathway data file (JSON)
            
        Returns:
            Tuple[Dict[str, float], Optional[Dict[str, float]]]: 
                Pathway scores and gene scores (if available)
        """
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            
            # Check if we have the expected format
            if 'pathway_details' in data:
                pathway_scores = {
                    pathway: details.get('impact_score', 0) 
                    for pathway, details in data['pathway_details'].items()
                }
                
                # Extract gene scores if available
                gene_scores = None
                if 'variant_contributions' in data:
                    gene_scores = {}
                    for variant in data['variant_contributions']:
                        if 'gene' in variant and 'impact_score' in variant:
                            gene = variant['gene']
                            score = variant['impact_score']
                            gene_scores[gene] = score
                
                return pathway_scores, gene_scores
            else:
                # Simple format where the file is just a dictionary of pathway -> score
                return data, None
                
        except Exception as e:
            st.error(f"Error loading pathway data: {e}")
            return {}, None
    
    def load_treatment_data(self, file_path: str) -> Dict:
        """
        Load treatment outcome data from a file.
        
        Args:
            file_path (str): Path to the treatment data file (JSON)
            
        Returns:
            Dict: Dictionary with treatment outcome data
        """
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            return data
        except Exception as e:
            st.error(f"Error loading treatment data: {e}")
            return {}
    
    def render_variant_section(
        self, 
        variants_df: pd.DataFrame,
        section_title: str = "Genetic Variant Analysis"
    ):
        """
        Render the variant analysis section of the dashboard.
        
        Args:
            variants_df (pd.DataFrame): DataFrame with variant information
            section_title (str): Title for the section
        """
        st.header(section_title)
        
        if variants_df.empty:
            st.warning("No variant data available.")
            return
        
        # Create tabs for different visualizations
        variant_tabs = st.tabs([
            "Variant Summary", 
            "Distribution", 
            "Impact Scores", 
            "Lollipop Plot",
            "Comparison"
        ])
        
        with variant_tabs[0]:
            st.subheader("Variant Summary Table")
            
            # Get default columns to display
            summary_df = self.variant_viz._get_default_variant_columns(variants_df)
            st.dataframe(summary_df, use_container_width=True)
            
            # Show download button for the data
            csv_data = summary_df.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="Download Data",
                data=csv_data,
                file_name="variant_summary.csv",
                mime="text/csv",
                key="download-variant-summary"
            )
        
        with variant_tabs[1]:
            st.subheader("Variant Distribution")
            
            # Select grouping column
            col_options = [col for col in variants_df.columns 
                         if variants_df[col].nunique() < len(variants_df) / 2]
            if not col_options:
                col_options = ["chrom"] if "chrom" in variants_df.columns else variants_df.columns[:1].tolist()
                
            group_col = st.selectbox(
                "Group variants by:",
                options=col_options,
                index=0
            )
            
            # Create distribution plot
            fig = self.variant_viz.plot_variant_distribution(
                variants_df=variants_df,
                groupby_column=group_col,
                title=f"Variant Distribution by {group_col.title()}",
                interactive=True
            )
            
            if fig is not None:
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.warning(f"Unable to create distribution plot using column: {group_col}")
        
        with variant_tabs[2]:
            st.subheader("Variant Impact Scores")
            
            # Find score columns
            score_columns = [col for col in variants_df.columns 
                           if 'score' in col.lower() or 'impact' in col.lower()]
            if not score_columns:
                st.warning("No score columns found in the data.")
            else:
                # Select score column
                score_col = st.selectbox(
                    "Select score column:",
                    options=score_columns,
                    index=0
                )
                
                # Optional threshold
                use_threshold = st.checkbox("Use threshold")
                threshold = None
                
                if use_threshold:
                    min_val = float(variants_df[score_col].min())
                    max_val = float(variants_df[score_col].max())
                    threshold = st.slider(
                        "Threshold value:",
                        min_value=min_val,
                        max_value=max_val,
                        value=(min_val + max_val) / 2
                    )
                
                # Create score plot
                fig = self.variant_viz.plot_variant_scores(
                    variants_df=variants_df,
                    score_column=score_col,
                    threshold=threshold,
                    title=f"{score_col.replace('_', ' ').title()} Distribution",
                    interactive=True
                )
                
                if fig is not None:
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.warning(f"Unable to create score plot using column: {score_col}")
        
        with variant_tabs[3]:
            st.subheader("Variant Position Lollipop Plot")
            
            # Check if we have the necessary columns
            required_cols = ["pos"]
            score_columns = [col for col in variants_df.columns 
                           if 'score' in col.lower() or 'impact' in col.lower()]
            
            if "pos" not in variants_df.columns:
                st.warning("Position column 'pos' not found in data.")
            elif not score_columns:
                st.warning("No score columns found in the data.")
            else:
                # Select gene filter if available
                gene_col = None
                gene = None
                
                for col in ["gene", "gene_name", "gene_id"]:
                    if col in variants_df.columns:
                        gene_col = col
                        break
                
                if gene_col:
                    genes = sorted(variants_df[gene_col].unique().tolist())
                    gene = st.selectbox(
                        "Filter by gene:",
                        options=["All genes"] + genes,
                        index=0
                    )
                    
                    if gene == "All genes":
                        gene = None
                
                # Select score column
                score_col = st.selectbox(
                    "Select score column for lollipop plot:",
                    options=score_columns,
                    index=0 if score_columns else 0
                )
                
                # Create lollipop plot
                fig = self.variant_viz.plot_variant_lollipop(
                    variants_df=variants_df,
                    gene_column=gene_col or "gene",
                    score_column=score_col,
                    position_column="pos",
                    gene=gene,
                    interactive=True
                )
                
                if fig is not None:
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.warning("Unable to create lollipop plot with the current selections.")
        
        with variant_tabs[4]:
            st.subheader("Variant Effect Comparison")
            
            # Find score columns for comparison
            score_columns = [col for col in variants_df.columns 
                           if 'score' in col.lower() or 'impact' in col.lower() 
                           or 'effect' in col.lower()]
            
            if len(score_columns) < 2:
                st.warning("Need at least 2 score columns for comparison.")
            else:
                # Select score columns to compare
                selected_scores = st.multiselect(
                    "Select scores to compare:",
                    options=score_columns,
                    default=score_columns[:min(5, len(score_columns))]
                )
                
                if len(selected_scores) < 2:
                    st.warning("Please select at least 2 score columns.")
                else:
                    # Determine variant ID column
                    id_col = None
                    for col in ["id", "variant_id", "var_id"]:
                        if col in variants_df.columns:
                            id_col = col
                            break
                    
                    # Limit to a manageable number of variants
                    max_variants = st.slider(
                        "Maximum number of variants to display:",
                        min_value=3,
                        max_value=20,
                        value=8
                    )
                    
                    # Create subset for visualization
                    if len(variants_df) > max_variants:
                        # Try to select variants with highest overall scores
                        subset_df = variants_df.iloc[:max_variants].copy()
                    else:
                        subset_df = variants_df.copy()
                    
                    # Create comparison plot
                    fig = self.variant_viz.plot_variant_effect_comparison(
                        variants_df=subset_df,
                        effect_columns=selected_scores,
                        variant_id_column=id_col or "id",
                        title="Comparison of Variant Effect Scores",
                        interactive=True
                    )
                    
                    if fig is not None:
                        st.plotly_chart(fig, use_container_width=True)
                    else:
                        st.warning("Unable to create comparison plot with the current selections.")
    
    def render_pathway_section(
        self,
        pathway_scores: Dict[str, float],
        gene_scores: Optional[Dict[str, float]] = None,
        section_title: str = "Pathway Analysis"
    ):
        """
        Render the pathway analysis section of the dashboard.
        
        Args:
            pathway_scores (Dict[str, float]): Dictionary of pathway names and scores
            gene_scores (Dict[str, float], optional): Dictionary of gene names and scores
            section_title (str): Title for the section
        """
        st.header(section_title)
        
        if not pathway_scores:
            st.warning("No pathway data available.")
            return
        
        # Create tabs for different visualizations
        pathway_tabs = st.tabs([
            "Pathway Scores", 
            "Gene Contributions",
            "Network Visualization"
        ])
        
        with pathway_tabs[0]:
            st.subheader("Pathway Impact Scores")
            
            # Create pathway score visualization
            fig = self.pathway_viz.plot_emt_pathway(
                pathway_scores=pathway_scores,
                gene_scores=None,  # Show gene scores in another tab
                title="Pathway Activation/Inhibition Scores",
                interactive=True
            )
            
            if fig is not None:
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.warning("Unable to create pathway score visualization.")
            
            # Show data table
            st.subheader("Pathway Data")
            pathway_df = pd.DataFrame({
                'Pathway': list(pathway_scores.keys()),
                'Impact Score': list(pathway_scores.values())
            }).sort_values('Impact Score', key=abs, ascending=False)
            
            st.dataframe(pathway_df, use_container_width=True)
        
        with pathway_tabs[1]:
            st.subheader("Gene Contributions")
            
            if not gene_scores:
                st.info("No gene contribution data available.")
            else:
                # Create gene contribution visualization
                gene_items = sorted(gene_scores.items(), key=lambda x: abs(x[1]), reverse=True)
                gene_names = [item[0] for item in gene_items]
                gene_values = [item[1] for item in gene_items]
                
                # Limit to top genes
                max_genes = st.slider(
                    "Maximum number of genes to display:",
                    min_value=5,
                    max_value=50,
                    value=20
                )
                
                if len(gene_names) > max_genes:
                    gene_names = gene_names[:max_genes]
                    gene_values = gene_values[:max_genes]
                
                # Create bar chart
                gene_fig = px.bar(
                    x=gene_names,
                    y=gene_values,
                    labels={'x': 'Gene', 'y': 'Contribution Score'},
                    title="Top Contributing Genes",
                    color=gene_values,
                    color_continuous_scale="RdBu_r"
                )
                
                gene_fig.update_layout(
                    xaxis_tickangle=45,
                    template="plotly_white" if self.theme == "light" else "plotly_dark"
                )
                
                st.plotly_chart(gene_fig, use_container_width=True)
                
                # Show data table
                st.subheader("Gene Data")
                gene_df = pd.DataFrame({
                    'Gene': gene_names,
                    'Contribution Score': gene_values
                })
                
                st.dataframe(gene_df, use_container_width=True)
        
        with pathway_tabs[2]:
            st.subheader("Network Visualization")
            st.info("To visualize pathway networks, you need to provide node and edge data files.")
            
            # Example of how to use the network visualization if data is available
            st.markdown("""
            ### Sample Network Visualization Code:
            ```python
            # Load node and edge data
            nodes_df = pd.read_csv('path/to/nodes.csv')
            edges_df = pd.read_csv('path/to/edges.csv')
            
            # Create network visualization
            network_fig = pathway_viz.plot_pathway_network(
                nodes=nodes_df,
                edges=edges_df,
                node_id_col="id",
                source_col="source",
                target_col="target",
                node_type_col="type",
                title="EMT Pathway Network"
            )
            
            # Display the network
            st.plotly_chart(network_fig, use_container_width=True)
            ```
            """)
    
    def render_treatment_section(
        self,
        treatment_data: Dict,
        section_title: str = "Treatment Outcome Analysis"
    ):
        """
        Render the treatment outcome analysis section of the dashboard.
        
        Args:
            treatment_data (Dict): Dictionary with treatment outcome data
            section_title (str): Title for the section
        """
        st.header(section_title)
        
        if not treatment_data:
            st.warning("No treatment outcome data available.")
            return
        
        # Create tabs for different visualizations
        treatment_tabs = st.tabs([
            "Treatment Comparison", 
            "Adverse Events", 
            "Survival Analysis"
        ])
        
        with treatment_tabs[0]:
            st.subheader("Treatment Outcome Comparison")
            
            # Format data for visualization
            treatments_df = pd.DataFrame()
            
            # Check if data is organized by treatment type
            if isinstance(treatment_data, dict) and 'radiation' in treatment_data:
                # Data is organized by treatment type
                formatted_data = {}
                for treatment, outcomes in treatment_data.items():
                    if 'response_probability' in outcomes:
                        if 'treatment' not in formatted_data:
                            formatted_data['treatment'] = []
                        formatted_data['treatment'].append(treatment)
                        
                        for key, value in outcomes.items():
                            if key != 'adverse_event_risk' and key != 'contributing_variants' and key != 'pathway_impacts':
                                if key not in formatted_data:
                                    formatted_data[key] = []
                                formatted_data[key].append(value)
                        
                        # Handle adverse_event_risk which is a nested dict
                        if 'adverse_event_risk' in outcomes and 'overall' in outcomes['adverse_event_risk']:
                            if 'adverse_event_risk_overall' not in formatted_data:
                                formatted_data['adverse_event_risk_overall'] = []
                            formatted_data['adverse_event_risk_overall'].append(outcomes['adverse_event_risk']['overall'])
                
                if formatted_data:
                    treatments_df = pd.DataFrame(formatted_data)
            else:
                st.info("Please provide treatment outcome data in the expected format.")
            
            if not treatments_df.empty:
                # Create treatment comparison visualization
                fig = self.treatment_viz.plot_treatment_comparison(
                    treatments_data=treatments_df,
                    response_col='response_probability',
                    survival_col='survival_probability',
                    adverse_col='adverse_event_risk_overall',
                    treatment_col='treatment',
                    title="Treatment Outcome Comparison",
                    interactive=True
                )
                
                if fig is not None:
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.warning("Unable to create treatment comparison visualization.")
                
                # Show data table
                st.subheader("Treatment Data")
                st.dataframe(treatments_df, use_container_width=True)
                
                # Let user download the data
                csv_data = treatments_df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label="Download Treatment Data",
                    data=csv_data,
                    file_name="treatment_outcomes.csv",
                    mime="text/csv",
                    key="download-treatment-data"
                )
        
        with treatment_tabs[1]:
            st.subheader("Adverse Event Analysis")
            
            # Check if we have adverse event data
            has_adverse_data = False
            adverse_events_data = {}
            
            if isinstance(treatment_data, dict):
                for treatment, outcomes in treatment_data.items():
                    if 'adverse_event_risk' in outcomes and isinstance(outcomes['adverse_event_risk'], dict):
                        adverse_events_data[treatment] = outcomes['adverse_event_risk']
                        has_adverse_data = True
            
            if has_adverse_data:
                # Create adverse event visualization
                fig = self.treatment_viz.plot_detailed_adverse_events(
                    adverse_events_data=adverse_events_data,
                    title="Adverse Event Risk by Treatment Type",
                    interactive=True
                )
                
                if fig is not None:
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.warning("Unable to create adverse event visualization.")
                
                # Show data table
                st.subheader("Adverse Event Data")
                
                # Format for display
                adverse_data_list = []
                for treatment, events in adverse_events_data.items():
                    for event, risk in events.items():
                        adverse_data_list.append({
                            'Treatment': treatment,
                            'Adverse Event': event,
                            'Risk': risk
                        })
                
                adverse_df = pd.DataFrame(adverse_data_list)
                if not adverse_df.empty:
                    pivot_df = adverse_df.pivot(index='Adverse Event', columns='Treatment', values='Risk')
                    st.dataframe(pivot_df, use_container_width=True)
            else:
                st.info("No detailed adverse event data available.")
        
        with treatment_tabs[2]:
            st.subheader("Survival Analysis")
            
            # Check if we have survival curve data
            has_survival_data = False
            
            if isinstance(treatment_data, dict) and 'survival_curves' in treatment_data:
                # Display pre-computed survival curves
                survival_data = treatment_data['survival_curves']
                has_survival_data = True
            else:
                st.info("""
                To display survival curves, provide data in the format:
                ```
                {
                    "survival_curves": {
                        "time_points": [0, 6, 12, 18, 24, ...],
                        "treatments": {
                            "radiation": [1.0, 0.95, 0.85, 0.75, 0.65, ...],
                            "chemotherapy": [1.0, 0.92, 0.81, 0.70, 0.60, ...],
                            ...
                        },
                        "confidence_intervals": {
                            "radiation": [(1.0, 1.0), (0.92, 0.98), ...],
                            ...
                        }
                    }
                }
                ```
                """)
                
                # Generate example data for demo
                st.write("Generating example survival curves for demonstration:")
                
                # Create example time points and survival probabilities
                time_points = list(range(0, 61, 6))  # 0 to 60 months by 6
                
                survival_data = {
                    "time_points": time_points,
                    "treatments": {
                        "Radiation": [1.0] + [max(0, 1.0 - 0.05*i - 0.01*i**2) for i in range(1, len(time_points))],
                        "Chemotherapy": [1.0] + [max(0, 1.0 - 0.08*i - 0.005*i**2) for i in range(1, len(time_points))],
                        "Immunotherapy": [1.0] + [max(0, 1.0 - 0.03*i - 0.008*i**2) for i in range(1, len(time_points))]
                    }
                }
                has_survival_data = True
            
            if has_survival_data:
                # Extract data
                time_points = survival_data.get('time_points', [])
                survival_probs = survival_data.get('treatments', {})
                confidence_intervals = survival_data.get('confidence_intervals', {})
                
                if time_points and survival_probs:
                    # Create survival curve visualization
                    fig = self.treatment_viz.plot_survival_curve(
                        time_points=time_points,
                        survival_probs=survival_probs,
                        confidence_intervals=confidence_intervals,
                        title="Predicted Survival Curves",
                        interactive=True
                    )
                    
                    if fig is not None:
                        st.plotly_chart(fig, use_container_width=True)
                    else:
                        st.warning("Unable to create survival curve visualization.")
                else:
                    st.warning("Incomplete survival data: need time points and probabilities.")


def create_dashboard(theme: str = "light") -> Dashboard:
    """
    Create a new dashboard instance.
    
    Args:
        theme (str): Color theme for visualizations ("light" or "dark")
        
    Returns:
        Dashboard: New dashboard instance
    """
    return Dashboard(theme=theme)


def run_dashboard():
    """Run the Streamlit dashboard application."""
    st.set_page_config(
        page_title="Evo2_Clinical Dashboard",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    st.title("Evo2_Clinical Analysis Dashboard")
    st.sidebar.header("Evo2_Clinical")
    
    # Load theme preference
    theme = st.sidebar.selectbox(
        "Theme:",
        options=["Light", "Dark"],
        index=0
    ).lower()
    
    # Create dashboard instance
    dashboard = Dashboard(theme=theme)
    
    # Create tabs for different analysis types
    tabs = st.tabs(["Variant Analysis", "Pathway Analysis", "Treatment Outcomes"])
    
    with tabs[0]:
        st.sidebar.header("Variant Data")
        st.header("Variant Analysis")
        
        # Allow user to upload or select variant data
        variant_source = st.sidebar.radio(
            "Variant data source:",
            options=["Upload file", "Use example data"]
        )
        
        variants_df = pd.DataFrame()
        
        if variant_source == "Upload file":
            uploaded_file = st.sidebar.file_uploader("Upload variant data:", type=["csv", "json"])
            if uploaded_file is not None:
                try:
                    if uploaded_file.name.endswith('.csv'):
                        variants_df = pd.read_csv(uploaded_file)
                    elif uploaded_file.name.endswith('.json'):
                        variants_df = pd.read_json(uploaded_file)
                    
                    st.success(f"Loaded variant data with {len(variants_df)} variants.")
                except Exception as e:
                    st.error(f"Error loading file: {e}")
        else:
            # Generate example data
            st.info("Using example variant data")
            
            # Create synthetic data
            variants_df = pd.DataFrame({
                'id': [f"rs{i}" for i in range(10000, 10020)],
                'chrom': np.random.choice(['1', '2', '3', '4', '5', '6', 'X'], 20),
                'pos': np.random.randint(1000000, 5000000, 20),
                'ref': np.random.choice(['A', 'C', 'G', 'T'], 20),
                'alt': np.random.choice(['A', 'C', 'G', 'T'], 20),
                'gene': np.random.choice(['BRCA1', 'TP53', 'EGFR', 'KRAS', 'IL6'], 20),
                'impact_score': np.random.beta(2, 5, 20),
                'functional_score': np.random.normal(0.5, 0.2, 20).clip(0, 1),
                'conservation_score': np.random.beta(5, 2, 20)
            })
            
            # Ensure ref != alt
            for i in range(len(variants_df)):
                if variants_df.loc[i, 'ref'] == variants_df.loc[i, 'alt']:
                    available = [b for b in ['A', 'C', 'G', 'T'] if b != variants_df.loc[i, 'ref']]
                    variants_df.loc[i, 'alt'] = np.random.choice(available)
        
        # Render the variant section if we have data
        if not variants_df.empty:
            dashboard.render_variant_section(variants_df)
        else:
            st.info("Please upload variant data or select 'Use example data'.")
    
    with tabs[1]:
        st.sidebar.header("Pathway Data")
        st.header("Pathway Analysis")
        
        # Allow user to upload or select pathway data
        pathway_source = st.sidebar.radio(
            "Pathway data source:",
            options=["Upload file", "Use example data"],
            key="pathway_source_radio"
        )
        
        pathway_scores = {}
        gene_scores = None
        
        if pathway_source == "Upload file":
            uploaded_file = st.sidebar.file_uploader("Upload pathway data:", type=["json"], key="pathway_uploader")
            if uploaded_file is not None:
                try:
                    data = json.load(uploaded_file)
                    
                    if 'pathway_details' in data:
                        pathway_scores = {
                            pathway: details.get('impact_score', 0) 
                            for pathway, details in data['pathway_details'].items()
                        }
                        
                        # Extract gene scores if available
                        if 'variant_contributions' in data:
                            gene_scores = {}
                            for variant in data['variant_contributions']:
                                if 'gene' in variant and 'impact_score' in variant:
                                    gene = variant['gene']
                                    score = variant['impact_score']
                                    gene_scores[gene] = score
                    else:
                        # Simple format where the file is just a dictionary of pathway -> score
                        pathway_scores = data
                    
                    st.success(f"Loaded pathway data with {len(pathway_scores)} pathways.")
                except Exception as e:
                    st.error(f"Error loading file: {e}")
        else:
            # Generate example data
            st.info("Using example pathway data")
            
            # Create synthetic pathway data
            pathways = [
                'EMT', 'TGF_beta', 'VEGF', 'hypoxia', 'WNT', 
                'NOTCH', 'P53', 'PI3K', 'MAPK', 'JAK_STAT'
            ]
            
            pathway_scores = {
                path: np.random.uniform(-0.8, 0.8) for path in pathways
            }
            
            # Create synthetic gene data
            genes = [
                'SNAI1', 'SNAI2', 'TWIST1', 'CDH1', 'VIM', 'ZEB1', 'ZEB2',
                'TGFB1', 'SMAD2', 'SMAD3', 'SMAD4', 'CTNNB1', 'NOTCH1',
                'HIF1A', 'VEGFA', 'VEGFB', 'VEGFC', 'FLT1', 'KDR'
            ]
            
            gene_scores = {
                gene: np.random.uniform(-0.9, 0.9) for gene in genes
            }
        
        # Render the pathway section if we have data
        if pathway_scores:
            dashboard.render_pathway_section(pathway_scores, gene_scores)
        else:
            st.info("Please upload pathway data or select 'Use example data'.")
    
    with tabs[2]:
        st.sidebar.header("Treatment Outcome Data")
        st.header("Treatment Outcome Analysis")
        
        # Allow user to upload or select treatment data
        treatment_source = st.sidebar.radio(
            "Treatment data source:",
            options=["Upload file", "Use example data"],
            key="treatment_source_radio"
        )
        
        treatment_data = {}
        
        if treatment_source == "Upload file":
            uploaded_file = st.sidebar.file_uploader("Upload treatment data:", type=["json"], key="treatment_uploader")
            if uploaded_file is not None:
                try:
                    treatment_data = json.load(uploaded_file)
                    st.success(f"Loaded treatment outcome data.")
                except Exception as e:
                    st.error(f"Error loading file: {e}")
        else:
            # Generate example data
            st.info("Using example treatment outcome data")
            
            # Create synthetic treatment data
            treatment_types = ['radiation', 'chemotherapy', 'immunotherapy']
            treatment_data = {}
            
            for treatment in treatment_types:
                # Base response and survival probabilities
                base_response = 0.5
                base_survival = 0.7
                
                if treatment == 'radiation':
                    response_prob = base_response + 0.2
                    survival_prob = base_survival + 0.05
                    adverse_events = {
                        'pneumonitis': np.random.beta(2, 8),
                        'fibrosis': np.random.beta(1.5, 10),
                        'esophagitis': np.random.beta(2, 7),
                        'overall': np.random.beta(2, 5)
                    }
                elif treatment == 'chemotherapy':
                    response_prob = base_response + 0.1
                    survival_prob = base_survival - 0.1
                    adverse_events = {
                        'neutropenia': np.random.beta(3, 7),
                        'neuropathy': np.random.beta(2, 6),
                        'nausea': np.random.beta(4, 6),
                        'overall': np.random.beta(3, 6)
                    }
                else:  # immunotherapy
                    response_prob = base_response - 0.1
                    survival_prob = base_survival + 0.1
                    adverse_events = {
                        'colitis': np.random.beta(1, 9),
                        'pneumonitis': np.random.beta(1, 10),
                        'hepatitis': np.random.beta(1, 11),
                        'endocrinopathy': np.random.beta(2, 8),
                        'overall': np.random.beta(1.5, 8)
                    }
                
                treatment_data[treatment] = {
                    'response_probability': float(response_prob),
                    'survival_probability': float(survival_prob),
                    'adverse_event_risk': adverse_events,
                    'confidence_score': float(np.random.uniform(0.6, 0.9))
                }
            
            # Add survival curve data
            time_points = list(range(0, 61, 6))  # 0 to 60 months by 6
            
            survival_curves = {
                'time_points': time_points,
                'treatments': {
                    'radiation': [1.0] + [max(0, 1.0 - 0.05*i - 0.01*i**2) for i in range(1, len(time_points))],
                    'chemotherapy': [1.0] + [max(0, 1.0 - 0.08*i - 0.005*i**2) for i in range(1, len(time_points))],
                    'immunotherapy': [1.0] + [max(0, 1.0 - 0.03*i - 0.008*i**2) for i in range(1, len(time_points))]
                },
                'confidence_intervals': {
                    'radiation': [(1.0, 1.0)] + [(max(0, p-0.05), min(1.0, p+0.05)) 
                                               for p in [max(0, 1.0 - 0.05*i - 0.01*i**2) for i in range(1, len(time_points))]],
                }
            }
            
            treatment_data['survival_curves'] = survival_curves
        
        # Render the treatment section if we have data
        if treatment_data:
            dashboard.render_treatment_section(treatment_data)
        else:
            st.info("Please upload treatment outcome data or select 'Use example data'.")
    
    # Display footer
    st.markdown("---")
    st.markdown(
        "Evo2_Clinical Dashboard | "
        "Built with Streamlit | "
        "© 2025"
    )


if __name__ == "__main__":
    run_dashboard()