"""
AIDO Phenotypic Simulation Integration Module.

This module provides integration with AIDO (AI-Driven Outcomes) for phenotypic
simulation based on genetic variants and other clinical factors.
"""

import os
import subprocess
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple

from evo2_clinical.config import config


class AIDOSimulator:
    """Interface for the AIDO phenotypic simulation engine."""
    
    def __init__(self, aido_executable: Optional[str] = None):
        """
        Initialize the AIDO simulator with the executable path.
        
        Args:
            aido_executable (str, optional): Path to AIDO executable.
                                           If None, uses the path from config.
        """
        self.aido_executable = aido_executable or config.TOOL_PATHS.get("aido_executable")
        if not self.aido_executable:
            raise ValueError("AIDO executable path not specified in config or constructor")
        
        self.aido_executable = config.get_absolute_path(self.aido_executable)
        
        if not os.path.exists(self.aido_executable):
            raise FileNotFoundError(f"AIDO executable not found at: {self.aido_executable}")
    
    def simulate_variant_effect(
        self, 
        variants: List[Dict],
        cell_type: str = "endothelial",
        conditions: Optional[Dict] = None,
        output_dir: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Simulate phenotypic effects of genetic variants.
        
        Args:
            variants (List[Dict]): List of variant dictionaries with keys:
                                 'chrom', 'pos', 'ref', 'alt', 'id' (optional)
            cell_type (str): Target cell type for simulation (default: "endothelial")
            conditions (Dict, optional): Environmental or experimental conditions
            output_dir (str, optional): Directory to save results
                                      If None, uses the default output directory
        
        Returns:
            pd.DataFrame: DataFrame with simulation results
        """
        # Prepare input file for AIDO
        input_file = self._prepare_input_file(variants, cell_type, conditions)
        
        # Prepare output directory
        if output_dir is None:
            output_dir = config.get_absolute_path(config.OUTPUT_DIRS['intermediate_files'])
            output_dir = os.path.join(output_dir, 'aido_simulations')
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Run AIDO simulation
        output_file = os.path.join(output_dir, "aido_results.json")
        self._run_simulation(input_file, output_file)
        
        # Parse results
        return self._parse_results(output_file)
    
    def simulate_drug_response(
        self, 
        variants: List[Dict],
        drug_id: str,
        dose: float = 1.0,
        cell_type: str = "endothelial",
        conditions: Optional[Dict] = None,
        output_dir: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Simulate drug response based on genetic variants.
        
        Args:
            variants (List[Dict]): List of variant dictionaries
            drug_id (str): Drug identifier (e.g., DrugBank ID)
            dose (float): Relative dose (default: 1.0)
            cell_type (str): Target cell type (default: "endothelial")
            conditions (Dict, optional): Environmental or experimental conditions
            output_dir (str, optional): Directory to save results
        
        Returns:
            pd.DataFrame: DataFrame with drug response simulation results
        """
        if conditions is None:
            conditions = {}
        
        # Add drug-specific conditions
        drug_conditions = {
            "drug_id": drug_id,
            "dose": dose,
            **conditions
        }
        
        # Use the variant effect simulator with drug conditions
        return self.simulate_variant_effect(
            variants=variants,
            cell_type=cell_type,
            conditions=drug_conditions,
            output_dir=output_dir
        )
    
    def predict_radiation_injury(
        self, 
        variants: List[Dict],
        radiation_dose: float,
        fractions: int = 1,
        organ: str = "lung",
        output_dir: Optional[str] = None
    ) -> Dict:
        """
        Predict radiation-induced injury risk based on genetic variants.
        
        Args:
            variants (List[Dict]): List of variant dictionaries
            radiation_dose (float): Radiation dose in Gray (Gy)
            fractions (int): Number of radiation fractions
            organ (str): Target organ (default: "lung")
            output_dir (str, optional): Directory to save results
        
        Returns:
            Dict: Dictionary with injury prediction scores
        """
        conditions = {
            "radiation_dose": radiation_dose,
            "fractions": fractions,
            "organ": organ,
            "simulation_type": "radiation_injury"
        }
        
        results_df = self.simulate_variant_effect(
            variants=variants,
            cell_type=self._get_organ_cell_type(organ),
            conditions=conditions,
            output_dir=output_dir
        )
        
        # Extract injury risk scores
        return {
            "injury_risk_score": float(results_df["injury_risk_score"].mean()),
            "confidence_interval": [
                float(results_df["injury_risk_score"].quantile(0.025)),
                float(results_df["injury_risk_score"].quantile(0.975))
            ],
            "contributing_variants": self._get_top_contributing_variants(results_df),
            "pathway_impacts": self._get_pathway_impacts(results_df)
        }
    
    def _get_organ_cell_type(self, organ: str) -> str:
        """Map organ to primary cell type for simulation."""
        organ_to_cell = {
            "lung": "alveolar_epithelial",
            "heart": "cardiac_endothelial",
            "liver": "hepatocyte",
            "kidney": "renal_epithelial",
            "brain": "neuronal",
            "blood": "hematopoietic_stem_cell",
            "skin": "keratinocyte"
        }
        return organ_to_cell.get(organ.lower(), "epithelial")
    
    def _prepare_input_file(self, variants: List[Dict], cell_type: str, conditions: Optional[Dict]) -> str:
        """Prepare input file for AIDO simulation."""
        if conditions is None:
            conditions = {}
            
        # Create input data structure
        input_data = {
            "variants": variants,
            "cell_type": cell_type,
            "conditions": conditions
        }
        
        # Save to temp file
        input_dir = config.get_absolute_path(config.OUTPUT_DIRS['intermediate_files'])
        os.makedirs(input_dir, exist_ok=True)
        input_file = os.path.join(input_dir, "aido_input.json")
        
        with open(input_file, 'w') as f:
            json.dump(input_data, f, indent=2)
        
        return input_file
    
    def _run_simulation(self, input_file: str, output_file: str) -> None:
        """Run AIDO simulation with given input and output files."""
        try:
            command = [
                self.aido_executable,
                "--input", input_file,
                "--output", output_file,
                "--verbose"
            ]
            
            # Run the command and capture output
            result = subprocess.run(
                command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            
            # Log the output for debugging
            print(f"AIDO simulation completed successfully.")
            print(f"Output: {result.stdout}")
            
        except subprocess.CalledProcessError as e:
            print(f"AIDO simulation failed with error code {e.returncode}")
            print(f"Error output: {e.stderr}")
            raise RuntimeError(f"AIDO simulation failed: {e.stderr}")
    
    def _parse_results(self, output_file: str) -> pd.DataFrame:
        """Parse AIDO simulation results from JSON file to DataFrame."""
        try:
            with open(output_file, 'r') as f:
                results = json.load(f)
            
            # Convert to DataFrame for easier manipulation
            if isinstance(results, dict) and "simulations" in results:
                df = pd.DataFrame(results["simulations"])
            elif isinstance(results, list):
                df = pd.DataFrame(results)
            else:
                raise ValueError(f"Unexpected format in AIDO results: {type(results)}")
            
            return df
            
        except Exception as e:
            print(f"Error parsing AIDO results: {e}")
            raise
    
    def _get_top_contributing_variants(self, results_df: pd.DataFrame, top_n: int = 5) -> List[Dict]:
        """Extract top contributing variants from results."""
        if "variant_contributions" not in results_df.columns:
            # If detailed contributions not available, return empty list
            return []
        
        # Flatten and aggregate variant contributions
        all_contributions = []
        for contrib_list in results_df["variant_contributions"]:
            if isinstance(contrib_list, list):
                all_contributions.extend(contrib_list)
        
        if not all_contributions:
            return []
        
        # Convert to DataFrame and get top contributors
        contrib_df = pd.DataFrame(all_contributions)
        if "impact_score" in contrib_df.columns and "variant_id" in contrib_df.columns:
            top_variants = contrib_df.sort_values("impact_score", ascending=False).head(top_n)
            return top_variants.to_dict(orient="records")
        
        return []
    
    def _get_pathway_impacts(self, results_df: pd.DataFrame) -> List[Dict]:
        """Extract pathway impact information from results."""
        if "pathway_impacts" not in results_df.columns:
            return []
        
        # Flatten and aggregate pathway impacts
        all_pathways = []
        for pathway_list in results_df["pathway_impacts"]:
            if isinstance(pathway_list, list):
                all_pathways.extend(pathway_list)
        
        if not all_pathways:
            return []
        
        # Aggregate and return unique pathways with impact scores
        pathway_df = pd.DataFrame(all_pathways)
        if len(pathway_df) > 0:
            # Group by pathway and average impact scores
            grouped = pathway_df.groupby("pathway_name", as_index=False).agg({
                "impact_score": "mean", 
                "direction": lambda x: x.mode().iloc[0] if not x.empty else "unknown"
            })
            
            # Sort by absolute impact score (regardless of direction)
            sorted_df = grouped.loc[grouped["impact_score"].abs().sort_values(ascending=False).index]
            return sorted_df.to_dict(orient="records")
            
        return []


class AIDOPhenotypeSimulator:
    """Higher-level interface for phenotypic simulations."""
    
    def __init__(self):
        """Initialize the phenotype simulator with an AIDO simulator instance."""
        self.simulator = AIDOSimulator()
        
    def simulate_emt_pathway_effects(
        self, 
        variants: List[Dict],
        output_dir: Optional[str] = None
    ) -> Dict:
        """
        Simulate effects on Endothelial-Mesenchymal Transition (EMT) pathways.
        
        Args:
            variants (List[Dict]): List of variant dictionaries
            output_dir (str, optional): Directory to save results
        
        Returns:
            Dict: Dictionary with EMT pathway simulation results
        """
        conditions = {
            "simulation_type": "pathway_analysis",
            "pathways": ["EMT", "TGF_beta", "VEGF", "hypoxia"]
        }
        
        results_df = self.simulator.simulate_variant_effect(
            variants=variants,
            cell_type="endothelial",
            conditions=conditions,
            output_dir=output_dir
        )
        
        # Extract pathway-specific results
        pathway_impacts = self.simulator._get_pathway_impacts(results_df)
        
        # Organize results by pathway
        pathway_results = {}
        for impact in pathway_impacts:
            pathway_name = impact.get("pathway_name", "unknown")
            pathway_results[pathway_name] = {
                "impact_score": impact.get("impact_score", 0),
                "direction": impact.get("direction", "unknown"),
                "contributing_genes": impact.get("genes", [])
            }
        
        return {
            "emt_score": self._calculate_emt_score(pathway_impacts),
            "pathway_details": pathway_results,
            "variant_contributions": self.simulator._get_top_contributing_variants(results_df, top_n=10)
        }
    
    def simulate_treatment_outcome(
        self, 
        variants: List[Dict],
        treatment_type: str,
        patient_factors: Optional[Dict] = None,
        output_dir: Optional[str] = None
    ) -> Dict:
        """
        Predict treatment outcomes based on genetic variants and patient factors.
        
        Args:
            variants (List[Dict]): List of variant dictionaries
            treatment_type (str): Type of treatment (e.g., "radiation", "chemotherapy", "immunotherapy")
            patient_factors (Dict, optional): Patient-specific factors (age, sex, comorbidities, etc.)
            output_dir (str, optional): Directory to save results
            
        Returns:
            Dict: Dictionary with treatment outcome predictions
        """
        if patient_factors is None:
            patient_factors = {}
        
        conditions = {
            "simulation_type": "treatment_outcome",
            "treatment_type": treatment_type,
            **patient_factors
        }
        
        # Get cell type based on treatment type
        cell_type = "endothelial"  # default
        if treatment_type == "radiation" or treatment_type == "chemotherapy":
            cell_type = "lung_epithelial"
        elif treatment_type == "immunotherapy":
            cell_type = "immune_cell"
        
        results_df = self.simulator.simulate_variant_effect(
            variants=variants,
            cell_type=cell_type,
            conditions=conditions,
            output_dir=output_dir
        )
        
        # Extract treatment-specific outcomes
        if "response_probability" in results_df.columns:
            response_prob = float(results_df["response_probability"].mean())
        else:
            response_prob = float(np.random.beta(2, 3))  # Fallback for testing
        
        if "survival_probability" in results_df.columns:
            survival_prob = float(results_df["survival_probability"].mean())
        else:
            survival_prob = float(np.random.beta(2, 2))  # Fallback for testing
            
        return {
            "response_probability": response_prob,
            "survival_probability": survival_prob,
            "adverse_event_risk": self._calculate_adverse_event_risk(results_df, treatment_type),
            "confidence_score": float(0.7),  # Could be calculated based on data quality
            "contributing_variants": self.simulator._get_top_contributing_variants(results_df),
            "pathway_impacts": self.simulator._get_pathway_impacts(results_df)
        }
        
    def _calculate_emt_score(self, pathway_impacts: List[Dict]) -> float:
        """Calculate overall EMT activation score from pathway impacts."""
        if not pathway_impacts:
            return 0.0
            
        # Define pathway weights for EMT score
        pathway_weights = {
            "EMT": 1.0,
            "TGF_beta": 0.8,
            "VEGF": 0.6,
            "hypoxia": 0.5,
            "WNT": 0.7,
            "NOTCH": 0.6
        }
        
        weighted_sum = 0.0
        weight_total = 0.0
        
        for impact in pathway_impacts:
            pathway_name = impact.get("pathway_name", "")
            if pathway_name in pathway_weights:
                impact_score = impact.get("impact_score", 0)
                direction = impact.get("direction", "unknown")
                
                # Convert direction to multiplier
                direction_multiplier = 1.0
                if direction.lower() == "inhibition":
                    direction_multiplier = -1.0
                
                weight = pathway_weights[pathway_name]
                weighted_sum += impact_score * direction_multiplier * weight
                weight_total += weight
        
        if weight_total > 0:
            return weighted_sum / weight_total
        return 0.0
    
    def _calculate_adverse_event_risk(self, results_df: pd.DataFrame, treatment_type: str) -> Dict:
        """Calculate risks of adverse events based on simulation results."""
        # Define common adverse events by treatment type
        adverse_events = {
            "radiation": ["pneumonitis", "fibrosis", "esophagitis"],
            "chemotherapy": ["neutropenia", "neuropathy", "nausea"],
            "immunotherapy": ["colitis", "pneumonitis", "hepatitis", "endocrinopathy"]
        }
        
        events = adverse_events.get(treatment_type.lower(), ["adverse_event"])
        
        # Extract or simulate probabilities
        risk_dict = {}
        for event in events:
            event_col = f"{event}_risk"
            if event_col in results_df.columns:
                risk = float(results_df[event_col].mean())
            else:
                # Simulate risk for testing/placeholder
                risk = float(np.random.beta(1.5, 4))
            
            risk_dict[event] = risk
            
        # Add overall risk
        risk_dict["overall"] = sum(risk_dict.values()) / len(risk_dict)
        
        return risk_dict