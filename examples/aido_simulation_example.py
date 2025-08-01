#!/usr/bin/env python
"""
Example script demonstrating the use of AIDO phenotypic simulation capabilities.

This example shows how to:
1. Simulate variant effects on EMT pathways
2. Predict radiation injury risk based on variants
3. Simulate treatment outcomes

Note: This script assumes that the AIDO executable path is properly set in the configuration.
"""

import os
import json
from pathlib import Path
import pandas as pd

from evo2_clinical.analysis import AIDOSimulator, AIDOPhenotypeSimulator
from evo2_clinical.config import Config

# Load configuration
config_path = Path(__file__).parent.parent / "evo2_clinical_config.yaml"
config = Config(str(config_path))

def load_test_variants():
    """Load example variants for demonstration."""
    # These are example variants - replace with real variants from your data
    return [
        {"chrom": "1", "pos": 11856378, "ref": "G", "alt": "A", "id": "rs12345"},
        {"chrom": "2", "pos": 21228720, "ref": "C", "alt": "T", "id": "rs67890"},
        {"chrom": "7", "pos": 55249063, "ref": "T", "alt": "G", "id": "rs55554"},
        {"chrom": "17", "pos": 7676154, "ref": "G", "alt": "C", "id": "rs77777"},
    ]

def example_1_basic_simulation():
    """Example 1: Basic variant effect simulation."""
    print("\n=== Example 1: Basic Variant Effect Simulation ===")
    
    # Initialize the AIDO simulator
    try:
        simulator = AIDOSimulator()
        print("AIDO simulator initialized successfully.")
    except ValueError as e:
        print(f"Error: {e}")
        print("Make sure AIDO executable path is set in your configuration.")
        return
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Make sure AIDO executable exists at the configured path.")
        return
        
    # Load test variants
    variants = load_test_variants()
    print(f"Loaded {len(variants)} test variants.")
    
    # Set up simulation parameters
    cell_type = "endothelial"
    conditions = {
        "temperature": 37.0,  # Celsius
        "oxygen_level": "normoxia",  # or "hypoxia"
        "growth_factors": ["VEGF", "FGF"],
    }
    
    # Run simulation
    print(f"Running simulation for {cell_type} cells...")
    try:
        results = simulator.simulate_variant_effect(
            variants=variants,
            cell_type=cell_type,
            conditions=conditions
        )
        print("Simulation completed successfully.")
        print(f"Result shape: {results.shape}")
        print("\nFirst few rows:")
        print(results.head())
    except Exception as e:
        print(f"Simulation failed: {e}")

def example_2_emt_pathway_analysis():
    """Example 2: EMT pathway analysis."""
    print("\n=== Example 2: EMT Pathway Analysis ===")
    
    # Initialize the phenotype simulator
    phenotype_sim = AIDOPhenotypeSimulator()
    print("AIDO phenotype simulator initialized.")
    
    # Load test variants
    variants = load_test_variants()
    
    # Run EMT pathway simulation
    print("Simulating effects on EMT pathways...")
    try:
        results = phenotype_sim.simulate_emt_pathway_effects(variants)
        
        print("\nEMT Activation Score:", results["emt_score"])
        print("\nPathway Details:")
        for pathway, details in results["pathway_details"].items():
            print(f"  {pathway}:")
            print(f"    Impact Score: {details['impact_score']:.3f}")
            print(f"    Direction: {details['direction']}")
        
        print("\nTop Contributing Variants:")
        for i, variant in enumerate(results["variant_contributions"][:3]):
            print(f"  {i+1}. {variant.get('variant_id', 'unknown')}: {variant.get('impact_score', 0):.3f}")
    except Exception as e:
        print(f"EMT pathway analysis failed: {e}")

def example_3_radiation_injury_prediction():
    """Example 3: Radiation injury prediction."""
    print("\n=== Example 3: Radiation Injury Prediction ===")
    
    # Initialize the AIDO simulator
    simulator = AIDOSimulator()
    
    # Load test variants
    variants = load_test_variants()
    
    # Set radiation parameters
    radiation_dose = 20.0  # Gy
    fractions = 5
    organ = "lung"
    
    # Run prediction
    print(f"Predicting radiation injury for {radiation_dose}Gy in {fractions} fractions to {organ}...")
    try:
        results = simulator.predict_radiation_injury(
            variants=variants,
            radiation_dose=radiation_dose,
            fractions=fractions,
            organ=organ
        )
        
        print("\nInjury Risk Score:", f"{results['injury_risk_score']:.3f}")
        print(f"95% Confidence Interval: [{results['confidence_interval'][0]:.3f}, {results['confidence_interval'][1]:.3f}]")
        
        print("\nContributing Variants:")
        for i, variant in enumerate(results["contributing_variants"][:3]):
            print(f"  {i+1}. {variant.get('variant_id', 'unknown')}: {variant.get('impact_score', 0):.3f}")
        
        print("\nAffected Pathways:")
        for i, pathway in enumerate(results["pathway_impacts"][:3]):
            print(f"  {i+1}. {pathway.get('pathway_name', 'unknown')}: {pathway.get('impact_score', 0):.3f} ({pathway.get('direction', 'unknown')})")
    except Exception as e:
        print(f"Radiation injury prediction failed: {e}")

def example_4_treatment_outcome_prediction():
    """Example 4: Treatment outcome prediction."""
    print("\n=== Example 4: Treatment Outcome Prediction ===")
    
    # Initialize the phenotype simulator
    phenotype_sim = AIDOPhenotypeSimulator()
    
    # Load test variants
    variants = load_test_variants()
    
    # Set treatment parameters
    treatment_type = "radiation"
    patient_factors = {
        "age": 65,
        "sex": "female",
        "smoking_history": True,
        "comorbidities": ["hypertension"],
    }
    
    # Run prediction
    print(f"Predicting outcomes for {treatment_type} treatment...")
    try:
        results = phenotype_sim.simulate_treatment_outcome(
            variants=variants,
            treatment_type=treatment_type,
            patient_factors=patient_factors
        )
        
        print("\nTreatment Response Probability:", f"{results['response_probability']:.3f}")
        print("Survival Probability:", f"{results['survival_probability']:.3f}")
        print("\nAdverse Event Risks:")
        
        for event, risk in results["adverse_event_risk"].items():
            if event != "overall":
                print(f"  {event.capitalize()}: {risk:.3f}")
        print(f"  Overall Risk: {results['adverse_event_risk']['overall']:.3f}")
        
        print(f"\nConfidence Score: {results['confidence_score']:.3f}")
    except Exception as e:
        print(f"Treatment outcome prediction failed: {e}")

def save_results_example():
    """Example: Saving simulation results to file."""
    print("\n=== Example: Saving Simulation Results ===")
    
    # Initialize the phenotype simulator
    phenotype_sim = AIDOPhenotypeSimulator()
    
    # Load test variants
    variants = load_test_variants()
    
    # Create output directory
    output_dir = Path("examples/output")
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Run simulation for different treatment types
    treatment_types = ["radiation", "chemotherapy", "immunotherapy"]
    
    results = {}
    
    for treatment in treatment_types:
        print(f"Simulating {treatment} treatment outcomes...")
        try:
            outcome = phenotype_sim.simulate_treatment_outcome(
                variants=variants,
                treatment_type=treatment
            )
            results[treatment] = outcome
            
            # Create summary for this treatment
            summary = {
                "treatment": treatment,
                "response_probability": outcome["response_probability"],
                "survival_probability": outcome["survival_probability"],
                "adverse_event_risk_overall": outcome["adverse_event_risk"]["overall"],
                "confidence": outcome["confidence_score"]
            }
            
            # Save individual result
            with open(output_dir / f"{treatment}_results.json", "w") as f:
                json.dump(outcome, f, indent=2)
                
        except Exception as e:
            print(f"Error simulating {treatment} treatment: {e}")
    
    # Create comparison table
    if results:
        comparison_data = []
        for treatment, outcome in results.items():
            comparison_data.append({
                "Treatment": treatment.capitalize(),
                "Response Probability": f"{outcome['response_probability']:.3f}",
                "Survival Probability": f"{outcome['survival_probability']:.3f}",
                "Adverse Event Risk": f"{outcome['adverse_event_risk']['overall']:.3f}"
            })
        
        comparison_df = pd.DataFrame(comparison_data)
        print("\nTreatment Comparison:")
        print(comparison_df)
        
        # Save comparison as CSV
        comparison_df.to_csv(output_dir / "treatment_comparison.csv", index=False)
        print(f"Results saved to {output_dir}")

if __name__ == "__main__":
    print("AIDO Simulation Examples")
    print("=======================")
    
    # Check if AIDO path is configured
    if not config.TOOL_PATHS.get("aido_executable"):
        print("Warning: AIDO executable path is not set in the configuration.")
        print("Setting a placeholder path for demonstration purposes.")
        # For demonstration only, set a dummy path
        config.TOOL_PATHS["aido_executable"] = "/tmp/aido"
        with open("/tmp/aido", "w") as f:
            f.write("#!/bin/sh\necho '{\"simulations\": [{\"injury_risk_score\": 0.42, \"pathway_impacts\": [{\"pathway_name\": \"EMT\", \"impact_score\": 0.8, \"direction\": \"activation\"}]}]}' > $4")
        os.chmod("/tmp/aido", 0o755)
    
    # Run examples
    try:
        example_1_basic_simulation()
        example_2_emt_pathway_analysis()
        example_3_radiation_injury_prediction()
        example_4_treatment_outcome_prediction()
        save_results_example()
    except Exception as e:
        print(f"Example execution failed: {e}")
    
    print("\nExamples completed. Check output directory for saved results.")