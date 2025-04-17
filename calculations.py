"""
Python Package: Molecular Analysis of CTEPH and lncRNAs

Description:
This package provides a set of tools and modules to analyze genomic and transcriptomic data,
with a focus on Chronic Thromboembolic Pulmonary Hypertension (CTEPH) and the role of
long non-coding RNAs (lncRNAs) in disease mechanisms. It incorporates functionalities for
GWAS data analysis, RNA sequence analysis, and disease association studies.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy.stats import fisher_exact, mannwhitneyu
import warnings
import pandas as pd

# Suppress future warnings, especially from statsmodels
warnings.simplefilter(action='ignore', category=FutureWarning)

# =============================================================================
# 1. GWAS Data Analysis Module
# =============================================================================

class GWASAnalyzer:
    """
    Analyzes Genome-Wide Association Study (GWAS) data to identify significant genetic variants
    associated with a disease, particularly CTEPH.

    Key Functionality:
    - Loads and preprocesses GWAS summary statistics data.
    - Filters variants based on significance thresholds (p-value).
    - Identifies variants located within or near specified genomic regions (e.g., lncRNA loci).
    - Performs basic statistical analysis (e.g., allele frequency calculation).
    """
    def __init__(self, gwas_data: pd.DataFrame):
        """
        Initializes the GWASAnalyzer with GWAS data.

        Args:
            gwas_data (pd.DataFrame): Pandas DataFrame containing GWAS summary statistics.
                Required columns: 'SNP' (variant identifier), 'CHR' (chromosome),
                'POS' (position), 'P' (p-value), 'A1' (allele 1), 'A2' (allele 2).
        """
        self.gwas_data = gwas_data.copy()  # Create a copy to avoid modifying the original DataFrame
        self._validate_gwas_data()

    def _validate_gwas_data(self):
        """
        Validates the input GWAS data for required columns.
        Raises ValueError if any required column is missing.
        """
        required_columns = ['SNP', 'CHR', 'POS', 'P', 'A1', 'A2']
        missing_columns = [col for col in required_columns if col not in self.gwas_data.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns in GWAS data: {', '.join(missing_columns)}")

    def filter_significant_variants(self, p_threshold: float = 0.05) -> pd.DataFrame:
        """
        Filters GWAS variants based on a p-value threshold.

        Args:
            p_threshold (float, optional): Significance threshold for the p-value. Defaults to 0.05.

        Returns:
            pd.DataFrame: DataFrame containing variants with p-values below the threshold.
        """
        significant_variants = self.gwas_data[self.gwas_data['P'] < p_threshold].copy()
        return significant_variants

    def identify_variants_in_region(self, chromosome: str, start: int, end: int) -> pd.DataFrame:
        """
        Identifies GWAS variants located within a specified genomic region.

        Args:
            chromosome (str): Chromosome name (e.g., '1', 'chr1', 'X').
            start (int): Start position of the region.
            end (int): End position of the region.

        Returns:
            pd.DataFrame: DataFrame containing variants within the specified region.
        """
        region_variants = self.gwas_data[
            (self.gwas_data['CHR'] == chromosome) &
            (self.gwas_data['POS'] >= start) &
            (self.gwas_data['POS'] <= end)
            ].copy()
        return region_variants

    def calculate_allele_frequencies(self, population_allele_freq: Optional[Dict[str, float]] = None) -> pd.DataFrame:
        """
        Calculates allele frequencies from the GWAS data.  If population allele
        frequencies are provided, it adds a column 'MAF' (Minor Allele Frequency)
        representing the *population* MAF.  If not, it calculates MAF from the GWAS
        sample itself.

        Args:
            population_allele_freq (Optional[Dict[str, float]], optional):
                Dictionary of population allele frequencies, keyed by SNP identifier.
                Defaults to None.

        Returns:
            pd.DataFrame: DataFrame with allele frequencies.  If population_allele_freq
                is provided, includes a 'MAF' column.
        """

        def _calculate_maf(row):
            """Calculates minor allele frequency from allele counts."""
            alleles = [row['A1'], row['A2']]
            counts = Counter(alleles)
            return min(counts.values()) / sum(counts.values())

        gwas_data_with_freq = self.gwas_data.copy() # Avoid modifying original
        if population_allele_freq:
            gwas_data_with_freq['MAF'] = gwas_data_with_freq['SNP'].map(population_allele_freq)
            gwas_data_with_freq['MAF'] = gwas_data_with_freq['MAF'].fillna(0) # Assume 0 if not found.
        else:
            gwas_data_with_freq['MAF'] = gwas_data_with_freq.apply(_calculate_maf, axis=1)
        return gwas_data_with_freq

    def perform_regression_analysis(self, phenotype_data: pd.DataFrame,
                                   covariates: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Performs regression analysis to assess the association between genetic variants
        and a quantitative phenotype.

        Args:
            phenotype_data (pd.DataFrame): DataFrame containing phenotype data.
                Required columns: 'SUBJECT_ID' (subject identifier), 'PHENOTYPE' (phenotype value).
            covariates (Optional[List[str]], optional): List of covariate names to include in the model.
                Defaults to None.  These covariates must be columns in *both*
                `self.gwas_data` *and* `phenotype_data`.

        Returns:
            pd.DataFrame: DataFrame with regression analysis results (beta, standard error, p-value).
        """
        # Merge GWAS and phenotype data on a common identifier (e.g., subject ID)
        merged_data = pd.merge(self.gwas_data, phenotype_data, left_on='SUBJECT_ID', right_on='SUBJECT_ID')

        # Prepare the regression formula
        if covariates:
            formula = f"PHENOTYPE ~ SNP + {' + '.join(covariates)}"
            # Check that covariates are in the merged_data
            missing_covariates = [c for c in covariates if c not in merged_data.columns]
            if missing_covariates:
                raise ValueError(f"Covariates not found in merged data: {', '.join(missing_covariates)}")
        else:
            formula = "PHENOTYPE ~ SNP"

        # Perform the regression analysis
        model = ols(formula, merged_data).fit()
        results_df = pd.DataFrame({
            'SNP': [model.params.index[1]],  # Assuming SNP is the second term in the model
            'beta': [model.params[1]],
            'std_err': [model.bse[1]],
            'p_value': [model.pvalues[1]]
        })
        return results_df

# =============================================================================
# 2. RNA Sequence Analysis Module
# =============================================================================

class RNAAnalyzer:
    """
    Analyzes RNA sequence data, including lncRNA sequences.

    Key Functionality:
    - Loads RNA sequence data from FASTA files.
    - Calculates basic sequence statistics (e.g., length, GC content).
    - (Placeholder) Predicts RNA secondary structure (integration with external tools may be needed).
    - (Placeholder) Predicts lncRNA essentiality.
    """
    def __init__(self, fasta_file: str):
        """
        Initializes the RNAAnalyzer with the path to an RNA sequence FASTA file.

        Args:
            fasta_file (str): Path to the FASTA file containing RNA sequences.
        """
        self.fasta_file = fasta_file
        self.sequences = self._load_sequences()

    def _load_sequences(self) -> List[SeqRecord]:
        """
        Loads RNA sequences from the FASTA file.

        Returns:
            List[SeqRecord]: List of SeqRecord objects representing the RNA sequences.
        Raises:
            FileNotFoundError: If the fasta file does not exist
        """
        try:
            sequences = list(SeqIO.parse(self.fasta_file, "fasta"))
            return sequences
        except FileNotFoundError:
            raise FileNotFoundError(f"Fasta file not found at {self.fasta_file}")

    def calculate_sequence_statistics(self) -> pd.DataFrame:
        """
        Calculates sequence length and GC content for each RNA sequence.

        Returns:
            pd.DataFrame: DataFrame containing sequence statistics (ID, length, GC content).
        """
        data = []
        for record in self.sequences:
            sequence = str(record.seq).upper()  # Ensure uppercase for GC calculation
            length = len(sequence)
            gc_content = (sequence.count('G') + sequence.count('C')) / length * 100 if length > 0 else 0
            data.append({'ID': record.id, 'length': length, 'GC_content': gc_content})
        return pd.DataFrame(data)

    def predict_rna_structure(self, sequence_id: str) -> str:
        """
        (Placeholder) Predicts the secondary structure of an RNA sequence.
        This is a placeholder; a real implementation would likely involve calling
        an external tool or web service (e.g., RNAfold, ViennaRNA).

        Args:
            sequence_id (str): ID of the RNA sequence to analyze.

        Returns:
            str: Predicted RNA secondary structure (placeholder, currently returns "No prediction").
        """
        # Placeholder:  In a real implementation, you would:
        # 1.  Extract the sequence using the sequence_id.
        # 2.  Call an external program or web service (e.g., via subprocess.run()
        #      or a library like requests) to predict the structure.
        # 3.  Parse the output of the external tool and return the structure.
        #
        #  Example using a hypothetical external tool (replace with actual tool):
        #  try:
        #      result = subprocess.run(['RNAfold', '-i', sequence], capture_output=True, text=True, check=True)
        #      structure = result.stdout.strip()
        #      return structure
        #  except subprocess.CalledProcessError as e:
        #      print(f"Error predicting structure: {e}")
        #      return "Error"
        return "No prediction"  # Placeholder

    def predict_lncrna_essentiality(self, sequence_id: str) -> float:
        """
        (Placeholder) Predicts the essentiality of an lncRNA.  This is a placeholder;
        a real implementation would require a sophisticated model, potentially
        based on sequence features, structural information, and comparative genomics.

        Args:
            sequence_id (str): ID of the lncRNA sequence to analyze.

        Returns:
            float: Predicted essentiality score (placeholder, range 0-1, higher = more essential).
        """
        # Placeholder:  A real implementation would involve:
        # 1.  Extracting the sequence.
        # 2.  Calculating relevant features (e.g., conservation, structural features).
        # 3.  Applying a pre-trained model (e.g., a machine learning model)
        #      to predict essentiality.
        return 0.5  # Placeholder (return a default value)

# =============================================================================
# 3. Disease Association Module
# =============================================================================

class DiseaseAssociationAnalyzer:
    """
    Analyzes the association between genetic variants/lncRNA features and disease phenotypes,
    specifically focusing on CTEPH.

    Key Functionality:
    - Integrates GWAS data with phenotype data.
    - Performs statistical tests to assess associations (e.g., logistic regression, t-tests).
    - (Placeholder) Implements Bayesian Prioritisation for variant/lncRNA ranking.
    """
    def __init__(self, gwas_analyzer: GWASAnalyzer, rna_analyzer: RNAAnalyzer, phenotype_data: pd.DataFrame):
        """
        Initializes the DiseaseAssociationAnalyzer with GWAS data, RNA sequence data, and phenotype data.

        Args:
            gwas_analyzer (GWASAnalyzer):  Pre-initialized GWASAnalyzer object.
            rna_analyzer (RNAAnalyzer): Pre-initialized RNAAnalyzer object.
            phenotype_data (pd.DataFrame): DataFrame containing phenotype data.
                Required columns: 'SUBJECT_ID' (subject identifier), 'DISEASE_STATUS' (binary: 0=control, 1=case).
        """
        self.gwas_analyzer = gwas_analyzer
        self.rna_analyzer = rna_analyzer # currently not used, but kept for potential future use
        self.phenotype_data = phenotype_data.copy() # Avoid modifying
        self._validate_phenotype_data()

    def _validate_phenotype_data(self):
        """
        Validates the input phenotype data.
        Raises ValueError if required columns are missing or if the disease status
        is not binary.
        """
        required_columns = ['SUBJECT_ID', 'DISEASE_STATUS']
        missing_columns = [col for col in required_columns if col not in self.phenotype_data.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns in phenotype data: {', '.join(missing_columns)}")
        if not set(self.phenotype_data['DISEASE_STATUS']).issubset({0, 1}):
            raise ValueError("DISEASE_STATUS must be binary (0 or 1)")

    def assess_variant_association(self, p_threshold: float = 0.05) -> pd.DataFrame:
        """
        Assesses the association between significant GWAS variants and disease status (CTEPH)
        using logistic regression.

        Args:
            p_threshold (float, optional): Significance threshold for GWAS variants. Defaults to 0.05.

        Returns:
            pd.DataFrame: DataFrame containing results of the logistic regression analysis
                         for each significant variant.
        """
        significant_variants = self.gwas_analyzer.filter_significant_variants(p_threshold)
        if significant_variants.empty:
            print("No significant variants found based on the given p-value threshold.")
            return pd.DataFrame()  # Return empty DataFrame

        # Merge GWAS data with phenotype data
        merged_data = pd.merge(significant_variants, self.phenotype_data, on='SUBJECT_ID')

        results = []
        for index, row in significant_variants.iterrows():
            snp = row['SNP']
            # Ensure that the SNP column exists in the merged data.
            if snp not in merged_data.columns:
                print(f"SNP {snp} not found in merged data. Skipping.")
                continue
            # Prepare data for the current SNP
            X = merged_data[[snp]]  # Use double brackets to select a DataFrame
            X = sm.add_constant(X)  # Add a constant term for the intercept
            y = merged_data['DISEASE_STATUS']

            # Perform logistic regression
            model = sm.Logit(y, X)
            result = model.fit(disp=False)  # disp=False suppresses verbose output

            # Extract results
            results.append({
                'SNP': snp,
                'odds_ratio': result.params[snp],
                'p_value': result.pvalues[snp],
                'conf_lower': result.conf_int().loc[snp, 0],
                'conf_upper': result.conf_int().loc[snp, 1]
            })
        return pd.DataFrame(results)

    def compare_lncrna_expression(self, expression_data: pd.DataFrame,
                                  lncrna_id: str,
                                  group_col: str = 'DISEASE_STATUS') -> pd.DataFrame:
        """
        Compares the expression levels of a specified lncRNA between different groups
        (e.g., disease vs. control) using the Mann-Whitney U test.

        Args:
            expression_data (pd.DataFrame): DataFrame containing lncRNA expression data.
                Required columns: 'SUBJECT_ID', and a column for the lncRNA (e.g., 'ANRIL').
            lncrna_id (str): The ID of the lncRNA to compare (e.g., 'ANRIL').  This should
                           match a column name in `expression_data`.
            group_col (str, optional): The column name in `expression_data` that indicates the group
                             (e.g., 'DISEASE_STATUS'). Defaults to 'DISEASE_STATUS'.

        Returns:
            pd.DataFrame: DataFrame containing the results of the Mann-Whitney U test.
        Raises:
            ValueError: If the lncrna_id or group_col is not found in expression_data.
        """
        if lncrna_id not in expression_data.columns:
            raise ValueError(f"lncRNA ID '{lncrna_id}' not found in expression data.")
        if group_col not in expression_data.columns:
            raise ValueError(f"Group column '{group_col}' not found in expression data.")

        # Separate expression data for the two groups
        group1_data = expression_data[expression_data[group_col] == 0][lncrna_id]
        group2_data = expression_data[expression_data[group_col] == 1][lncrna_id]

        # Perform Mann-Whitney U test
        u_statistic, p_value = mannwhitneyu(group1_data, group2_data)

        # Calculate group means for context
        mean_group1 = group1_data.mean()
        mean_group2 = group2_data.mean()

        results_df = pd.DataFrame({
            'lncrna_id': [lncrna_id],
            'u_statistic': [u_statistic],
            'p_value': [p_value],
            'mean_group_0': [mean_group1],  # Added for clarity
            'mean_group_1': [mean_group2]   # Added for clarity
        })
        return results_df

    def bayesian_prioritization(self, variant_lncrna_pairs: List[Dict[str, str]],
                               priors: Dict[str, float]) -> pd.DataFrame:
        """
        (Placeholder) Performs Bayesian prioritization of variant-lncRNA pairs based on
        prior probabilities (priors).  This is a placeholder; a real implementation
        would require a more detailed model and potentially external data sources.

        Args:
            variant_lncrna_pairs (List[Dict[str, str]]): List of dictionaries, where each dictionary
                represents a variant-lncRNA pair.  Each dictionary should have keys
                'variant_id' and 'lncrna_id'.
            priors (Dict[str, float]): Dictionary of prior probabilities for each variant-lncRNA pair,
                keyed by a unique identifier for the pair (e.g., a combination of
                variant_id and lncrna_id).

        Returns:
            pd.DataFrame: DataFrame containing variant-lncRNA pairs and their posterior probabilities.
        """
        # Placeholder:  A real implementation would involve:
        # 1.  Calculating likelihoods for each pair based on available data
        #     (e.g., GWAS results, expression data, functional predictions).
        # 2.  Applying Bayes' theorem to calculate posterior probabilities:
        #        P(pair | data) = [P(data | pair) * P(pair)] / P(data)
        # 3.  Ranking the pairs based on their posterior probabilities.

        results = []
        for pair in variant_lncrna_pairs:
            variant_id = pair['variant_id']
            lncrna_id = pair['lncrna_id']
            pair_id = f"{variant_id}-{lncrna_id}"  # Create a unique ID
            prior = priors.get(pair_id, 0.01)  # Default prior if not found
            # Placeholder likelihood (replace with actual calculation)
            likelihood = 0.5  # Example likelihood
            # Placeholder evidence (replace with actual calculation)
            evidence = 0.1 # Example evidence
            posterior = (likelihood * prior) / evidence # Bayes Theorem
            results.append({
                'variant_id': variant_id,
                'lncrna_id': lncrna_id,
                'prior': prior,
                'likelihood': likelihood, # added likelihood
                'posterior': posterior
            })
        return pd.DataFrame(results).sort_values(by='posterior', ascending=False)

if __name__ == "__main__":

    # Minimal example GWAS data
    # example_gwas_df = pd.DataFrame({
    #     "SNP": ["rs123", "rs456"],
    #     "CHR": ["1", "1"],
    #     "POS": [10100, 20200],
    #     "P": [0.01, 0.15],
    #     "A1": ["A", "C"],
    #     "A2": ["G", "T"]
    # })
    example_gwas_df = pd.read_csv("data/gwas/gwas-catalog-associations.tsv", sep="\t")  # Load example GWAS data from a TSV file
    print(example_gwas_df.columns)
    example_gwas_df = example_gwas_df.loc[:, ["SNPS", "CHR_ID", "CHR_POS", "P-VALUE", "DISEASE/TRAIT", "PUBMEDID", "STRONGEST SNP-RISK ALLELE"]]  # Select relevant columns
    # rename to standard names
    example_gwas_df.rename(columns={"SNPS": "SNP", "CHR_ID": "CHR", "CHR_POS": "POS", "P-VALUE": "P", "STRONGEST SNP-RISK ALLELE": "A1"}, inplace=True)
    example_gwas_df["A2"] = example_gwas_df["A1"].apply(lambda x: "G" if x == "A" else "A")  # Dummy allele 2 for example
    print(example_gwas_df.head())  # Display the first few rows of the GWAS data

    # Instantiate the class (example)
    gwas_analyzer = GWASAnalyzer(example_gwas_df)
    significant = gwas_analyzer.filter_significant_variants(1e-93)
    print("Significant variants:\n", significant)