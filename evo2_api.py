#!/usr/bin/env python3
"""
NVIDIA Evo2 API Integration Module

This module provides functionality to interact with the NVIDIA Evo2 cloud service API
for variant scoring and genomic analysis. It handles API authentication, request/response
formatting, and error handling.

Author: Shabab Khan
Date: April 1, 2025
"""

import requests
import json
import logging
import os
import time
from typing import Dict, List, Optional, Union, Any
import pandas as pd
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("Evo2API")

class Evo2API:
    """Client for interacting with the NVIDIA Evo2 cloud API."""
    
    def __init__(self, api_key: Optional[str] = None, base_url: str = "https://api.nvidia.com/evo2/"):
        """Initialize the Evo2 API client.
        
        Args:
            api_key: NVIDIA API key for authentication. If None, will try to load from environment
            base_url: Base URL for the NVIDIA Evo2 API
        """
        # Try to get API key from environment if not provided
        self.api_key = api_key or os.environ.get("NVIDIA_EVO2_API_KEY")
        if not self.api_key:
            logger.warning("No API key provided. Set NVIDIA_EVO2_API_KEY environment variable or pass key directly.")
            
        self.base_url = base_url
        self.session = requests.Session()
        self.session.headers.update({
            "Content-Type": "application/json",
            "Authorization": f"Bearer {self.api_key}",
            "User-Agent": "Evo2Clinical/1.0"
        })
        
        logger.info(f"Initialized Evo2API client with base URL: {base_url}")
    
    def score_variants(self, variants: Union[pd.DataFrame, List[Dict]], 
                      context_window: int = 1000000,
                      genome_build: str = "GRCh38") -> pd.DataFrame:
        """Score variants using the Evo2 cloud API.
        
        Args:
            variants: DataFrame or list of variant dictionaries to score
            context_window: Size of genomic context window in bp
            genome_build: Reference genome build (GRCh37 or GRCh38)
            
        Returns:
            DataFrame with variants and their Evo2 scores
        """
        # Convert DataFrame to list of dictionaries if needed
        if isinstance(variants, pd.DataFrame):
            variant_list = variants.to_dict(orient='records')
        else:
            variant_list = variants
            
        # Format request payload
        payload = {
            "variants": variant_list,
            "parameters": {
                "context_window": context_window,
                "genome_build": genome_build
            }
        }
        
        # Make API request
        logger.info(f"Sending {len(variant_list)} variants to Evo2 API for scoring")
        try:
            response = self.session.post(
                f"{self.base_url}score",
                json=payload
            )
            response.raise_for_status()
            
            # Parse response
            result = response.json()
            
            if "error" in result:
                logger.error(f"API returned error: {result['error']}")
                return pd.DataFrame()
                
            # Convert results to DataFrame
            if isinstance(result.get("scored_variants"), list):
                return pd.DataFrame(result["scored_variants"])
            else:
                logger.error("API response missing scored_variants list")
                return pd.DataFrame()
                
        except requests.exceptions.RequestException as e:
            logger.error(f"API request failed: {str(e)}")
            return pd.DataFrame()
            
    def analyze_gene(self, gene_id: str, population: str = "ALL") -> Dict:
        """Analyze variants in a specific gene using the Evo2 cloud API.
        
        Args:
            gene_id: Ensembl or HGNC gene identifier
            population: Population code (ALL, AFR, AMR, EAS, EUR, SAS)
            
        Returns:
            Dictionary with gene analysis results
        """
        payload = {
            "gene_id": gene_id,
            "population": population
        }
        
        logger.info(f"Requesting gene analysis for {gene_id} in {population} population")
        try:
            response = self.session.post(
                f"{self.base_url}analyze/gene",
                json=payload
            )
            response.raise_for_status()
            return response.json()
            
        except requests.exceptions.RequestException as e:
            logger.error(f"Gene analysis request failed: {str(e)}")
            return {"error": str(e)}
            
    def analyze_lncrna(self, lncrna_id: str) -> Dict:
        """Analyze a long non-coding RNA using the Evo2 cloud API.
        
        Args:
            lncrna_id: Identifier for the lncRNA
            
        Returns:
            Dictionary with lncRNA analysis results
        """
        payload = {
            "lncrna_id": lncrna_id
        }
        
        logger.info(f"Requesting lncRNA analysis for {lncrna_id}")
        try:
            response = self.session.post(
                f"{self.base_url}analyze/lncrna",
                json=payload
            )
            response.raise_for_status()
            return response.json()
            
        except requests.exceptions.RequestException as e:
            logger.error(f"lncRNA analysis request failed: {str(e)}")
            return {"error": str(e)}
    
    def get_encode_data(self, region: str, cell_types: Optional[List[str]] = None) -> Dict:
        """Get ENCODE data for a genomic region from the Evo2 cloud API.
        
        Args:
            region: Genomic region in format 'chrX:start-end'
            cell_types: Optional list of cell types to filter for
            
        Returns:
            Dictionary with ENCODE data for the region
        """
        payload = {
            "region": region
        }
        
        if cell_types:
            payload["cell_types"] = cell_types
            
        logger.info(f"Requesting ENCODE data for region {region}")
        try:
            response = self.session.post(
                f"{self.base_url}data/encode",
                json=payload
            )
            response.raise_for_status()
            return response.json()
            
        except requests.exceptions.RequestException as e:
            logger.error(f"ENCODE data request failed: {str(e)}")
            return {"error": str(e)}
    
    def submit_batch_job(self, variants: List[Dict], job_name: Optional[str] = None) -> str:
        """Submit a batch scoring job to the Evo2 cloud API.
        
        Args:
            variants: List of variant dictionaries to score
            job_name: Optional name for the batch job
            
        Returns:
            Job ID string for checking status
        """
        payload = {
            "variants": variants,
            "name": job_name or f"evo2_batch_{int(time.time())}"
        }
        
        logger.info(f"Submitting batch job with {len(variants)} variants")
        try:
            response = self.session.post(
                f"{self.base_url}jobs/batch",
                json=payload
            )
            response.raise_for_status()
            result = response.json()
            
            if "job_id" in result:
                logger.info(f"Batch job submitted with ID: {result['job_id']}")
                return result["job_id"]
            else:
                logger.error("No job_id returned from API")
                return ""
                
        except requests.exceptions.RequestException as e:
            logger.error(f"Batch job submission failed: {str(e)}")
            return ""
    
    def get_job_status(self, job_id: str) -> Dict:
        """Check the status of a batch job.
        
        Args:
            job_id: ID of the job to check
            
        Returns:
            Dictionary with job status information
        """
        logger.info(f"Checking status of job {job_id}")
        try:
            response = self.session.get(
                f"{self.base_url}jobs/status/{job_id}"
            )
            response.raise_for_status()
            return response.json()
            
        except requests.exceptions.RequestException as e:
            logger.error(f"Job status check failed: {str(e)}")
            return {"status": "error", "error": str(e)}
    
    def get_job_results(self, job_id: str) -> pd.DataFrame:
        """Get results of a completed batch job.
        
        Args:
            job_id: ID of the completed job
            
        Returns:
            DataFrame with job results
        """
        logger.info(f"Retrieving results for job {job_id}")
        try:
            response = self.session.get(
                f"{self.base_url}jobs/results/{job_id}"
            )
            response.raise_for_status()
            result = response.json()
            
            if "results" in result and isinstance(result["results"], list):
                return pd.DataFrame(result["results"])
            else:
                logger.error("Invalid results format returned from API")
                return pd.DataFrame()
                
        except requests.exceptions.RequestException as e:
            logger.error(f"Job results retrieval failed: {str(e)}")
            return pd.DataFrame()


def main():
    """Example usage of the Evo2API client."""
    # This requires a valid API key to be set as NVIDIA_EVO2_API_KEY environment variable
    api = Evo2API()
    
    # Example variant for scoring
    variants = [
        {
            "chrom": "19",
            "pos": 45412079,
            "ref": "T",
            "alt": "C",
            "id": "rs429358"
        }
    ]
    
    # Score variants
    result = api.score_variants(variants)
    print(result)
    
    # Analyze a gene
    gene_result = api.analyze_gene("APOE")
    print(gene_result)


if __name__ == "__main__":
    main()