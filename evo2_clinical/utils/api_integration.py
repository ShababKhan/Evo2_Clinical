"""
API Integration Module.

This module provides automated data fetching from external biomedical databases
including ENCODE, GWAS Catalog, OpenTargets Platform, and 1000 Genomes Project.
"""

import os
import json
import time
import logging
import requests
from typing import Dict, List, Optional, Union, Any, Tuple
import pandas as pd
from pathlib import Path
from urllib.parse import urljoin

from evo2_clinical.config import config


class APIClient:
    """Base class for API clients."""
    
    def __init__(self, base_url: str, api_key: Optional[str] = None):
        """
        Initialize the API client.
        
        Args:
            base_url (str): Base URL for the API
            api_key (str, optional): API key for authentication
        """
        self.base_url = base_url
        self.api_key = api_key
        self.session = requests.Session()
        
        # Set up logging
        self.logger = logging.getLogger(__name__)
        
        # Configure session with default headers
        self.session.headers.update({
            "Content-Type": "application/json",
            "Accept": "application/json"
        })
        
        if api_key:
            self.session.headers.update({"Authorization": f"Bearer {api_key}"})
    
    def _handle_response(self, response: requests.Response) -> Dict:
        """
        Handle API response and raise appropriate exceptions.
        
        Args:
            response (requests.Response): Response object from requests
            
        Returns:
            Dict: JSON response data
            
        Raises:
            requests.HTTPError: If the response status code indicates an error
        """
        try:
            response.raise_for_status()
            return response.json()
        except requests.HTTPError as e:
            error_info = {}
            try:
                error_info = response.json()
            except:
                error_info = {"text": response.text}
            
            self.logger.error(f"API error: {e} - {error_info}")
            raise
        except json.JSONDecodeError:
            self.logger.error("Failed to decode JSON response")
            return {"text": response.text}
    
    def get(self, endpoint: str, params: Optional[Dict] = None) -> Dict:
        """
        Make a GET request to the API.
        
        Args:
            endpoint (str): API endpoint
            params (Dict, optional): Query parameters
            
        Returns:
            Dict: Response data
        """
        url = urljoin(self.base_url, endpoint)
        response = self.session.get(url, params=params)
        return self._handle_response(response)
    
    def post(self, endpoint: str, data: Dict) -> Dict:
        """
        Make a POST request to the API.
        
        Args:
            endpoint (str): API endpoint
            data (Dict): Request data
            
        Returns:
            Dict: Response data
        """
        url = urljoin(self.base_url, endpoint)
        response = self.session.post(url, json=data)
        return self._handle_response(response)


class EncodeClient(APIClient):
    """Client for the ENCODE API."""
    
    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize the ENCODE API client.
        
        Args:
            api_key (str, optional): API key for authentication
        """
        super().__init__("https://www.encodeproject.org/", api_key)
    
    def search_experiments(
        self, 
        assay_type: Optional[str] = None,
        cell_type: Optional[str] = None,
        status: str = "released",
        limit: int = 10
    ) -> List[Dict]:
        """
        Search for experiments in ENCODE.
        
        Args:
            assay_type (str, optional): Type of assay (e.g., 'ChIP-seq')
            cell_type (str, optional): Cell type
            status (str): Status of the experiment (default: 'released')
            limit (int): Maximum number of results to return
            
        Returns:
            List[Dict]: List of experiment data
        """
        params = {
            "type": "Experiment",
            "status": status,
            "limit": limit,
            "format": "json"
        }
        
        if assay_type:
            params["assay_title"] = assay_type
        
        if cell_type:
            params["biosample_ontology.term_name"] = cell_type
        
        response = self.get("search/", params)
        
        # Extract just the experiment data
        if "@graph" in response:
            return response["@graph"]
        return []
    
    def get_experiment(self, experiment_id: str) -> Dict:
        """
        Get details of a specific experiment.
        
        Args:
            experiment_id (str): ENCODE experiment ID
            
        Returns:
            Dict: Experiment data
        """
        return self.get(f"experiments/{experiment_id}/?format=json")
    
    def get_files(
        self, 
        experiment_id: str, 
        file_type: Optional[str] = None,
        output_type: Optional[str] = None
    ) -> List[Dict]:
        """
        Get files associated with an experiment.
        
        Args:
            experiment_id (str): ENCODE experiment ID
            file_type (str, optional): Type of file (e.g., 'bam', 'bigWig')
            output_type (str, optional): Output type (e.g., 'alignments', 'signal p-value')
            
        Returns:
            List[Dict]: List of file data
        """
        params = {
            "type": "File",
            "dataset": experiment_id,
            "format": "json"
        }
        
        if file_type:
            params["file_type"] = file_type
            
        if output_type:
            params["output_type"] = output_type
        
        response = self.get("search/", params)
        
        if "@graph" in response:
            return response["@graph"]
        return []
    
    def download_file(self, file_id: str, output_dir: str) -> str:
        """
        Download a file from ENCODE.
        
        Args:
            file_id (str): ENCODE file ID
            output_dir (str): Directory to save the file
            
        Returns:
            str: Path to the downloaded file
        """
        file_info = self.get(f"files/{file_id}/?format=json")
        if "href" not in file_info:
            raise ValueError(f"File {file_id} does not have a download URL")
        
        file_url = urljoin(self.base_url, file_info["href"])
        file_name = file_info.get("title", file_id)
        
        output_path = os.path.join(output_dir, file_name)
        os.makedirs(output_dir, exist_ok=True)
        
        # Download the file
        response = self.session.get(file_url, stream=True)
        response.raise_for_status()
        
        with open(output_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        return output_path


class GWASCatalogClient(APIClient):
    """Client for the GWAS Catalog API."""
    
    def __init__(self):
        """Initialize the GWAS Catalog API client."""
        super().__init__("https://www.ebi.ac.uk/gwas/rest/api/")
    
    def search_studies(
        self, 
        trait: Optional[str] = None,
        pubmed_id: Optional[str] = None,
        limit: int = 20
    ) -> List[Dict]:
        """
        Search for studies in the GWAS Catalog.
        
        Args:
            trait (str, optional): Disease/trait to search for
            pubmed_id (str, optional): PubMed ID of the publication
            limit (int): Maximum number of results to return
            
        Returns:
            List[Dict]: List of study data
        """
        params = {}
        
        if trait:
            # Use the studies endpoint with a trait filter
            endpoint = f"studies/search/findByDiseaseTrait?trait={trait}"
        elif pubmed_id:
            # Search by PubMed ID
            endpoint = f"studies/search/findByPublicationIdPubmedId?pubmedId={pubmed_id}"
        else:
            # Get all studies
            endpoint = "studies"
            params["size"] = limit
        
        response = self.get(endpoint, params)
        
        # Extract the study data
        if "_embedded" in response and "studies" in response["_embedded"]:
            return response["_embedded"]["studies"][:limit]
        return []
    
    def get_study(self, study_id: str) -> Dict:
        """
        Get details of a specific study.
        
        Args:
            study_id (str): GWAS study accession ID
            
        Returns:
            Dict: Study data
        """
        return self.get(f"studies/{study_id}")
    
    def get_associations(self, study_id: str) -> List[Dict]:
        """
        Get associations from a study.
        
        Args:
            study_id (str): GWAS study accession ID
            
        Returns:
            List[Dict]: List of association data
        """
        response = self.get(f"studies/{study_id}/associations")
        
        if "_embedded" in response and "associations" in response["_embedded"]:
            return response["_embedded"]["associations"]
        return []
    
    def search_variants(self, rsid: str) -> List[Dict]:
        """
        Search for variants by rsID.
        
        Args:
            rsid (str): rsID of the variant (e.g., 'rs12345')
            
        Returns:
            List[Dict]: List of variant data
        """
        response = self.get(f"singleNucleotidePolymorphisms/{rsid}/associations")
        
        if "_embedded" in response and "associations" in response["_embedded"]:
            return response["_embedded"]["associations"]
        return []
    
    def search_genes(self, gene_name: str) -> List[Dict]:
        """
        Search for associations by gene name.
        
        Args:
            gene_name (str): Gene name/symbol
            
        Returns:
            List[Dict]: List of association data
        """
        response = self.get(f"genes/{gene_name}/associations")
        
        if "_embedded" in response and "associations" in response["_embedded"]:
            return response["_embedded"]["associations"]
        return []


class ThousandGenomesClient(APIClient):
    """Client for the 1000 Genomes Project API."""
    
    def __init__(self):
        """Initialize the 1000 Genomes API client."""
        super().__init__("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/")
        # Add IGSR API for metadata
        self.igsr_api = "https://www.internationalgenome.org/api/beta/"
    
    def get_populations(self) -> pd.DataFrame:
        """
        Get all populations in the 1000 Genomes Project.
        
        Returns:
            pd.DataFrame: DataFrame with population data
        """
        response = requests.get(f"{self.igsr_api}population/all")
        response.raise_for_status()
        data = response.json()
        
        # Convert to DataFrame
        df = pd.DataFrame(data["_embedded"]["populations"])
        return df
    
    def get_samples(
        self, 
        population: Optional[str] = None,
        superpopulation: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Get samples from the 1000 Genomes Project.
        
        Args:
            population (str, optional): Population code (e.g., 'GBR')
            superpopulation (str, optional): Superpopulation code (e.g., 'EUR')
            
        Returns:
            pd.DataFrame: DataFrame with sample data
        """
        endpoint = "sample/all"
        params = {}
        
        if population:
            params["population"] = population
        if superpopulation:
            params["superpopulation"] = superpopulation
        
        response = requests.get(f"{self.igsr_api}{endpoint}", params=params)
        response.raise_for_status()
        data = response.json()
        
        # Convert to DataFrame
        df = pd.DataFrame(data["_embedded"]["samples"])
        return df
    
    def get_variant_data(
        self,
        chromosome: str,
        start: int,
        end: int,
        phase: int = 3,
        cache_dir: Optional[str] = None
    ) -> str:
        """
        Get variant data from a specific region.
        
        Args:
            chromosome (str): Chromosome name (e.g., '1', 'X')
            start (int): Start position
            end (int): End position
            phase (int): Phase of the 1000 Genomes Project (default: 3)
            cache_dir (str, optional): Directory to cache downloaded files
            
        Returns:
            str: Path to the VCF file
        """
        # Use the configured data path or a default
        if cache_dir is None:
            cache_dir = config.get_data_path("1000genomes")
        
        os.makedirs(cache_dir, exist_ok=True)
        
        # Determine the URL based on phase
        if phase == 3:
            base_path = "phase3/data"
        else:
            raise ValueError(f"Unsupported phase: {phase}")
        
        # For large regions, download the whole chromosome VCF
        if end - start > 10_000_000:
            vcf_url = f"{self.base_url}{base_path}/ALL.chr{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
            output_path = os.path.join(cache_dir, os.path.basename(vcf_url))
            
            if not os.path.exists(output_path):
                self.logger.info(f"Downloading VCF for chromosome {chromosome}...")
                response = requests.get(vcf_url, stream=True)
                response.raise_for_status()
                
                with open(output_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
            
            return output_path
        else:
            # For smaller regions, we could use the Tabix API, but as a fallback, download the chromosome VCF
            vcf_url = f"{self.base_url}{base_path}/ALL.chr{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
            output_path = os.path.join(cache_dir, os.path.basename(vcf_url))
            
            if not os.path.exists(output_path):
                self.logger.info(f"Downloading VCF for chromosome {chromosome}...")
                response = requests.get(vcf_url, stream=True)
                response.raise_for_status()
                
                with open(output_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
            
            return output_path


class OpenTargetsClient(APIClient):
    """Client for the Open Targets Platform API."""
    
    def __init__(self):
        """Initialize the Open Targets API client."""
        super().__init__("https://api.platform.opentargets.org/api/v4/graphql")
    
    def _execute_query(self, query: str, variables: Optional[Dict] = None) -> Dict:
        """
        Execute a GraphQL query.
        
        Args:
            query (str): GraphQL query
            variables (Dict, optional): Query variables
            
        Returns:
            Dict: Query result
        """
        payload = {"query": query}
        if variables:
            payload["variables"] = variables
        
        response = self.session.post(self.base_url, json=payload)
        response.raise_for_status()
        return response.json()
    
    def search_disease(self, disease_term: str, limit: int = 10) -> List[Dict]:
        """
        Search for diseases.
        
        Args:
            disease_term (str): Disease search term
            limit (int): Maximum number of results to return
            
        Returns:
            List[Dict]: List of disease data
        """
        query = """
        query SearchDisease($searchTerm: String!, $limit: Int!) {
          search(queryString: $searchTerm, entityNames: [DISEASE], page: {size: $limit}) {
            total
            aggregations {
              entities {
                name
                count
              }
            }
            hits {
              id
              entity {
                id
                name
                ... on Disease {
                  description
                }
              }
            }
          }
        }
        """
        
        variables = {
            "searchTerm": disease_term,
            "limit": limit
        }
        
        result = self._execute_query(query, variables)
        
        if "data" in result and "search" in result["data"] and "hits" in result["data"]["search"]:
            return result["data"]["search"]["hits"]
        return []
    
    def get_disease_associations(self, disease_id: str, limit: int = 100) -> List[Dict]:
        """
        Get target associations for a disease.
        
        Args:
            disease_id (str): Disease ID (e.g., 'EFO_0000253')
            limit (int): Maximum number of results to return
            
        Returns:
            List[Dict]: List of target association data
        """
        query = """
        query DiseaseAssociations($diseaseId: String!, $limit: Int!) {
          disease(efoId: $diseaseId) {
            id
            name
            associatedTargets(page: {size: $limit}) {
              count
              rows {
                target {
                  id
                  approvedSymbol
                  approvedName
                }
                score
                datatypeScores {
                  id
                  score
                }
              }
            }
          }
        }
        """
        
        variables = {
            "diseaseId": disease_id,
            "limit": limit
        }
        
        result = self._execute_query(query, variables)
        
        if "data" in result and "disease" in result["data"] and "associatedTargets" in result["data"]["disease"]:
            return result["data"]["disease"]["associatedTargets"]["rows"]
        return []
    
    def get_target_details(self, target_id: str) -> Dict:
        """
        Get details for a target.
        
        Args:
            target_id (str): Target ID (e.g., 'ENSG00000197386')
            
        Returns:
            Dict: Target data
        """
        query = """
        query TargetDetails($targetId: String!) {
          target(ensemblId: $targetId) {
            id
            approvedSymbol
            approvedName
            biotype
            genomicLocation {
              chromosome
              start
              end
              strand
            }
            subcellularLocations {
              location
              source
            }
            proteinAnnotations {
              id
              accession
              functions {
                term
                code
              }
            }
            hallmarks {
              attributes {
                name
                description
              }
              cancer {
                name
                promote
                suppress
                description
              }
              biologicalProcesses {
                name
                description
              }
            }
          }
        }
        """
        
        variables = {
            "targetId": target_id
        }
        
        result = self._execute_query(query, variables)
        
        if "data" in result and "target" in result["data"]:
            return result["data"]["target"]
        return {}
    
    def get_evidence(self, disease_id: str, target_id: str, limit: int = 20) -> List[Dict]:
        """
        Get evidence for a disease-target pair.
        
        Args:
            disease_id (str): Disease ID (e.g., 'EFO_0000253')
            target_id (str): Target ID (e.g., 'ENSG00000197386')
            limit (int): Maximum number of results to return
            
        Returns:
            List[Dict]: List of evidence data
        """
        query = """
        query EvidenceQuery($targetId: String!, $diseaseId: String!, $limit: Int!) {
          evidence(
            ensemblIds: [$targetId]
            efoIds: [$diseaseId]
            size: $limit
          ) {
            count
            rows {
              disease {
                id
                name
              }
              target {
                id
                approvedSymbol
              }
              score
              datasource
              datatypeId
              literature
            }
          }
        }
        """
        
        variables = {
            "targetId": target_id,
            "diseaseId": disease_id,
            "limit": limit
        }
        
        result = self._execute_query(query, variables)
        
        if "data" in result and "evidence" in result["data"] and "rows" in result["data"]["evidence"]:
            return result["data"]["evidence"]["rows"]
        return []