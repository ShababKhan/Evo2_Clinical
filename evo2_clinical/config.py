"""
Configuration module for the Evo2_Clinical package.

This module provides configuration settings for data paths, tool paths, and output directories.
"""

import os
import yaml
from pathlib import Path


class Config:
    """Configuration class for the Evo2_Clinical package."""
    
    def __init__(self, config_file=None):
        """
        Initialize configuration with default values or from a config file.
        
        Args:
            config_file (str, optional): Path to a YAML configuration file.
        """
        # Default configuration
        self.DATA_PATHS = {
            "gwas_catalog": "",
            "1000_genomes_vcf": "",
            "hapmap_data": "",
            "encode_data": "",
            "endothelial_genes_list": "",
            "lncrna_list": "",
            "emt_genes_list": "",
            "cteph_gwas_loci": "",
            "pah_known_variants": "",
            "mesothelioma_variants": "",
        }
        
        self.TOOL_PATHS = {
            "evo2_executable": "",
            "aido_executable": "",
        }
        
        self.OUTPUT_DIRS = {
            "databases": "output/databases",
            "intermediate_files": "output/intermediate",
            "logs": "output/logs",
        }
        
        # Load config from file if provided
        if config_file:
            self.load_config(config_file)
        
    def load_config(self, config_file):
        """
        Load configuration from a YAML file.
        
        Args:
            config_file (str): Path to a YAML configuration file.
        
        Returns:
            bool: True if successful, False otherwise.
        """
        try:
            with open(config_file, 'r') as f:
                config_data = yaml.safe_load(f)
                
            # Update configurations if present in the file
            if 'DATA_PATHS' in config_data:
                self.DATA_PATHS.update(config_data['DATA_PATHS'])
            if 'TOOL_PATHS' in config_data:
                self.TOOL_PATHS.update(config_data['TOOL_PATHS'])
            if 'OUTPUT_DIRS' in config_data:
                self.OUTPUT_DIRS.update(config_data['OUTPUT_DIRS'])
                
            return True
        except Exception as e:
            print(f"Error loading configuration: {e}")
            return False
    
    def save_config(self, config_file):
        """
        Save the current configuration to a YAML file.
        
        Args:
            config_file (str): Path to save the YAML configuration file.
            
        Returns:
            bool: True if successful, False otherwise.
        """
        try:
            config_data = {
                'DATA_PATHS': self.DATA_PATHS,
                'TOOL_PATHS': self.TOOL_PATHS,
                'OUTPUT_DIRS': self.OUTPUT_DIRS
            }
            
            # Create directory if it doesn't exist
            os.makedirs(os.path.dirname(config_file), exist_ok=True)
            
            with open(config_file, 'w') as f:
                yaml.dump(config_data, f, default_flow_style=False)
                
            return True
        except Exception as e:
            print(f"Error saving configuration: {e}")
            return False
    
    def get_absolute_path(self, path):
        """
        Convert a relative path to absolute path based on package location.
        
        Args:
            path (str): Relative path
            
        Returns:
            str: Absolute path
        """
        if os.path.isabs(path):
            return path
        
        # Get the package root directory
        package_root = Path(__file__).parent.parent
        return os.path.join(package_root, path)
    
    def ensure_output_dirs(self):
        """
        Create output directories if they don't exist.
        
        Returns:
            bool: True if all directories were created or already exist, False otherwise.
        """
        try:
            for dir_path in self.OUTPUT_DIRS.values():
                abs_path = self.get_absolute_path(dir_path)
                os.makedirs(abs_path, exist_ok=True)
            return True
        except Exception as e:
            print(f"Error creating output directories: {e}")
            return False


# Create a default configuration instance
config = Config()