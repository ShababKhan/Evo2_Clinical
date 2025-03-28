#!/usr/bin/env python3
"""
Unit tests for the visualization module of the Evo2 pipeline.
"""

import unittest
import pandas as pd
import numpy as np
from pathlib import Path
import sys
import os
import matplotlib.pyplot as plt

# Add parent directory to path to import modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from visualization import VariantVisualizer

class TestVariantVisualizer(unittest.TestCase):
    """Test cases for the VariantVisualizer class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.visualizer = VariantVisualizer()
        
        # Create test data
        np.random.seed(42)  # For reproducibility
        self.test_variants = pd.DataFrame({
            'chrom': ['chr1'] * 50 + ['chr2'] * 50,
            'pos': np.random.randint(1, 1000000, 100),
            'gene': np.random.choice(['GATA2', 'EPAS1', 'KDR'], 100),
            'impact_score': np.random.random(100),
            'variant_type': np.random.choice(['SNP', 'INDEL', 'SV'], 100),
            'DNase': np.random.choice([0, 1], 100),
            'H3K27ac': np.random.choice([0, 1], 100),
            'H3K4me3': np.random.choice([0, 1], 100),
            'EUR': np.random.random(100),
            'EAS': np.random.random(100),
            'AFR': np.random.random(100)
        })
        
        # Create test output directory
        self.test_output_dir = Path("test_output")
        self.test_output_dir.mkdir(exist_ok=True)
        
    def tearDown(self):
        """Clean up after tests."""
        # Remove test output files
        for file in self.test_output_dir.glob("*"):
            file.unlink()
        self.test_output_dir.rmdir()
        
    def test_plot_variant_distribution(self):
        """Test variant distribution plotting."""
        output_path = self.test_output_dir / "variant_dist.png"
        self.visualizer.plot_variant_distribution(
            self.test_variants,
            save_path=str(output_path)
        )
        self.assertTrue(output_path.exists())
        
    def test_plot_variant_distribution_by_chromosome(self):
        """Test variant distribution plotting for a specific chromosome."""
        output_path = self.test_output_dir / "chr1_dist.png"
        self.visualizer.plot_variant_distribution(
            self.test_variants,
            chromosome="chr1",
            save_path=str(output_path)
        )
        self.assertTrue(output_path.exists())
        
    def test_plot_impact_heatmap(self):
        """Test impact heatmap plotting."""
        output_path = self.test_output_dir / "impact_heatmap.png"
        self.visualizer.plot_impact_heatmap(
            self.test_variants,
            genes=['GATA2', 'EPAS1', 'KDR'],
            save_path=str(output_path)
        )
        self.assertTrue(output_path.exists())
        
    def test_plot_population_frequencies(self):
        """Test population frequency plotting."""
        output_path = self.test_output_dir / "pop_freq.png"
        self.visualizer.plot_population_frequencies(
            self.test_variants,
            populations=['EUR', 'EAS', 'AFR'],
            save_path=str(output_path)
        )
        self.assertTrue(output_path.exists())
        
    def test_plot_encode_integration(self):
        """Test ENCODE integration plotting."""
        output_path = self.test_output_dir / "encode_features.png"
        self.visualizer.plot_encode_integration(
            self.test_variants,
            encode_features=['DNase', 'H3K27ac', 'H3K4me3'],
            save_path=str(output_path)
        )
        self.assertTrue(output_path.exists())
        
    def test_plot_lncrna_analysis(self):
        """Test lncRNA analysis plotting."""
        output_path = self.test_output_dir / "lncrna_analysis.html"
        self.visualizer.plot_lncrna_analysis(
            self.test_variants,
            lncrna_id="GATA2-AS1",
            save_path=str(output_path)
        )
        self.assertTrue(output_path.exists())
        
    def test_create_report_figures(self):
        """Test creation of all report figures."""
        self.visualizer.create_report_figures(
            self.test_variants,
            str(self.test_output_dir)
        )
        
        expected_files = [
            "variant_distribution.png",
            "impact_heatmap.png",
            "encode_integration.png",
            "gata2as1_analysis.html"
        ]
        
        for file in expected_files:
            self.assertTrue((self.test_output_dir / file).exists())


if __name__ == '__main__':
    unittest.main()