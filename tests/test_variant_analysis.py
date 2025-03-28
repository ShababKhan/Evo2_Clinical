#!/usr/bin/env python3
"""
Unit tests for the variant analysis module of the Evo2 pipeline.
"""

import unittest
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import os
from variant_analysis import VariantProcessor

class TestVariantProcessor(unittest.TestCase):
    """Test cases for the VariantProcessor class."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        cls.processor = VariantProcessor()
        
        # Create test data
        cls.test_variants = pd.DataFrame({
            'chrom': ['chr1', 'chr1', 'chr2'],
            'pos': [1000, 2000, 3000],
            'id': ['rs1', 'rs2', 'rs3'],
            'ref': ['A', 'C', 'G'],
            'alt': ['T', 'G', 'T'],
            'qual': [100, 200, 300],
            'filter': ['PASS', 'PASS', 'FAIL'],
            'gene': ['GATA2', 'EPAS1', 'GATA2-AS1'],
            'impact_score': [0.8, 0.5, 0.9],
            'is_endothelial': [True, True, False]
        })
        
    def test_initialization(self):
        """Test VariantProcessor initialization."""
        self.assertEqual(self.processor.reference_genome, "GRCh38")
        self.assertEqual(len(self.processor.variant_cache), 0)
        
    def test_score_variant_impact(self):
        """Test variant impact scoring."""
        variant_info = {
            'ref': 'A',
            'alt': 'T',
            'gene': 'GATA2'
        }
        score = self.processor.score_variant_impact(variant_info)
        self.assertIsInstance(score, float)
        self.assertTrue(0 <= score <= 1)
        
    def test_analyze_endothelial_variants(self):
        """Test endothelial variant analysis."""
        endothelial_genes = ['GATA2', 'EPAS1']
        result = self.processor.analyze_endothelial_variants(
            self.test_variants,
            endothelial_genes
        )
        
        self.assertEqual(len(result), 2)
        self.assertTrue(all(result['gene'].isin(endothelial_genes)))
        
    def test_analyze_lncrna_variants(self):
        """Test lncRNA variant analysis."""
        result = self.processor.analyze_lncrna_variants(
            self.test_variants,
            lncrna_id='GATA2-AS1'
        )
        
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]['gene'], 'GATA2-AS1')
        
    def test_annotate_variants(self):
        """Test variant annotation."""
        encode_data = pd.DataFrame({
            'chrom': ['chr1', 'chr2'],
            'pos': [1000, 3000],
            'DNase': [1, 0],
            'H3K27ac': [1, 1]
        })
        
        result = self.processor.annotate_variants(
            self.test_variants,
            encode_data
        )
        
        self.assertTrue('DNase' in result.columns)
        self.assertTrue('H3K27ac' in result.columns)
        
    def test_generate_variant_report(self):
        """Test report generation."""
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            self.processor.generate_variant_report(self.test_variants, tmp.name)
            self.assertTrue(os.path.exists(tmp.name))
            
            with open(tmp.name, 'r') as f:
                content = f.read()
                self.assertIn('total_variants', content)
                self.assertIn('high_impact_variants', content)
                self.assertIn('endothelial_variants', content)
            
            os.unlink(tmp.name)  # Clean up
            
    def test_analyze_population_frequencies(self):
        """Test population frequency analysis."""
        result = self.processor.analyze_population_frequencies(
            self.test_variants,
            population="ALL"
        )
        self.assertEqual(len(result), len(self.test_variants))


if __name__ == '__main__':
    unittest.main()