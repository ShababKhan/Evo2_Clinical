# Evo2 Clinical

A comprehensive bioinformatics pipeline for analyzing endothelial genetic variants in pulmonary diseases and cancer care.

## Overview

Evo2 Clinical is a Python-based computational pipeline that integrates various genomic analysis tools to study the role of endothelial cells in pulmonary disease progression, with a special focus on cancer treatment contexts.

## Installation

```bash
# Using pipenv (recommended)
pipenv install

# Using pip
pip install -r requirements.txt
```

## Key Features

### 1. Variant Analysis
- Analysis of common and rare variants from 1000 Genomes Project
- Focus on endothelial-specific genes
- Integration with GWAS catalog data
- Identification of linkage disequilibrium regions

### 2. LncRNA Analysis
- Specialized analysis of GATA2-AS1 and other endothelial lncRNAs
- Functional scoring system for lncRNA variants
- Integration with ENCODE data for cell-specific activity

### 3. Population Genetics
- Analysis of population-specific variants
- Special focus on adaptive variants (e.g., EPAS1 mutations in Tibetan Highlanders)
- Comparative population frequency analysis

### 4. Visualization Tools
- Population frequency distribution plots
- Variant distribution visualization
- ENCODE feature visualization
- Custom report generation

## Usage Examples

### Analyzing Population Frequencies
```python
from variant_analysis import analyze_population_frequencies
from visualization import plot_frequencies

# Analyze frequencies across populations
frequencies = analyze_population_frequencies("chr19")
```

### LncRNA Analysis
```python
from variant_analysis import analyze_lncrna_variants

# Analyze GATA2-AS1 variants
results = analyze_lncrna_variants("GATA2-AS1")
```

### Endothelial Variant Analysis
```python
from variant_analysis import analyze_endothelial_variants

# Analyze endothelial-specific variants
variants = analyze_endothelial_variants()
```

## Directory Structure

- `data/`: Raw and processed data files
- `examples/`: Example scripts demonstrating package usage
- `results/`: Output files and analysis results
- `tests/`: Unit tests
- Main pipeline modules:
  - `evo2_pipeline.py`: Core pipeline functionality
  - `variant_analysis.py`: Variant analysis functions
  - `visualization.py`: Data visualization tools

## Results Output

Analysis results are stored in the `results/` directory with the following structure:
- `chr19_analysis/`: Population-specific frequency analysis
- `lncrna_analysis/`: LncRNA variant analysis results
- Each analysis includes:
  - CSV files with detailed data
  - Visualization plots
  - HTML reports where applicable

## Testing

Run the test suite using:
```bash
pytest tests/
```

## Dependencies

Key dependencies include:
- numpy
- pandas
- scipy
- biopython
- pysam
- cyvcf2
- scikit-allel
- matplotlib
- seaborn
- plotly

## Contributing

Please read our contributing guidelines before submitting pull requests.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this software in your research, please cite:
[Citation information to be added]
