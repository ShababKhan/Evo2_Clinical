# 🧬 Evo2 Clinical

A powerful genomic analysis pipeline for studying genetic variants across genes and diseases.

## Overview

Evo2 Clinical is a comprehensive pipeline for analyzing genetic variants and their functional impacts. It integrates multiple data sources and computational methods to provide insights into variant effects across any genes or genomic regions of interest.

## Core Technology

The Evo2 Clinical pipeline leverages several key technologies:

1. **Evo2 Model**: A powerful AI system for analyzing genomic regions and predicting functional impacts of variants.

2. **Genomic Data Integration**:
   - 1000 Genomes Project data for population genetics
   - GWAS Catalog for disease associations
   - ENCODE data for regulatory elements across cell types

3. **Bioinformatics Libraries**:
   - CyVCF2 and Pysam for efficient VCF processing
   - Biopython for sequence analysis
   - PyBedTools for genomic interval operations

4. **Visualization Framework**:
   - Interactive visualizations with Plotly
   - Statistical plots with Matplotlib/Seaborn

5. **AIDO Simulations**:
   - Predictive modeling of phenotypic effects

## Features

- **Flexible Gene Selection**: Analyze any gene or gene set of interest
- **Variant Analysis Pipeline**: Process and annotate genetic variants from VCF files
- **Impact Scoring**: Predict functional impacts of variants using Evo2 and other methods
- **Multi-Disease Support**: Configure analyses for any disease or trait of interest
- **Population Genetics**: Compare variant frequencies across different populations
- **Cell-Type Specificity**: Filter analyses for specific cell types using ENCODE data
- **Interactive Visualizations**: Generate comprehensive reports with interactive plots
- **GWAS Integration**: Connect variants to disease associations from GWAS studies
- **Configurable Workflows**: Easily adapt pipelines for different research questions

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/Evo2_Clinical.git
cd Evo2_Clinical

# Install dependencies
pip install -r requirements.txt
```

## Usage / Guides

### 1. Command Line Interface

```bash
# Analyze variants in specific genes
python evo2_pipeline.py --data-dir data --output-dir results --genes BRCA1 BRCA2 TP53

# Analyze variants with cell type and trait filters
python evo2_pipeline.py --genes BRCA1 BRCA2 --cell-types BREAST MCF7 --traits "breast cancer"

# Analyze a single gene with region specification
python examples/analyze_gene_variants.py --gene-id BRCA1 --region chr17:41196312-41277500
```

### 2. Python API

```python
from evo2_pipeline import Evo2Pipeline

# Initialize the pipeline with custom configuration
config = {
    "data_dir": "data",
    "output_dir": "results/my_analysis"
}
pipeline = Evo2Pipeline(config)

# Define analysis parameters
genes = ["BRCA1", "BRCA2", "TP53"]
cell_types = ["BREAST", "MCF7"]
traits = ["breast cancer", "ovarian cancer"]

# Run analysis
pipeline.run_pipeline(
    genes=genes,
    cell_types=cell_types,
    traits=traits
)
```

### 3. Example Scripts

The `examples/` directory contains several scripts demonstrating common analysis workflows:

- `analyze_gene_variants.py`: Analyze variants in a single gene
- `analyze_gene_set.py`: Analyze variants across multiple related genes
- `analyze_population_frequencies.py`: Compare variant frequencies across populations

## Output

The pipeline generates comprehensive analysis results in the specified output directory:

```
results/
├── gene_analysis/
│   ├── GENE1_analysis.csv
│   ├── GENE1_distribution.png
│   └── GENE1_encode_features.png
├── population_frequencies/
│   ├── AFR_frequencies.csv
│   ├── EUR_frequencies.csv
│   └── population_frequencies.png
└── report/
    ├── variant_distribution.png
    └── impact_scores.png
```

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use Evo2 Clinical in your research, please cite:

```bibtex
@software{evo2_clinical_2025,
    title = {Evo2 Clinical: A Comprehensive Genomic Analysis Pipeline},
    author = {Khan, Shabab},
    year = {2025},
    url = {https://github.com/yourusername/Evo2_Clinical}
}
```

## Interactive Web Interface

The project includes an interactive web interface built with Streamlit and PyGWalker for exploring genomic variants:

```bash
# Install dependencies
pip install -r requirements.txt

# Run the Streamlit app
streamlit run app.py
```

### Features of the Web Interface

1. **Data Upload**
   - Support for VCF files (.vcf, .vcf.gz)
   - Automatic variant processing and annotation

2. **Analysis Configuration**
   - Gene selection (multiple genes supported)
   - Cell type filtering for ENCODE data
   - Disease/trait association analysis

3. **Interactive Visualization**
   - PyGWalker interface for custom analysis
   - Preset visualizations:
     - Variant distribution plots
     - Impact score analysis
     - Population frequency comparisons
     - ENCODE feature integration

4. **Data Export**
   - Download processed variant data
   - Export visualization results

### System Requirements

- Python 3.8 or higher
- Modern web browser
- Minimum 4GB RAM recommended
