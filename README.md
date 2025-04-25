# Evo2_Clinical: Decoding Endothelial Influence on Pulmonary Disease

![Version](https://img.shields.io/badge/version-0.1.0-blue)
![Python](https://img.shields.io/badge/python-%3E%3D3.8-green)
![License](https://img.shields.io/badge/license-MIT-yellow)

**Evo2_Clinical** is a computational framework for investigating endothelial influence on pulmonary disease, integrating computational insights with molecular biology. This project aims to understand the role of the endothelium in pulmonary disease progression throughout cancer treatment.

## 🧬 Project Overview

Evo2_Clinical focuses on the endothelium (endothelial cells) and their role in lung inflammation and fibrosis during cancer care. The framework integrates multiple data sources, scoring algorithms, and analytical methods to provide insights into genetic variants and lncRNA functionality in endothelial cells.

### 🔍 Disease Areas

- Chronic Thromboembolic Pulmonary Hypertension (CTEPH)
- Pulmonary Arterial Hypertension (PAH)
- Radiation/Drug-induced Pneumonitis and Fibrosis
- Malignant Pleural Mesothelioma
- Pulmonary Vascular Diseases

### 🧪 Specific Molecules/Genes of Interest

- **lncRNAs**: ANRIL, PIRAT, LRAC (LINC01899), GATA2-AS1
- Genes involved in Hypoxia Regulation (e.g., EPAS1 SNPs)
- Genes involved in Endothelial-Mesenchymal Transition (EMT)
- Epigenetic Mediators
- Endothelial genes

## ⚙️ Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/Evo2_Clinical.git
cd Evo2_Clinical

# Install the package
pip install -e .

# Install dependencies
pip install -r requirements.txt
```

## 🚀 Quick Start

### Using the Command Line Interface (CLI)

```bash
# Run the basic pipeline
evo2-clinical run --config path/to/config.yml

# Score variants
evo2-clinical score-variants --vcf path/to/variants.vcf --output results/

# Analyze lncRNAs
evo2-clinical analyze-lncrna --input path/to/lncrna_list.txt --output results/
```

### Using the Streamlit GUI

```bash
# Launch the Streamlit interface
streamlit run evo2_clinical_app.py
```

### Using the Python API

```python
from evo2_clinical.pipeline.main import Pipeline
from evo2_clinical.config import Config

# Initialize the pipeline with config
config = Config("path/to/config.yml")
pipeline = Pipeline(config)

# Run the full pipeline
pipeline.run()

# Or run specific analyses
variants_df = pipeline.score_variants(vcf_path="path/to/variants.vcf")
lncrna_df = pipeline.analyze_lncrnas(lncrna_list=["GATA2-AS1", "ANRIL"])
```

## 📊 Database Outputs

The framework generates several databases containing critical information:

1. **Endothelial Variants Database**: Common and rare genetic variants from the 1000 Genomes project found within endothelial genes.
2. **lncRNA Functionality Database**: Functionality scores for endothelial lncRNAs.
3. **GATA2-AS1 Predictions Database**: Predictions of variant effects specifically for the lncRNA GATA2-AS1.
4. **Epigenetic Mediator Scores Database**: Scores for epigenetic mediators with functionally relevant SNPs.
5. **Functional Variants Database**: Focused on EMT pathway genes, CTEPH GWAS, PAH variants, and mesothelioma variants.

## 🧠 Core Components

### Data Integration

- GWAS Catalogs integration (trait-associated loci)
- HapMap/1000 Genomes Data processing
- ENCODE Data integration for cell-specific activity filtering (ENDOS)

### Analysis Tools

- Variant & Function Scoring using the Evo2 methodology
- lncRNA functionality scoring
- Phenotypic simulation
- Translational bridging between animal models and human contexts

## 📁 Project Structure

```
Evo2_Clinical/
├── evo2_clinical/            # Core package
│   ├── __init__.py
│   ├── cli.py                # Command-line interface
│   ├── config.py             # Configuration utilities
│   ├── analysis/             # Analysis modules
│   │   ├── variant_scoring.py
│   │   └── lncrna_scoring.py
│   ├── data/                 # Data handling modules
│   │   └── io.py
│   ├── database/             # Database management
│   │   └── manager.py
│   ├── pipeline/             # Pipeline orchestration
│   │   └── main.py
│   └── utils/                # Utility functions
│       └── helpers.py
├── evo2_clinical_app.py      # Streamlit GUI application
├── data/                     # Example/reference data
├── examples/                 # Example scripts
├── tests/                    # Unit tests
├── output/                   # Output directory
│   ├── databases/            # Generated databases
│   └── logs/                 # Log files
├── requirements.txt          # Python dependencies
├── setup.py                  # Package setup file
└── README.md                 # This file
```

## 🔧 Configuration

The Evo2_Clinical framework uses YAML configuration files to specify data paths, tool settings, and output directories. Example:

```yaml
DATA_PATHS:
  gwas_catalog: "data/gwas/gwas-catalog-associations.tsv"
  1000_genomes_vcf: "data/1000genomes/chr19_variants.vcf.gz"
  reference_genome: "data/reference/GRCh38.chr19.fa.gz"
  
TOOL_PATHS:
  evo2_executable: "bin/evo2_tool"
  
OUTPUT_DIRS:
  databases: "output/databases"
  logs: "logs"
```

## 📚 Documentation

For detailed documentation, please refer to the [Wiki](https://github.com/yourusername/Evo2_Clinical/wiki) or the [Documentation](https://evo2-clinical.readthedocs.io/).

## 👥 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🙏 Acknowledgments
* Jeffrey Man, MD/PhD