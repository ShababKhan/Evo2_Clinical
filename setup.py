# install all data science and bioinformatics libraries for endothelial research program
import os 

os.system('pip install numpy pandas matplotlib seaborn scikit-learn scipy' \
' statsmodels tensorflow keras nltk spacy gensim xgboost lightgbm catboost' \
' opencv-python dlib pytesseract' \
          # Bioinformatics libraries
          'biopython pysam pyvcf pybedtools hgvs cyvcf2' \
          # Genomics analysis
          'plink scikit-allel PyVCF3 gtfparse pyensembl' \
          # Data processing and visualization for genomics
          'plotly dash anndata scanpy leidenalg' \
          # Machine learning extensions
          'shap lime eli5' \
          # Statistical genetics
          'statsmodels pandas-plink' \
          # API clients for genomic databases
          'biomart intermine gget')


# Core data science
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn
import scipy
import statsmodels
import tensorflow as tf
import keras

# NLP and text processing
import nltk
import spacy
import gensim

# Advanced ML
import xgboost
import lightgbm
import catboost
import shap  # For model interpretability

# Computer vision
import cv2
import dlib
import pytesseract

# Bioinformatics and genomics
try:
    import Bio  # Biopython for sequence analysis
    import pysam  # For SAM/BAM/VCF files
    import vcf  # PyVCF for VCF files
    import pybedtools  # For BED file manipulation
    import hgvs  # For genomic variant notation
    import cyvcf2  # Fast VCF parsing
    import allel  # Scikit-allel for genomic data structures
    from gtfparse import read_gtf  # Parse GTF files
    import pyensembl  # Access Ensembl data
    import plotly.express as px  # Advanced visualization
    
    print("All required libraries for the Endothelial Research Program successfully imported.")
    print("Ready for GWAS analysis, variant scoring with Evo2, and genomic data processing.")
except ImportError as e:
    print(f"Warning: Some bioinformatics libraries not properly installed: {e}")
    print("Run this script again to attempt reinstallation.") 