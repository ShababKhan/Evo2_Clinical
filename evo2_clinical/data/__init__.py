"""
Data loading and processing module for the Evo2_Clinical package.
"""

from .io import (
    load_gwas_catalog_data,
    load_vcf_file,
    load_1000_genomes_data,
    filter_by_encode_data,
    load_gene_list,
    load_lncrna_list,
    load_emt_genes,
    load_gtf_file
)