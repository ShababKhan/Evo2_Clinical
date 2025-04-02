#!/usr/bin/env python3
"""
Smart Database Download Script for Evo2 Pipeline

Downloads and organizes genomic databases needed for variant analysis:
- Gene annotations (RefSeq) for the entire genome
- Reference genome (GRCh38 complete genome or selected chromosomes)
- 1000 Genomes Project VCF files for all or selected chromosomes
- ENCODE data for various cell types
- GWAS Catalog

Smart downloading features:
- Checks if files already exist before downloading
- Only downloads missing files
- Supports full genome or selective chromosome-by-chromosome downloads
- Tracks download status and logs results
- Provides option to force re-download if needed

Author: Shabab Khan
Date: April 1, 2025
"""

import os
import sys
import logging
import gzip
import shutil
import urllib.request
from pathlib import Path
from tqdm import tqdm
import argparse
import concurrent.futures

class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)

class SmartDatabaseDownloader:
    def __init__(self, force_download=False, chromosomes=None, max_threads=4):
        """
        Initialize the downloader.
        
        Args:
            force_download: Whether to force download even if files exist
            chromosomes: List of chromosomes to download (e.g., ["1", "2", "X"])
                         or None for all chromosomes
            max_threads: Maximum number of parallel downloads
        """
        self.force_download = force_download
        self.max_threads = max_threads
        
        # Set chromosomes to download
        self.all_chromosomes = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
        self.chromosomes = chromosomes if chromosomes else self.all_chromosomes
        
        self.setup_logging()
        self.setup_directories()
        
    def setup_logging(self):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.StreamHandler(),
                logging.FileHandler('database_download.log')
            ]
        )
        self.logger = logging.getLogger("SmartDatabaseDownloader")

    def setup_directories(self):
        """Create necessary directories if they don't exist."""
        dirs = [
            'data/reference',
            'data/1000genomes',
            'data/encode',
            'data/gwas'
        ]
        
        # Add chromosome-specific directories
        for chrom in self.chromosomes:
            dirs.append(f'data/reference/chr{chrom}')
            dirs.append(f'data/1000genomes/chr{chrom}')
        
        for dir_path in dirs:
            Path(dir_path).mkdir(parents=True, exist_ok=True)

    def download_with_progress(self, url, output_path):
        """Download a file with progress bar if it doesn't exist or force_download is True."""
        # Check if file already exists
        if os.path.exists(output_path) and not self.force_download:
            self.logger.info(f"File already exists (skipping download): {output_path}")
            return True
            
        # File doesn't exist or force_download is True, so download it
        try:
            self.logger.info(f"Downloading {url} to {output_path}...")
            with DownloadProgressBar(unit='B', unit_scale=True,
                               miniters=1, desc=f"Downloading {os.path.basename(output_path)}") as t:
                urllib.request.urlretrieve(url, filename=output_path,
                                     reporthook=t.update_to)
            self.logger.info(f"Successfully downloaded: {output_path}")
            return True
        except Exception as e:
            self.logger.error(f"Error downloading {url}: {str(e)}")
            return False

    def decompress_file(self, gz_path):
        """Decompress a gzipped file if the decompressed version doesn't exist."""
        output_path = str(gz_path)[:-3]  # Remove .gz extension
        
        # Check if decompressed file already exists
        if os.path.exists(output_path) and not self.force_download:
            self.logger.info(f"Decompressed file already exists (skipping): {output_path}")
            return True
            
        try:
            self.logger.info(f"Decompressing {gz_path}...")
            with gzip.open(gz_path, 'rb') as f_in:
                with open(output_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            self.logger.info(f"Decompressed {gz_path} to {output_path}")
            return True
        except Exception as e:
            self.logger.error(f"Error decompressing {gz_path}: {str(e)}")
            return False

    def download_gene_annotations(self):
        """Download and process gene annotations for the entire genome."""
        self.logger.info("Processing gene annotations for the entire genome...")
        url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz"
        output_path = "data/reference/GRCh38.genes.gtf.gz"
        
        if self.download_with_progress(url, output_path):
            return output_path
        return None

    def download_reference_genome_full(self):
        """Download complete reference genome."""
        self.logger.info("Processing complete reference genome...")
        url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
        output_path = "data/reference/GRCh38.fa.gz"
        
        if self.download_with_progress(url, output_path):
            return output_path
        return None
        
    def download_reference_genome_chromosome(self, chromosome):
        """Download reference sequence for a specific chromosome."""
        self.logger.info(f"Processing reference genome for chromosome {chromosome}...")
        url = f"https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr{chromosome}.fa.gz"
        output_path = f"data/reference/chr{chromosome}/GRCh38.chr{chromosome}.fa.gz"
        
        if self.download_with_progress(url, output_path):
            return output_path
        return None

    def download_reference_genome(self):
        """Download reference genome - full or by chromosome."""
        all_files = []
        
        # If all chromosomes are requested, download the full genome
        if set(self.chromosomes) == set(self.all_chromosomes):
            if full_genome := self.download_reference_genome_full():
                all_files.append(full_genome)
        else:
            # Otherwise download individual chromosomes
            with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_threads) as executor:
                futures = []
                for chrom in self.chromosomes:
                    futures.append(
                        executor.submit(self.download_reference_genome_chromosome, chrom)
                    )
                
                for future in concurrent.futures.as_completed(futures):
                    if result := future.result():
                        all_files.append(result)
        
        return all_files

    def download_1000g_data_chromosome(self, chromosome):
        """Download 1000 Genomes Project data for a specific chromosome."""
        self.logger.info(f"Processing 1000 Genomes Project data for chromosome {chromosome}...")
        url = f"https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/ALL.chr{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
        output_path = f"data/1000genomes/chr{chromosome}/chr{chromosome}_variants.vcf.gz"
        
        if self.download_with_progress(url, output_path):
            return output_path
        return None

    def download_1000g_data(self):
        """Download 1000 Genomes Project data for all selected chromosomes."""
        all_files = []
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_threads) as executor:
            futures = []
            for chrom in self.chromosomes:
                futures.append(
                    executor.submit(self.download_1000g_data_chromosome, chrom)
                )
            
            for future in concurrent.futures.as_completed(futures):
                if result := future.result():
                    all_files.append(result)
        
        return all_files

    def download_encode_data(self):
        """Download ENCODE data for various cell types."""
        self.logger.info("Processing ENCODE data...")
        
        # Expanded collection of ENCODE datasets
        dataset_info = [
            # HUVEC (Human Umbilical Vein Endothelial Cells) datasets
            {
                "url": "https://www.encodeproject.org/files/ENCFF017XLW/@@download/ENCFF017XLW.bed.gz",
                "output": "data/encode/HUVEC_H3K27ac_ChIPseq.bed.gz",
                "description": "HUVEC H3K27ac ChIP-seq peaks"
            },
            {
                "url": "https://www.encodeproject.org/files/ENCFF721JMB/@@download/ENCFF721JMB.bed.gz",
                "output": "data/encode/HUVEC_DNase_seq.bed.gz",
                "description": "HUVEC DNase-seq peaks"
            },
            # HepG2 (Liver) datasets
            {
                "url": "https://www.encodeproject.org/files/ENCFF798RSS/@@download/ENCFF798RSS.bed.gz",
                "output": "data/encode/HepG2_H3K4me3_ChIPseq.bed.gz",
                "description": "HepG2 H3K4me3 ChIP-seq peaks"
            },
            # GM12878 (Lymphoblastoid) datasets
            {
                "url": "https://www.encodeproject.org/files/ENCFF417EGU/@@download/ENCFF417EGU.bed.gz",
                "output": "data/encode/GM12878_CTCF_ChIPseq.bed.gz",
                "description": "GM12878 CTCF ChIP-seq peaks"
            }
        ]
        
        outputs = []
        for dataset in dataset_info:
            self.logger.info(f"Downloading {dataset['description']}...")
            if self.download_with_progress(dataset['url'], dataset['output']):
                outputs.append(dataset['output'])
        
        return outputs

    def download_gwas_catalog(self):
        """Download GWAS Catalog."""
        self.logger.info("Processing GWAS Catalog...")
        url = "https://www.ebi.ac.uk/gwas/api/search/downloads/full"
        output_path = "data/gwas/gwas-catalog-associations.tsv"
        
        if self.download_with_progress(url, output_path):
            return output_path
        return None
        
    def download_dataset(self, dataset_name):
        """Download a specific dataset by name."""
        self.logger.info(f"Downloading dataset: {dataset_name}")
        
        if dataset_name == "gene_annotations":
            return self.download_gene_annotations()
        elif dataset_name == "reference_genome":
            return self.download_reference_genome()
        elif dataset_name == "1000g_data":
            return self.download_1000g_data()
        elif dataset_name == "encode_data":
            return self.download_encode_data()
        elif dataset_name == "gwas_catalog":
            return self.download_gwas_catalog()
        else:
            self.logger.error(f"Unknown dataset name: {dataset_name}")
            return None

    def run(self, datasets=None):
        """Execute the database download process for specified datasets or all if none specified."""
        try:
            files_to_decompress = []
            
            # If no datasets specified, download all
            if not datasets:
                datasets = ["gene_annotations", "reference_genome", "1000g_data", "encode_data", "gwas_catalog"]
                
            self.logger.info(f"Starting download process for datasets: {', '.join(datasets)}")
            self.logger.info(f"Selected chromosomes: {', '.join(self.chromosomes)}")
            
            # Download specified datasets
            for dataset in datasets:
                result = self.download_dataset(dataset)
                
                # Add to decompression list if needed
                if dataset == "gene_annotations" and result:
                    files_to_decompress.append(result)
                elif dataset == "reference_genome" and result:
                    files_to_decompress.extend(result if isinstance(result, list) else [result])
                elif dataset == "1000g_data" and result:
                    files_to_decompress.extend(result if isinstance(result, list) else [result])
                elif dataset == "encode_data" and result:
                    files_to_decompress.extend(result if isinstance(result, list) else [result])
            
            # Decompress downloaded files (if they need to be decompressed)
            if files_to_decompress:
                self.logger.info(f"Decompressing {len(files_to_decompress)} downloaded files...")
                for gz_file in files_to_decompress:
                    self.decompress_file(gz_file)
            
            self.logger.info("All requested database downloads processed successfully")
            
        except Exception as e:
            self.logger.error(f"An error occurred during download process: {str(e)}")
            sys.exit(1)

def main():
    """Parse command line arguments and run the downloader."""
    parser = argparse.ArgumentParser(description="Smart downloader for genomic databases")
    parser.add_argument("--datasets", nargs="+", choices=["gene_annotations", "reference_genome", 
                                                       "1000g_data", "encode_data", "gwas_catalog"],
                     help="Specific datasets to download (default: all)")
    parser.add_argument("--force", action="store_true", help="Force download even if files exist")
    parser.add_argument("--chromosomes", nargs="+", 
                     help="Specific chromosomes to download (format: 1 2 3 ... X Y MT). Omit for all chromosomes.")
    parser.add_argument("--threads", type=int, default=4,
                     help="Maximum number of parallel downloads (default: 4)")
    
    args = parser.parse_args()
    
    # Validate chromosomes input if provided
    valid_chroms = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
    if args.chromosomes:
        for chrom in args.chromosomes:
            if chrom not in valid_chroms:
                print(f"Error: Invalid chromosome '{chrom}'. Valid options are: {', '.join(valid_chroms)}")
                sys.exit(1)
    
    # Initialize downloader
    downloader = SmartDatabaseDownloader(
        force_download=args.force,
        chromosomes=args.chromosomes,
        max_threads=args.threads
    )
    
    # Run the download process
    downloader.run(args.datasets)

if __name__ == "__main__":
    main()