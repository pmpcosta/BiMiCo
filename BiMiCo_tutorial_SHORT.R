
## Short tutorial for BiMiCo package basic functions

## The BiMiCo (Biomap Microbiome Core-tools) package offers single-step preprocessing of 16S microbial marker gene sequencing raw data files for Illumina platforms.
## The current tutorial and example data is aimed at the simple preprocessing of fastq data.


# Install CRAN dependencies

install.packages("BiocManager")
install.packages("ggplot2")

# Install Bioconductor dependencies

BiocManager::install("Biobase")
BiocManager::install("Biostrings")
BiocManager::install("phyloseq")
BiocManager::install("dada2")

# Load required packages
library(BiMiCo)
library(ggplot2)


# 

