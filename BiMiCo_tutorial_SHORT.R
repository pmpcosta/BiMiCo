
## Short tutorial for BiMiCo package basic functions

## The BiMiCo (Biomap Microbiome Core-tools) package offers single-step preprocessing of 16S microbial marker gene sequencing raw data files for 454 and Illumina platforms.
## The current tutorial and example data is aimed at the simple preprocessing of 454-generated fastq data. For more options and details, please see the "BiMiCo_tutorial_LONG" document.


# Install CRAN dependencies

install.packages("BiocManager")
install.packages("ggplot2")

# Install Bioconductor dependencies

BiocManager::install("Biobase")
BiocManager::install("Biostrings")
BiocManager::install("phyloseq")
BiocManager::install("dada2")

# Optional packages (required for "Long" tutorial/workflow)

install.packages("fastqcr")
install.packages("phangorn")
BiocManager::install("decontam")

# Load required packages
library(BiMiCo)
library(fastqcr)
library(ggplot2)


# Specify INPUTS:

# fastq files (directory containing raw fastq files intended for analysis exclusively, no other files): 
rawfqs <- "/media/peter1/Data03/WORK_2020/BIOMAP/Demo_data_MAARS/raw_fqs"

# sample phenodata (sample names in rows, sample variables in columns):
phedat <- read.csv("/media/peter1/Data03/WORK_2020/BIOMAP/Demo_data_MAARS/phenodata_maars_demo.csv", row.names=1)

# if not already contained in sample phenodata, please specify PRIMER TYPE (amplified region):
phedat$primer_type <- "V3-V4"

# taxonomy database to use with DADA2 (default is Greengenes v.13.5):
txset <- "/media/peter1/Data03/WORK_2020/BIOMAP/Demo_data_MAARS/gg_13_8_train_set_97.fa.gz"

# Specify OUTPUT folders:

# directory to write quality-filtered fastq files to:
filtered_fqs <- "/media/peter1/Data03/WORK_2020/BIOMAP/Demo_data_MAARS/filt_fqs"

# directory to write various outputs(ASV table, figures & graphs) to:
out_main <- "/media/peter1/Data03/WORK_2020/BIOMAP/Demo_data_MAARS/results/"


# Run single-command preprocessing for 454 reads. 
# Tip #1: trimming the length of 454 sequencing reads is highly recommended due to the drop-off in quality; please set trim_read_length accordingly, see "Long" tutorial for more details on QC/trimming
# Tip #2: Multithreading is enabled with the "mtthread=T" argument- not recommended when running from Rstudio.

toy_set <- bimico_454(rawfqs, txset, phedat, filtered_fqs, trim_read_length = 500, mtthread = F)

# 



