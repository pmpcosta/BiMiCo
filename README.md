---
# Short tutorial for BiMiCo package basic functions
---



The BiMiCo (Biomap Microbiome Core-tools) package offers single-step preprocessing of 16S microbial marker gene sequencing raw data files for 454 and Illumina platforms. \
The current tutorial and example data is aimed at the simple preprocessing of 454-generated fastq data. For more options and details, please see the "BiMiCo_tutorial_LONG" document. (pls contact me for appropriate demo dataset)

# Install dependencies

- Install CRAN dependencies

`install.packages("BiocManager")` \
`install.packages("ggplot2")`

- Install Bioconductor dependencies

`BiocManager::install("Biobase")` \
`BiocManager::install("Biostrings")` \
`BiocManager::install("phyloseq")` \
`BiocManager::install("dada2")`

# Load required packages

`library(BiMiCo)` \
`library(ggplot2)` \


# Specify INPUTS:

- fastq files (directory containing raw fastq files intended for analysis exclusively, no other files):

`rawfqs <- "./Demo_data/raw_fqs"`

- sample phenodata (sample names in rows, sample variables in columns):

`phedat <- read.csv("./Demo_data/phenodata_demo.csv", row.names=1)`

- if not already contained in sample phenodata, please specify **PRIMER TYPE** (amplified region):

`phedat$primer_type <- "V3-V4"`

- taxonomy database to use with DADA2 (default is Greengenes v.13.5):

`txset <- "./Demo_data/gg_13_8_train_set_97.fa.gz"`

# Specify OUTPUT folders:

- directory to write quality-filtered fastq files to:

`filtered_fqs <- "./Demo_data/filt_fqs"`

- directory to write various outputs(ASV table, figures & graphs) to:

`out_main <- "./Demo_data/results/"`


# Generate a summary report based on the desired variable column of phenodata (case/control, sample source, sequencing batch...)

- UNDER DEVELOPMENT: specify phyloseq object, variable of interest, output file basename- PDF of ordination and top 4 ASVs, extended in the future

`bimico_summary(toy_set, "clinical_group", "ex1")`




