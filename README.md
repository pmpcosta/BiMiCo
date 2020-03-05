# BiMiCo

## Example for Microbiome Core methods

16S microbiome analysis toolbox

## Introduction

Analysis pipeline based around DADA2 and Phyloseq R packages for the processing of metagenome dataset  from raw read files to exploratory and comparative statistics. 

## Installation:


## Input files:

* directory of raw sequencing reads in fastq format
* sample/patient phenotype data (e.g. xls, csv table)
* taxonomic database for 16S rRNA-based taxonomic assignment, currently preferred: Greengenes (https://zenodo.org/record/158955/files/gg_13_8_train_set_97.fa.gz?download=1)

## Brief outline of the workflow:

1. Specify your path to fastq files to either
  + std_454_preprocess
  + std_illum_SE_preprocess
  + std_illum_PE_preprocess
    functions, according to sequencing platform (SE = single-end, PE = paired-end) *Currently extension must be "*.f*q"
    A subdirectory of input dir will be created, named using arg. 2 of "std_preprocess" functions, storing quality-filtered fastqs.
  
2. Inspect read quality with "insp_reads", takes as argument the path to raw fastq, path to filtered fastq, and index of file to inspect (e.g. "2" for 2nd fastq file). Generates two separate plots of before-after quality scores

3. Create ASV table from filtered fastq reads. Platform-dependent, so choose either
  + std_454_ASVtab
  + std_illum_SE_ASVtab
  + std_illum_PE_ASVtab
  function. Input is path to filtered fastqs, output chimera and duplicate-filtered ASV table with ASVs in ROWS / samples in COLS.
  
4. Identify taxonomy of ASVs in "asgn_tax" using the supplied dada2-compatible taxonomic database, e.g. RDP, SILVA or Greengenes (currently latter is preferred, since PicRust functional inference can be later used with Greengenes taxonomy). Input ASV table from step 3 and Greengenes database path (dl link in Input files above)

5. Store results in a Phyloseq object for further visualization and analysis, "create_phylo". Specify ASV table (step 3), taxonomy table (step 4) and sample phenotype data - object will be created without phenodata, but comparative statistics cannot be done. *Optionally add phylogenetic tree

6. Filter contamination using the "decontam" package (optional)

7. Various analyses can use Phyloseq functions directly or from Shiny-phyloseq, e.g. PCoA, alpha - beta diversity plots, species abundance barplots. In some cases ASV counts wil need to be transformed to relative abundance. Later even phyloseq-Deseq2 based differential relative abundance can be added  
