


## Short tutorial for BiMiCo package

# fastq files: "D:/WORK_2020/BIOMAP/BiMiCo_01/BiMiCo/data/fastqs"

# sample data: "D:/WORK_2020/BIOMAP/BiMiCo_01/BiMiCo/data/exmpl_phenodata.csv"

# taxonomy database: "D:/WORK_2020/BIOMAP/WP2_core_methods/exampl_data_HMP/gg_13_8_train_set_97.fa.gz"


# INSTALLATION


library(BiMiCo)

##### A - SEQUENCE PREPROCESSING #####


# 1. Filter 454 reads, write subdir of fastq path

std_454_preprocess("D:/WORK_2020/BIOMAP/BiMiCo_01/BiMiCo/data/fastqs",
                   "filt_fastqs",
                   400)

# 2. inspect reads (e.g. check 3rd file in dir)

insp_reads("D:/WORK_2020/BIOMAP/BiMiCo_01/BiMiCo/data/fastqs",
           "D:/WORK_2020/BIOMAP/BiMiCo_01/BiMiCo/data/fastqs/filt_fastqs",
           3)


# 3. create ASV table

asvtab <- std_454_ASVtab("D:/WORK_2020/BIOMAP/BiMiCo_01/BiMiCo/data/fastqs/filt_fastqs")

# 4. Assign taxonomy

taxtab <- asgntax(asvtab, "D:/WORK_2020/BIOMAP/WP2_core_methods/exampl_data_HMP/gg_13_8_train_set_97.fa.gz")

# 5. Store data in Phyloseq object

# import example phenodata

exmpl_phenodata <- read.csv("D:/WORK_2020/BIOMAP/BiMiCo_01/BiMiCo/data/exmpl_phenodata.csv", row.names=1)

myphy <- create_phylo(asvtab, taxtab, exmpl_phenodata)

# 6. filter non-bacterial sequences

myphy <- rm_nonbac(myphy)





# n+1. decontamination (Optional if neg. controls are available)


##### B - STATISTICAL ANALYSES #####

# 1. Create ordination plot

pcoa_ph(myphy, "Sex") # color by patient gender

pcoa_ph(myphy, "Geo") # color by patient geographic location



##### more phyloseq plots...




