
## Short tutorial for BiMiCo package

# fastq files: 

rawfqs <- "/media/peter1/Data03/WORK_2020/BIOMAP/Demo_data_MAARS/"

# sample data: "/media/peter1/Data03/WORK_2020/BIOMAP/BiMiCo/demo_pheno.csv"

demo_pheno <- read.csv("/media/peter1/Data03/WORK_2020/BIOMAP/BiMiCo/demo_pheno.csv", row.names=1)

# taxonomy database: 

txset <- "/media/peter1/Data03/WORK_2020/BIOMAP/BiMiCo/gg_13_8_train_set_97.fa.gz"

#txset2 <- "/media/peter1/Data03/WORK_2020/BIOMAP/BiMiCo/rdp_species_assignment_14.fa.gz"

# INSTALLATION



library(BiMiCo)

##### A - SEQUENCE PREPROCESSING #####

library(ggplot2)

# 1. Filter 454 reads, write to subdir

std_454_preprocess(rawfqs,
                   "filt_fastqs")

# 2. inspect reads (e.g. check 3rd file in dir)

insp_reads(rawfqs,
           "/media/peter1/Data03/WORK_2020/BIOMAP/Demo_data_MAARS/filt_fastqs",
           3)


# 3. create ASV table

asvtab <- std_454_ASVtab("/media/peter1/Data03/WORK_2020/BIOMAP/Demo_data_MAARS/filt_fastqs")

# 4. Assign taxonomy

taxtab <- asgntax(asvtab, txset, revcomp = T)

#spectab <- assignSpecies(rownames(asvtab), txset2)

# 5. Store data in Phyloseq object

# import example phenodata

#exmpl_phenodata <- read.csv("D:/WORK_2020/BIOMAP/BiMiCo_01/BiMiCo/data/exmpl_phenodata.csv", row.names=1)

tmp1 <- data.frame(rownames(demo_pheno[12:15, ]), colnames(asvtab))

rownames(demo_pheno) <- colnames(asvtab)

myphy <- create_phylo(asvtab, taxtab, demo_pheno)

# 6. filter non-bacterial sequences

myphy <- rm_nonbac(myphy)

# n+1. decontamination (Optional if neg. controls are available)


##### B - STATISTICAL ANALYSES #####

# 1. Filter low abundance and prevalence

#myphy_filt 
  
prevalencedf = data.frame(Prevalence = (apply(X = otu_table(myphy),
                                             MARGIN = 1,
                                             FUN = function(x){sum(x > 0)})),
                          TotalAbundance = taxa_sums(myphy),
                          tax_table(myphy))

plot(prevalencedf$TotalAbundance, prevalencedf$Prevalence)


# 2. Create ordination plot

pcoa_ph(myphy, "Type") # color by patient gender

pcoa_ph(myphy, "PlateNum") # color by patient geographic location

# 3. Barplot

plot_bar(myphy, fill = "Genus") + theme(legend.position="bottom")

plot_bar(myphy, fill = "Order", x="Type", facet_grid = "PlateNum") + theme(legend.position="bottom")

plot_bar(myphy, fill = "Species") + theme(legend.position="bottom")

##### more phyloseq plots...




