
## Short tutorial for BiMiCo package

# INSTALLATION

##### later

## BLA

# Load required packages
library(BiMiCo)
library(fastqcr)
library(ggplot2)
library(decontam)


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

# directory to write quality control reports by FastQC, pre- and post filtering to:
qcdir_pre <- "/media/peter1/Data03/WORK_2020/BIOMAP/Demo_data_MAARS/qc_reports_pre"
qcdir_post <- "/media/peter1/Data03/WORK_2020/BIOMAP/Demo_data_MAARS/qc_reports_post"

# directory to write various outputs(ASV table, figures & graphs) to:
out_main <- "/media/peter1/Data03/WORK_2020/BIOMAP/Demo_data_MAARS/results/"


# INPUT REQ: PRIMER TYPE; Batches?
# Outs A: filtered fastqs; qc reports; norm. matrix; phylos. object?
# Outs B: ordination, boxplots, diff. abundance etc.



##### PART A - SEQUENCE PREPROCESSING #####

# 1. Filter 454 reads, write to subdir

# 1.a inspect reads with fastqcr pre-filtering (individual fastq reports)
fastqc(fq.dir = rawfqs,
       qc.dir = qcdir_pre,
       threads = 1)
# (summary table)
fqc_pre <- qc_aggregate(qcdir_pre)

# 1.b Preprocessing reads with convenience functions (either 454 or Illumina single/paired-end fastqs)

# Excessive duplication and GC-content bias are expected in 16S sequencing data.
# As read quality deteriorates in the second half of 454 reads, we also apply trimming by length to avoid spurious annotations
prep_454(rawfqs,
         filtered_fqs, 
         mtthread = F,
         trim5end=10,
         trim_read_length = 500)

# 1.c inspect quality post-filtering 
# Note: MultiQC command-line tool is ideal for graphically summarizing a folder of FastQC reports
fastqc(fq.dir = filtered_fqs,
       qc.dir = qcdir_post,
       threads = 1)
# (summary table)
fqc_post <- qc_aggregate(qcdir_post)

# 3. create ASV table
asvtab <- asvtab_454(filtered_fqs, mtthread=F)

# clean ASV table from too short sequences
asvtab <- asvtab[ nchar(rownames(asvtab))>=50, ]

# 4. Assign taxonomy based on the dataset specified in 'txset'. Default is Greengenes 13.5 database
taxtab <- asgntax(asvtab, txset, revcomp = T, mtthread = F)

# 5. Save preprocessed data
# match phenodata rownames
rownames(phedat) <- gsub(".fastq", "", rownames(phedat))
colnames(asvtab) <- gsub("_filt", "", colnames(asvtab))
asvtab <- asvtab[ , order(colnames(asvtab))]
phedat <- phedat[ order(rownames(phedat)), ]
summary(rownames(phedat)==colnames(asvtab))

# Save all relevant data files
# Store data in Phyloseq object
phs <- create_phylo(asvtab, taxtab, phedat)

# save(phs, file="./phs.R")

##### PART B - STATISTICAL ANALYSES #####


# 2. Filter low abundance and prevalence

prevalencedf = data.frame(Prevalence = (apply(X = otu_table(phs),
                                             MARGIN = 1,
                                             FUN = function(x){sum(x > 0)})),
                          TotalAbundance = taxa_sums(phs),
                          tax_table(phs))

plot(prevalencedf$TotalAbundance, prevalencedf$Prevalence)


# 2. Create ordination plot

pcoa_ph(phs, "clinical_group") + geom_text_repel(aes(label=rownames(phs@sam_data)), size=5, alpha=0.5, hjust=2 ) # color by patient gender

pcoa_ph(phs, "Institution") # color by patient geographic location

# 3. Barplot

plot_bar(myphy, fill = "Genus") + theme(legend.position="bottom")

plot_bar(myphy, fill = "Order", x="Type", facet_grid = "Source") + theme(legend.position="bottom")

plot_bar(myphy, fill = "Species") + theme(legend.position="bottom")

# 4. boxplots

## On Relative abundance?? boxes should be comparable

# single ASV

selbact = 1

ggplot(as.data.frame(myphy@otu_table@.Data[ selbact, ]),
       aes(x=as.matrix(myphy@sam_data[ , "Type"]),
           y=myphy@otu_table@.Data[ selbact, ]))+ 
  geom_boxplot() +
  geom_jitter(shape=16, position = position_jitter(0.2)) +
  ggtitle(paste(myphy@tax_table@.Data[ selbact, 1 ],
                 myphy@tax_table@.Data[ selbact, 5 ],
                 myphy@tax_table@.Data[ selbact, 6 ],
                 myphy@tax_table@.Data[ selbact, 7 ],
                sep = "__"))







##### more phyloseq plots...




