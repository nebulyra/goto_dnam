##############################################################
# SCRIPT 5: DNAmArray - Functional normalization
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(DNAmArray)
library(minfi)

##############################################################
# Load data
load('../GOTO_Data/Processing/GOTO_RGset-unfiltered.Rdata')
load('../GOTO_Data/Processing/GOTO_targets-unfiltered.Rdata')

# Load annotation
anno <- read_tsv(
  "/exports/molepi/users/ljsinke/LLS/Shared_Data/Manifests/EPIC.hg19.manifest.tsv.gz")

##############################################################
# Functional normalization
# Screeplot
pc <- screeplot(RGset)

# Normalize
GRset <- preprocessFunnorm(
  RGset, 
  nPCs=4)

# Save
save(
  RGset, 
  file = "../GOTO_Data/Processing/GOTO_RGset-unfiltered.Rdata")

save(
  GRset,
  file="../GOTO_Data/Processing/GOTO_GRset-unfiltered.Rdata"
)

##############################################################
# Summarized Experiment
# Betas
betas <- getBeta(
  RGset, 
  type="Illumina")

# Subset anno
anno <- anno %>% 
  dplyr::select(
    cpg = probeID,
    chr = CpG_chrm,
    start = CpG_beg,
    end = CpG_end) %>% 
  dplyr::filter(
    cpg %in% rownames(betas))
dim(anno)

# Save unfiltered
methData_unfiltered <- SummarizedExperiment(
  assays=SimpleList(betas=assay(GRset)), 
  colData=targets,
  rowData=anno)

save(
  methData_unfiltered,
  file="../GOTO_Data/Processing/GOTO_methData-unfiltered.Rdata")

targets <- as.data.frame(colData(methData_unfiltered))

save(
  targets,
  file="../GOTO_Data/Processing/GOTO_targets-unfiltered.Rdata"
)

##############################################################