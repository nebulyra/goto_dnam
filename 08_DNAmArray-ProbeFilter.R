##############################################################
# SCRIPT 8: DNAmArray - Probe filtering
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(DNAmArray)

# Load data
load("../GOTO_Data/Processing/GOTO_targets-unfiltered.Rdata")
load("../GOTO_Data/Processing/GOTO_RGset-unfiltered.Rdata")
load("../GOTO_Data/Processing/GOTO_GRset-unfiltered.Rdata")
load("../GOTO_Data/Processing/GOTO_methData-unfiltered.Rdata")

anno <- read_tsv(
  "../../LLS/Shared_Data/Manifests/EPIC.hg19.manifest.tsv.gz")
load("../../LLS/Shared_Data/ENCODE_Blacklist/ENCODEBlacklist_CpGomit-EPIC.RData")

##############################################################
# Filtering
# Low Intensity Probes
RGset <- probeFiltering(RGset)

# Detection P-value
betas <- reduce(GRset, RGset, what="beta")

# Filter anno
anno <- anno %>% 
  dplyr::select(
    cpg = probeID,
    chr = CpG_chrm,
    start = CpG_beg,
    end = CpG_end,
    everything()) %>% 
  dplyr::filter(
    cpg %in% rownames(betas))

# Sex Chromosomes
xy_cpgs <- (anno %>% 
              dplyr::filter(chr %in% c("X", "Y")))$cpg
betas <- betas[!rownames(betas) %in% xy_cpgs,]

# Encode Blacklist Regions
betas <- betas[!rownames(betas) %in% cpg_blacklist,]
dim(betas)

# Zhou probes
maskProbes <- anno[anno$MASK_general == TRUE,]$cpg
betas <- betas[!rownames(betas) %in% maskProbes,]

# 3IQR outlier Removal
iqr_dnam <- apply(betas, 1, function(x){
  iqr <- IQR(x, na.rm = TRUE)
  q1 <- quantile(x, na.rm=TRUE)[[2]]
  q3 <- quantile(x, na.rm=TRUE)[[4]]
  x <- ifelse((x <= q1 - (3*iqr) | x >= q3 + (3*iqr)), NA, x)
})
betas <- t(iqr_dnam)

# Missingness
perc_na <- rowSums(is.na(iqr_dnam))*100/ncol(iqr_dnam)
betas <- betas[,perc_na <= 95]

perc_na <- colSums(is.na(iqr_dnam))*100/nrow(iqr_dnam)
betas <- betas[perc_na <= 95,]

##############################################################
# Summarized Experiment
# Annotation
anno <- anno %>% 
  dplyr::filter(cpg %in% rownames(betas))

targets <- targets %>% 
  dplyr::filter(
    Basename %in% colnames(betas))

# Make object
methData <- SummarizedExperiment(
  assays=SimpleList(betas=betas), 
  rowData=anno, 
  colData=targets)

# Save filtered
RGset <- RGset[rownames(RGset) %in% rownames(methData),
               colnames(RGset) %in% colnames(methData)]
save(targets, file = "../GOTO_Data/GOTO_targets-filtered.Rdata")
save(methData, file="../GOTO_Data/GOTO_methData-filtered.Rdata")
save(RGset, file="../GOTO_Data/GOTO_RGset-filtered.Rdata")

##############################################################

