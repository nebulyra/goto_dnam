##############################################################
# SCRIPT 1: DNAmArray - MethylAid
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(MethylAid)
library(BiocParallel)
library(DNAmArray)
library(minfi)
library(omicsPrint)
library(irlba)
library(FDb.InfiniumMethylation.hg19)
library(ggfortify)
library(GenomicRanges)
library(SummarizedExperiment)
library(ggrepel)

##############################################################
# Load data
# Phenotype data
load("../GOTO_Data/Sample_Sheets/GOTO_wave1-targets.Rdata")

# CpG annotation
anno <- read_tsv(
  "../../LLS/Shared_Data/Manifests/EPIC.hg19.manifest.tsv.gz")

##############################################################
# Split Data
# Create a Basename which points to IDAT files
targets <- targets %>% 
  unite(Basename, 
        c(SentrixBarcode_A, SentrixPosition_A), 
        sep="_", 
        remove=FALSE)

# Split by study
targets <- split(targets, targets$study)

# Get targets for GOTO samples
targets_goto <- targets[[1]]

# Format
targets_goto <- targets_goto %>% 
  mutate(
    DNA_labnr = factor(DNA_labnr),
    IOP2_ID = factor(as.numeric(old_ID)), 
    HMU_ID = factor(ID),
    timepoint = factor(timepoint, 
                       levels = c("before", "after")),
    tissue = factor(tissue,
                    levels = c("fasted blood", "fat", "muscle")),
    op_status = factor(group,
                       levels = c("partner", "offspring")),
    sex = factor(old_sex, 
                 levels = c("male", "female")),
    plate = factor(Sample_Plate),
    well = factor(Sample_Well),
    array_n = factor(as.character(SentrixBarcode_A)),
    array_row = as.numeric(substr(SentrixPosition_A,3,3)))

# Save
save(targets_goto,
     file="../GOTO_Data/Processing/GOTO_targets-unfiltered.Rdata")

##############################################################
# MethylAid
# BiocParallel
BPPARAM <- MulticoreParam(8)

# Summarize
sData_goto <- summarize(targets_goto, 
                        batchSize=50, 
                        BPPARAM=BPPARAM,
                        force=TRUE)

# Visualize QC plots - outliers were saved as a csv and uploaded below
visualize(sData_goto) 

# Save sData object
save(sData_goto, file="../Processing/MethylAid/GOTO_sData-wave1.Rdata")

##############################################################
# Remove outliers
# Format
anno <- anno %>% 
  dplyr::select(
    cpg = probeID,
    chr = CpG_chrm,
    start = CpG_beg,
    end = CpG_end,
    strand = probe_strand,
    gene_HGNC,
    MASK_general) %>% 
  mutate(
    chr = substr(chr,4,5))

# Load in outliers
outliers <- read_csv("../Processing/MethylAid/GOTO_MethylAid-Outliers.csv")

# Remove outliers
outlier_pair <- (targets_goto %>% 
                   filter(Basename %in% outliers$ID) %>% 
                   mutate(
                     pair = paste0(IOP2_ID, tissue)
                   ))$pair

outlier_pair

targets_wave1 <- targets_goto %>% 
  mutate(
    pair = paste0(IOP2_ID, tissue)
  ) %>% 
  filter(
    !pair %in% outlier_pair) %>% 
  dplyr::select(-pair)

print(paste0("Data on ", nrow(targets_wave1), " samples remains..."))

##############################################################
# Resent samples
# Load phenotype data of resent samples
load("../Sample_Sheets/GOTO_wave2-targets.Rdata")

# Format
targets_wave2 <- targets %>% 
  mutate(
    DNA_labnr = as.factor(DNA_labnr),
    IOP2_ID = as.factor(IOP2_ID),
    HMU_ID = as.factor(HMU_ID),
    tissue = as.factor(tissue),
    timepoint = factor(timepoint, 
                       levels = c("before", "after")),
    sex = factor(sex, 
                 levels=c("male", "female")),
    op_status = factor(op_status, 
                       levels = c("partner", "offspring")),
    plate = factor(plate),
    well = as.factor(paste0(rep(LETTERS[1:8], 3), rep(c("09", "10", "11"), each=8))),
    Basename = paste0(chip_n, "_", chip_position),
    array_n = factor(chip_n),
    array_row = as.numeric(substr(chip_position,3,3)))

# Combine waves
targets <- rbind(targets_wave1, targets_wave2)

# Make pair ID
targets <- targets %>% 
  mutate(
    pair = paste0(IOP2_ID, tissue)
  ) 

# Remove non-paired samples
unpaired <- targets %>% 
  group_by(pair) %>% 
  dplyr::summarize(n=n()) %>% 
  ungroup()

unpaired <- (unpaired %>% 
               filter(n<2))$pair

unpaired

targets <- targets %>% filter(
  (!pair %in% unpaired)) %>% 
  dplyr::select(-pair)

print(paste0("Data on ", nrow(targets), " samples remains..."))

# Save
save(targets,
     file="../GOTO_Data/Processing/GOTO_targets-unfiltered.Rdata")

##############################################################
# MethylAid on full dataset
# Biocparallel
BPPARAM <- MulticoreParam(8)

# Summarize
sData_goto <- MethylAid::summarize(
  targets, 
  batchSize=50, 
  BPPARAM = BPPARAM,
  force=TRUE)

# Visualize
visualize(sData_goto)

# Save
save(sData_goto, 
     file="../Processing/MethylAid/GOTO_sData-all.Rdata")

##############################################################