##############################################################
# SCRIPT 3: DNAmArray - omicsPrint 
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(GenomicRanges)
library(SummarizedExperiment)
library(omicsPrint)
library(snpStats)

##############################################################
# CpGs with genotype information
# Annotation
snp_cpgs <- read_tsv(
  "../../LLS/Shared_Data/Manifests/EPIC.hg19.manifest.pop.tsv.gz")

# Format
snp_cpgs <- snp_cpgs %>% 
  dplyr::select(
    cpg = probeID,
    chr = CpG_chrm,
    start = CpG_beg,
    end = CpG_end,
    MASK_snp5_EUR) %>% 
  mutate(
    chr = substr(chr,4,5))

# Save SNPs
snp_cpgs <- (snp_cpgs %>% 
               dplyr::filter(MASK_snp5_EUR == TRUE))$cpg

print(paste0("There are ", length(snp_cpgs),
             " CpGs containing common European SNPs"))

##############################################################
# DNA methylation data
# Load
load("../GOTO_Data/Processing/GOTO_methData-unfiltered.Rdata")
methData_unfiltered

# Subset for SNPs
methData_snps <- methData_unfiltered[names(methData_unfiltered) %in% snp_cpgs,]

# Phenotype information
targets <- as.data.frame(colData(methData_snps))

# Check
table(rownames(targets) == colnames(methData_snps))

##############################################################
# omicsPrint - blood
# Subset
methData_blood <- methData_snps[ , 
                                 methData_snps$tissue == "fasted blood"]
methData_blood

# Phenotype information
targets_blood <- as.data.frame(colData(methData_blood))
dim(targets_blood)

# Beta values
betas_blood <- assay(methData_blood)
dim(betas_blood)

# Check
table(rownames(targets_blood) == colnames(betas_blood))

# Assumed relations
relations <- expand.grid(
  idx = colnames(methData_blood), 
  idy = colnames(methData_blood))

# Matching sample names
relations$sample_name_x <- targets_blood[relations$idx, "IOP2_ID"] 
relations$sample_name_y <- targets_blood[relations$idy, "IOP2_ID"]

# Relationship variable
relations$relation_type <- "unrelated"
relations$relation_type[relations$sample_name_x == relations$sample_name_y] <- "identical"
table(relations$relation_type)

# Genotype from omicsPrint
genotype <- beta2genotype(betas_blood, 
                          na.rm=TRUE,
                          minSep = 0.25,
                          minSize = 5,
                          centers = c(0.2, 0.5, 0.8),
                          assayName=NULL)

out <- alleleSharing(genotype, relations = relations)

# Identify mismatches
mismatches <- inferRelations(out, verbose=TRUE)

# Table of results
table(mismatches$relation, mismatches$predicted)

# Save objects
save(genotype, relations, out, mismatches, file = "../GOTO_Data/Processing/omicsPrint/GOTO_omicsPrint-blood.Rdata")

##############################################################
# omicsPrint - adipose
# Subset
methData_fat <- methData_snps[,methData_snps$tissue == "fat"]
methData_fat

# Phenotype information
targets_fat <- as.data.frame(colData(methData_fat))
dim(targets_fat)

# Beta values
betas_fat <- assay(methData_fat)
dim(betas_fat)

# Check
table(rownames(targets_fat) == colnames(betas_fat))

# Assumed relations
relations <- expand.grid(
  idx = colnames(methData_fat), 
  idy = colnames(methData_fat))

# Matching sample names
relations$sample_name_x <- targets_fat[relations$idx, "IOP2_ID"] 
relations$sample_name_y <- targets_fat[relations$idy, "IOP2_ID"]

# Relationship variable
relations$relation_type <- "unrelated"
relations$relation_type[relations$sample_name_x == relations$sample_name_y] <- "identical"
table(relations$relation_type)

# Genotype from omicsPrint
genotype <- beta2genotype(betas_fat, 
                          na.rm=TRUE,
                          minSep = 0.25,
                          minSize = 5,
                          centers = c(0.2, 0.5, 0.8),
                          assayName=NULL)

out <- alleleSharing(genotype, relations = relations)

# Identify mismatches
mismatches <- inferRelations(out, verbose=TRUE)

# Table of results
table(mismatches$relation, mismatches$predicted)

# Save objects
save(genotype, relations, out, mismatches, file = "../GOTO_Data/Processing/omicsPrint/GOTO_omicsPrint-fat.Rdata")

##############################################################
# omicsPrint - muscle
# Subset
methData_muscle <- methData_snps[,methData_snps$tissue == "muscle"]
methData_muscle

# Phenotype information
targets_muscle <- as.data.frame(colData(methData_muscle))
dim(targets_muscle)

# Beta values
betas_muscle <- assay(methData_muscle)
dim(betas_muscle)

# Check
table(rownames(targets_muscle) == colnames(betas_muscle))

# Assumed relations
relations <- expand.grid(
  idx = colnames(methData_muscle), 
  idy = colnames(methData_muscle))

# Matching sample names
relations$sample_name_x <- targets_muscle[relations$idx, "IOP2_ID"] 
relations$sample_name_y <- targets_muscle[relations$idy, "IOP2_ID"]

# Relationship variable
relations$relation_type <- "unrelated"
relations$relation_type[relations$sample_name_x == relations$sample_name_y] <- "identical"
table(relations$relation_type)

# Genotype from omicsPrint
genotype <- beta2genotype(betas_muscle, 
                          na.rm=TRUE,
                          minSep = 0.25,
                          minSize = 5,
                          centers = c(0.2, 0.5, 0.8),
                          assayName=NULL)

out <- alleleSharing(genotype, relations = relations)

# Identify mismatches
mismatches <- inferRelations(out, verbose=TRUE)

# Table of results
table(mismatches$relation, mismatches$predicted)

# Save objects
save(genotype, relations, out, mismatches, file = "../GOTO_Data/Processing/omicsPrint/GOTO_omicsPrint-muscle.Rdata")

##############################################################