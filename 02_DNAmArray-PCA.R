##############################################################
# SCRIPT 2: DNAmArray - PCA
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
save("../GOTO_Data/Processing/GOTO_targets-unfiltered.Rdata")

# CpG annotation
anno <- read_tsv(
  "../../LLS/Shared_Data/Manifests/EPIC.hg19.manifest.tsv.gz")

##############################################################
# RGSet object
# BiocParallel
register(MulticoreParam(6))

# Read in IDATs
RGset <- read.metharray.exp(
  targets, 
  base = ".",
  extended=TRUE)
RGset

# Beta values
betas <- getBeta(
  RGset, 
  type="Illumina")
dim(betas)

# Density plot
densityPlot(
  RGset, 
  main="Beta density plot", 
  xlab="Beta values")

##############################################################
# PCA
# Principal Components (PCs)
pc <- prcomp_irlba(
  t(betas), 
  n=6)
summary(pc)

# Add to targets
targets <- cbind(
  targets, 
  pc$x)

# PCA plot
pca_plot <- targets %>% 
  ggplot(aes(
    x=PC1,
    y=PC2,
    label = ifelse(
      (PC1 < -30 & tissue != "fasted blood") |
        (PC1 > -60 & tissue == "fasted blood"), 
      as.character(IOP2_ID), NA)
  )) 

pca_plot +
  geom_point(
    aes(
      colour = tissue)) +
  geom_label_repel(
    aes(color = tissue),
    fill = "white"
  ) +
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("PCA plot coloured by tissue") +
  theme(
    axis.text = element_text(
      size=9, 
      color="#1B2021"),
    axis.title = element_text(
      size=12, 
      hjust=0.5, 
      color="#1B2021"),
    plot.title = element_text(
      size=16, 
      hjust=0.5,
      face="bold", 
      color = "#548687"),
    panel.background = element_rect(
      fill="white", color="#1B2021"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(
      fill="white"))

# Remove outliers
outliers <- targets %>% 
  dplyr::filter(
      (IOP2_ID == "61479" & tissue == "fasted blood") |
      (IOP2_ID == "60871" & tissue == "fasted blood") |
      (IOP2_ID == "61093" & tissue == "fat") |
      (IOP2_ID == "60822" & tissue == "muscle"))

targets <- targets[!(targets$Basename %in% outliers$Basename),]

##############################################################
# Save
# Basename
rownames(targets) <- targets$Basename

# Betas
betas <- betas[,colnames(betas) %in% rownames(targets)]
# Check
xtabs(~colnames(betas) == rownames(targets))

# Annotation
anno <- anno %>% filter(cpg %in% rownames(betas))
betas <- betas[rownames(betas) %in% anno$cpg,]

# SummarizedExperiment object
methData_unfiltered <- SummarizedExperiment(
  assays=SimpleList(betas=betas), 
  colData=targets,
  rowData=anno)
methData_unfiltered

# Subset
RGset <- RGset[,colnames(RGset) %in% rownames(targets)]
RGset
xtabs(~colnames(RGset) == rownames(targets))
RGset <- RGset[ , colnames(RGset) %in% colnames(betas)]

# Save
save(
  methData_unfiltered,
  file="../Processing/GOTO_methData-unfiltered.Rdata")

save(
  targets, 
  file="../Processing/GOTO_targets-unfiltered.Rdata")

save(
  RGset, 
  file = "../Processing/GOTO_RGset-unfiltered.Rdata")

##############################################################