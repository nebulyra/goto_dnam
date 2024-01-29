##############################################################
# SCRIPT 29: Adipose - Nearby genes
##############################################################
# Setup
rm(list=ls())
library(edgeR)
library(cinaR)
library(limma)
library(tidyverse)
library(GenomicFeatures)
library(AnnotationHub)
library(SummarizedExperiment)

# Load data
load('../GOTO_Data/GOTO_results-top-fat.Rdata')

##############################################################
# Data on gene locations
# Load
load('../Shared_Data/geneRanges-75.Rdata')

##############################################################
# Data on CpG locations
# Make GRanges
cpg_ranges <- GRanges(seqnames=sub('chr', '', 
                                   top_cpgs$cpg_chr_hg19), 
                      IRanges(start=top_cpgs$cpg_start_hg19,
                              width=2,
                              names=top_cpgs$cpg))
cpg_ranges

##############################################################
# Genes within 100kb
# Initialize
distance.matrix <- matrix(NaN, 
                          nrow = length(gene_range), 
                          ncol = length(cpg_ranges), 
                          dimnames = list(names(gene_range), 
                                          names(cpg_ranges)))

# Store distances
for(i in 1:nrow(top_cpgs)){
  cpg <- cpg_ranges[i]
  calc.distance <- distance(cpg, gene_range, ignore.strand=T)
  distance.matrix[,i] <- calc.distance
}
dim(distance.matrix)

# Process
distance.matrix <- as.data.frame(distance.matrix)
colnames(distance.matrix) <- names(cpg_ranges)
distance.matrix$gene <- names(gene_range)

# List of genes
for(i in names(cpg_ranges)){
  df <- (distance.matrix %>% dplyr::select(gene, all_of(i)))
  df <- df[,2]<=100000
  out <- data.frame(
    gene=df$gene,
    cpg=i
  )
  
  if(i==names(cpg_ranges)[1]){
    res <- out
  } else {
    res <- rbind(res, out)
  }
}

# Save list
save(res, file='../GOTO_Data/eQTM/fat_gene_cpg_pairs.Rdata')

##############################################################
