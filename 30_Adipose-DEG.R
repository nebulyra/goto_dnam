##############################################################
# SCRIPT 30: Adipose - DEGs near CpGs
##############################################################
# Setup
rm(list=ls())
library(tidyverse)

# Load data
load('../GOTO_Data/eQTM/fat_gene_cpg_pairs.Rdata')

# How many genes and CpGs?
length(unique(res$gene))
length(unique(res$cpg))

# Store genes
genes <- unique(res$gene)

# Read in DEG analysis result
deg <- read.table('DEGs_fat_postprandial_combined_new.txt')

# Filter for genes close to differentially methylated CpGs
deg <- deg %>% filter(GENE.SYMBOL %in% genes)

# Adjust for multiple testing
deg$padj <- p.adjust(deg$pvalue, method='fdr')
deg <- deg %>% filter(padj <= 0.05)

# Save FDR DEGs
fdr_deg <- unique(deg$GENE.SYMBOL)

# Table
deg <- deg %>% 
  select(logFC, gene=GENE.SYMBOL, pvalue, padj)
out <- left_join(res, deg, by='gene')
out <- out %>% 
  select(cpg, gene, ens, logFC, pvalue, padj)

# Save
write_csv(out, file='../GOTO_Data/eQTM/Fat_DEG.csv')






