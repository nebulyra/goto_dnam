##############################################################
# SCRIPT 41: Blood - eQTM analysis
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
library(lme4)
library(lmerTest)

##############################################################
# Gene expression data
# Read in DEGs close to CpGs
deg <- read_csv('../GOTO_Data/eQTM/Blood_DEG.csv')
gene_list <- unique(deg$ens)

# Inspect
length(gene_list)
head(gene_list)

# RNAseq functions
source('../GOTO_Data/RNA/goto.rnaseq.functions.R')

# Paths
pathIN_dat <- "../GOTO_Data/RNA/merge.gene.counts_biopet_13052016.RData"
pathIN_cov <- "../GOTO_Data/RNA/datasheet_RNAseq_blood_V2.csv"

# Settings
filt.samp <- "tissue_blood|qc_sexswitch|qc_multdim2|qc_rep1|complete_pairs"

# Read in RNAseq
goto_exp <- read.gotornaseq(pathIN_dat = pathIN_dat, 
                        pathIN_cov, 
                        filt.samp = filt.samp, 
                        type = 'voom-export', 
                        quiet = FALSE)
goto_exp <- goto_exp[["dat"]]

# Select variables & DEGs
goto_exp <- goto_exp %>% 
  dplyr::select(IOP2_ID = sampID2, intervention, age, sex, 
                flowcell, all_of(gene_list))


ens2gene <- cinaR::grch38

##############################################################
# Analysis
# Run eQTM models
for(i in gene_list){
  models <- lmer(get(i) ~ intervention + age + sex + 
                   flowcell + (1|IOP2_ID), data = goto_exp)
  
  sym <- ens2gene$symbol[match(i, ens2gene$ensgene)]
  
  p <- coef(summary(models))[2,5]
  coef <- coef(summary(models))[2,1]
  
  df <- data.frame(
    gene = sym,
    ens = i,
    coef = coef,
    p = p
  )
  
  if (i == gene_list[1]){
    df_out <- df
  } else {
    df_out <- rbind(df_out, df)
    print(df_out)
  }
}

# Adjust p-values
df_out$padj <- p.adjust(df_out$p, method = 'fdr')

# Save
write_csv(df_out, file = '../GOTO_Data/eQTM/blood_eQTM.csv')

##############################################################
  
  
  
  