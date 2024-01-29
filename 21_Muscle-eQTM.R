##############################################################
# SCRIPT 21: Muscle - eQTM analysis
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
deg <- read_csv('../GOTO_Data/eQTM/Muscle_DEG.csv')
gene_list <- unique(deg$ens)

# Inspect
length(gene_list)
head(gene_list)

# RNAseq functions
source('../GOTO_Data/RNA/goto.rnaseq.functions.R')

# Paths
pathIN_dat <- "../GOTO_Data/RNA/merge.gene.counts_biopet_13052016.RData"
pathIN_cov <- "../GOTO_Data/RNA/muscle_QC_covariates_filesv2.csv"

# Settings
filt.samp <- "tissue_muscle|qc_sexswitch|qc_multdim2|qc_rep1|complete_pairs"

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
# Analyse
for(i in deg$cpg){
  # Keep overlapping genes
  keep <- gene_list %in% colnames(exp_df)
  gene_list <- gene_list[keep]
  
  # Create df for analysis
  lm_df <- cbind(targets,
                 betas %>% dplyr::select(i),
                 goto_exp)
  
  # Analyse
  models <- lapply(genes, function(x){
    lmer(substitute(i ~ j + age + sex +
                      plate + array_row + flowcell + 
                      (1|IOP2_ID),
                    list(i = as.name(i),
                         j = as.name(x))), data = lm_df)
  })
  
  if (length(models) > 0) {
    sym <- ens2gene$symbol[match(genes, ens2gene$ensgene)]
    
    p <- coef(summary(models[[1]]))[2,5]
    coef <- coef(summary(models[[1]]))[2,1]
    if (length(models) > 1) {
      for(k in 2:length(models)){
        p <- c(p, coef(summary(models[[k]]))[2,5])
        coef <- c(coef, coef(summary(models[[k]]))[2,1])
      }
      
    }
    
    df <- data.frame(
      cpg = i,
      gene = sym,
      ens = genes,
      coef = coef,
      p = p
    )
  }
  
  if (i == names(cpg_ranges)[1]){
    df_out <- df
  } else {
    df_out <- rbind(df_out, df)
  }
}
# Adjust p-values
df_out$padj <- p.adjust(df_out$p, method = 'fdr')

# Save
write_csv(df_out, file = '../GOTO_Data/eQTM/muscle_eQTM.csv')

##############################################################
  
  
  
  