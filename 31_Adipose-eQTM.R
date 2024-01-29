##############################################################
# SCRIPT 31: Adipose - eQTM analysis
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
deg <- read_csv('../GOTO_Data/eQTM/Fat_DEG.csv')
gene_list <- unique(deg$ens)

# Inspect
length(gene_list)
head(gene_list)

##############################################################
# Load DNAm
load('../GOTO_Data/GOTO_results-top-fat.Rdata')
load('../GOTO_Data/GOTO_methData-filtered.Rdata')

# Subset
methData <- methData[ , methData$tissue == 'fat']

# Phenotype data
targets <- as.data.frame(colData(methData)) %>% 
  dplyr::select(IOP2_ID, 
                timepoint,
                age, sex, plate,
                array_row) 

# ID list
ID_list <- targets %>% mutate(
  ID = paste0(IOP2_ID, '_', ifelse(timepoint=="before", 0, 1))
) %>% dplyr::select(ID)

# DNAm values
betas <- as.data.frame(t(assay(methData)))

##############################################################
# RNAseq
load('../GOTO_Data/RNA/expData_adipose.RData')

goto_exp <- expData %>%  
  dplyr::select(IOP2_ID, timepoint, final_plate,
                starts_with('ENS')) %>% 
  mutate(
    ID = str_c(IOP2_ID, '_', timepoint)
  )

exp_df <- left_join(ID_list, goto_exp, by = 'ID')

ens2gene <- cinaR::grch38

# Analyse
for(i in deg$cpg){
  # Keep overlapping genes
  keep <- gene_list %in% colnames(exp_df)
  gene_list <- gene_list[keep]
  
  # Create df for analysis
  lm_df <- cbind(targets,
                 betas %>% dplyr::select(i),
                 exp_df)
  
  # Analyse
  models <- lapply(genes, function(x){
    lmer(substitute(i ~ j + age + sex +
                      plate + array_row + final_plate + 
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
write_csv(df_out, file = '../GOTO_Data/eQTM/geneList-fat.csv')

##############################################################
  