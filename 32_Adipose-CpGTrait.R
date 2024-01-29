##############################################################
# SCRIPT 32: Adipose - Cpg trait associations
##############################################################
# Setup
rm(list=ls())
library(SummarizedExperiment)
library(tidyverse)
library(sva)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(minfi)
library(lattice)
library(ggrepel)
library(lme4)
library(lmerTest)

##############################################################
# DNA methylation
# Load DNAm data
load('../GOTO_Data/GOTO_methData-filtered.Rdata')

# Save for fat
methData <- methData[ , methData$tissue == 'fat']

# Load top CpGs
load('../GOTO_Data/GOTO_results-top-fat.Rdata')

# Phenotype data
targets <- as.data.frame(colData(methData))

# Betas from target CpGs
betas <- assay(methData)
betas <- betas[rownames(betas) %in% top_cpgs, ]
colnames(betas) <- targets$Basename

# Add to phenotype data
betas <- as.data.frame(t(betas))
rownames(betas) <- colnames(betas)
colnames(betas) <- rownames(betas)

betas$Basename <- rownames(betas)
targets <- left_join(targets, betas, by = 'Basename')
targets <- targets %>% 
  mutate(join = str_c(IOP2_ID, '_', 
                      ifelse(timepoint == "before", 0, 1))) %>% 
  dplyr::select(join, IOP2_ID, timepoint, sex, age, 
                plate, array_row, starts_with('cg'))

##############################################################
# Traits
traits <- c('bmi', 'wc', 'perc_body_fat', 'fasting_insulin',
            'leptin', 'adiponectin', 'interleukin_6', 'hdl_size',
            'hdl_c', 'systolic_bp')

# Analysis
for (i in traits){
  models <- lapply(
    top_cpgs, 
    function(x){
      lmer(substitute(j ~ get(i) + age + sex +
                        plate + array_row + (1|IOP2_ID),
                      list(j = as.name(x))),
           data=targets_blood)
    })
  
  out <- data.frame(
    trait = i, 
    cpg = top_cpgs[1],
    beta = summary(models[[1]])$coefficients[2,1],
    se = summary(models[[1]])$coefficients[2,2],
    pval = summary(models[[1]])$coefficients[2,5])
  
  for(k in 2:length(top_cpgs)){
    out <- rbind(
      out,
      data.frame(
        trait = trait,
        cpg = top_cpgs[k],
        beta = summary(models[[k]])$coefficients[2,1],
        se = summary(models[[k]])$coefficients[2,2],
        pval = summary(models[[k]])$coefficients[2,5]))
  }
  
  out$padj <- p.adjust(out$pval, method='fdr')
  
  if(i == files[1]){
    res <- out
  } else {
    res <- rbind(res, out)
  }
  
}

# Save
save(res, file='../GOTO_Data/Tables/ST14.Rdata')

##############################################################

