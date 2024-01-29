##############################################################
# SCRIPT 13: Muscle - Estimate change in cell types
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(Biobase)
library(MuSiC)
library(SummarizedExperiment)
library(GEOquery)
library(lme4)
library(lmerTest)

##############################################################
# Load
load('../GOTO_Data/GOTO_targets-filtered.Rdata')

cells <- c("cc_musc_slowSM", "cc_musc_FAP", 
           "cc_musc_fastSM", "cc_musc_SAT", 
           "cc_musc_TNKmacro", "cc_musc_endo", 
           "cc_musc_smooth")

targets <- targets %>% filter(tissue == 'muscle')

for(i in cells){
  base_mean <- mean((targets %>% 
    filter(timepoint == 0) %>% select(i))[,1], na.rm)
  
  base_sd <- sd((targets %>% 
    filter(timepoint == 0) %>% select(i))[,1], na.rm)
  
  fit <- lmer(get(i) ~ timepoint + age + sex + (1|IOP2_ID),
              targets)
  
  out <- data.frame(
    cell = i,
    mean = base_mean,
    sd = base_sd,
    es = summary(fit)$coefficients[2,1],
    se = summary(fit)$coefficients[2,2],
    p = summary(fit)$coefficients[2,5]
  )
  
  if(i == cells[1]){
    res <- out
  } else {
    res <- rbind(res, out)
  }
}

# Adjust p-values
res$padj <- p.adjust(res$p, method='fdr')

# Save
write_csv(res, file='../GOTO_Data/Tables/ST3.csv')

##############################################################

