##############################################################
# SCRIPT 35: Blood - Cell type changes
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

cells <- c("cc_blood_meas_eos", "cc_blood_meas_baso",    
           "cc_blood_meas_neut", "cc_blood_meas_lymph",
           "cc_blood_meas_mono", "cc_blood_music_CD4T", 
           "cc_blood_music_CD8T", "cc_blood_music_claM", 
           "cc_blood_music_cMOP", "cc_blood_music_CMP",
           "cc_blood_music_ery", "cc_blood_music_hMDP", 
           "cc_blood_music_immB", "cc_blood_music_interM",
           "cc_blood_music_matureN", "cc_blood_music_metaN",
           "cc_blood_music_myeN", "cc_blood_music_naiB",
           "cc_blood_music_nonM", "cc_blood_music_plasma", 
           "cc_blood_music_preM","cc_blood_music_regB",
           "cc_blood_music_toxiNK", "cc_blood_idol_CD8T",
           "cc_blood_idol_CD4T", "cc_blood_idol_NK",
           "cc_blood_idol_Bcell", "cc_blood_idol_Mono",
           "cc_blood_idol_Neu", "cc_blood_ext_Bas",
           "cc_blood_ext_Bmem", "cc_blood_ext_Bnv",
           "cc_blood_ext_CD4mem", "cc_blood_ext_CD4nv",
           "cc_blood_ext_CD8mem", "cc_blood_ext_CD8nv",
           "cc_blood_ext_Eos", "cc_blood_ext_Mono",
           "cc_blood_ext_Neu", "cc_blood_ext_NK", 
           "cc_blood_ext_Treg")

targets <- targets %>% filter(tissue == 'blood')

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
write_csv(res, file='ST15.csv')

##############################################################

