##############################################################
# SCRIPT 24: Adipose - Predict cell proportions with CIBERSORTx
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(Biobase)
library(MuSiC)

# Load data
load("../GOTO_Data/GOTO_targets-filtered.Rdata")
load("../GOTO_Data/GOTO_methData-filtered.Rdata")

##############################################################
# Load data from CIBERSORTx
load("../GOTO_Data/Cell_Counts/Fat/counts-targets.RData")

cells <- cells %>% 
  dplyr::select(IOP2_ID, timepoint, 
                cc_fat_adipo = adipocytes, 
                cc_fat_CD4T = cd4, 
                cc_fat_MVEC = MVEC, 
                cc_fat_m1 = m1_macro,
                cc_fat_m2 = m2_macro) %>% 
  mutate(IOP2_ID = as.factor(as.numeric(IOP2_ID)))

# Merge
targets <- left_join(targets, cells, by = c('IOP2_ID', 'timepoint'))
rownames(targets) <- targets$Basename

# Add to methData
check <- targets$Basename == colnames(methData)
xtabs(~check)

order <- colnames(methData)

targets <- targets[match(order, targets$Basename),]
rownames(targets) <- targets$Basename

colData(methData) <- DataFrame(targets)

# Save
save(targets, file="../GOTO_Data/GOTO_targets-filtered.Rdata")
save(methData, file="../GOTO_Data/GOTO_methData-filtered.Rdata")

##############################################################


