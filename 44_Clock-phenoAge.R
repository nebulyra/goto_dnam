##############################################################
# SCRIPT 44: Predicting bAge - phenoAge
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(dnaMethyAge)
library(SummarizedExperiment)

##############################################################
# Load data
load("../GOTO_Data/GOTO_targets-filtered.Rdata")
load("../GOTO_Data/Processing/GOTO_methData-unfiltered.Rdata")
load("../GOTO_Data/Processing/GOTO_RGset-unfiltered.Rdata")

# Subset
methData <- methData_unfiltered[ , methData_unfiltered$tissue == 'fasted blood']
targets <- targets %>% filter(tissue == "fasted blood")
betas <- assay(methData)

# cAge predictors
bage <- c('LevineM2018')

# Save phenotype data
info <- targets %>% 
  select(Sample=Basename, Age=age, Sex=sex) %>% 
  mutate(Sex = ifelse(Sex == 'male', 'Male', 'Female'))

# Calculate
for(i in bage){
  age <- methyAge(betas, 
                  clock=i, 
                  age_info=info, 
                  fit_method='Linear', 
                  do_plot=TRUE)
  age_blood$Basename <- age_blood$Sample 
  age_blood <- full_join(age, targets, by='Basename')
  
  save(age, file=paste0('../GOTO_Data/Clocks/',i,'_blood.Rdata'))
}

##############################################################