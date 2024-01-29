##############################################################
# SCRIPT 6: DNAmArray - Remove mismatches
##############################################################
# Setup
rm(list=ls())
library(SummarizedExperiment)
library(tidyverse)

##############################################################
# Muscle Sample Mixup
# Load data
load('../GOTO_Data/Processing/GOTO_targets-unfiltered.Rdata')
load('../GOTO_Data/Processing/GOTO_methData-unfiltered.Rdata')
load('../GOTO_Data/Processing/GOTO_RGset-unfiltered.Rdata')

# We want to swap `203548970088_R08C01` and `203548980011_R03C01`
# Inspect
targets[c('203548970088_R08C01','203548980011_R03C01'),]

# Save variables
person_60366 <- targets['203548970088_R08C01', 
                        c('IOP2_ID', 'sex', 'age', 'bmi', 'op_status')]
person_60365 <- targets['203548980011_R03C01', 
                        c('IOP2_ID', 'sex', 'age', 'bmi', 'op_status')]

# Swap
targets['203548970088_R08C01', 
        c('IOP2_ID', 'sex', 'age', 'bmi', 'op_status')] <- person_60365
targets['203548980011_R03C01', 
        c('IOP2_ID', 'sex', 'age', 'bmi', 'op_status')] <- person_60366

# Check
targets[c('203548970088_R08C01','203548980011_R03C01'),]

# Add to methData
colData(methData_unfiltered) <- DataFrame(targets)

##############################################################
# Fat Sample Removal
# Remove sample and pair
methData_unfiltered <- methData_unfiltered[,!colnames(methData_unfiltered) %in% c('203548970042_R06C01', '203548970042_R05C01')]
RGset <- RGset[,!colnames(RGset) %in% c('203548970042_R06C01', '203548970042_R05C01')]
targets <- targets %>% filter(!Basename %in% c('203548970042_R06C01', '203548970042_R05C01'))

##############################################################
# Remove non-compliers
# Save IDs
non_comply <- c("60227", "60220", "61792",
                "61311", "61283", "61284",
                "60474", "61763", "62238")
non_comply_basenames <- (targets %>% 
                           dplyr::filter(IOP2_ID %in% non_comply))$Basename

# Filter
targets <- targets %>% dplyr::filter(!IOP2_ID %in% non_comply)
methData_unfiltered <- methData_unfiltered[,!colnames(methData_unfiltered) %in% non_comply_basenames]
RGset <- RGset[,!colnames(RGset) %in% non_comply_basenames]

# Save
save(targets, file='../GOTO_Data/Processing/GOTO_targets-unfiltered.Rdata')
save(methData_unfiltered, file='../GOTO_Data/Processing/GOTO_methData-unfiltered.Rdata')
save(RGset, file='../GOTO_Data/Processing/GOTO_RGset-unfiltered.Rdata')

##############################################################