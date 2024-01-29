##############################################################
# SCRIPT 10: Non-response analysis
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(lme4)
library(lmerTest)

##############################################################
# Load data
load("../GOTO_Data/GOTO_targets-filtered.Rdata")

# Load metabolic trait data 
pheno_data <- read_csv('../GOTO_Data/Pheno/GOTO_PhenoData.csv')

##############################################################
# Combine
# Pheno variables to numeric
pheno_data$age <- as.numeric(pheno_data$age)
pheno_data$sex <- ifelse(pheno_data$sex == 'male', 0, 1)
pheno_data$bmi <- as.numeric(pheno_data$bmi)

##############################################################
# Subset IDs
# Entire GOTO
before <- pheno_data %>% 
  filter(timepoint == 0) %>% 
  arrange(IOP2_ID)

after <- pheno_data %>% 
  filter(timepoint == 1, IOP2_ID %in% before$IOP2_ID) %>% 
  arrange(IOP2_ID)

all_ids <- unique(after$IOP2_ID)

# Muscle
targets_muscle <- targets %>% filter(tissue == 'muscle')
muscle_ids <- as.numeric(as.character(unique(targets_muscle$IOP2_ID)))

# Adipose
targets_fat <- targets %>% filter(tissue == 'fat')
fat_ids <- as.numeric(as.character(unique(targets_fat$IOP2_ID)))

# Blood
targets_blood <- targets %>% filter(tissue == 'fasted blood')
blood_ids <- as.numeric(as.character(unique(targets_blood$IOP2_ID)))

# Overlap
overlap_ids <- muscle_ids[muscle_ids %in% fat_ids]
overlap_ids <- overlap_ids[overlap_ids %in% blood_ids]

traits <- colnames(pheno_data)[9:18]
traits

##############################################################
# Non-response analysis
# All v. Overlap
for(i in traits){
  change_in <- pheno_data %>% 
    pivot_wider(id_cols='IOP2_ID',
              names_from='timepoint',
              values_from=i) %>% 
    mutate(delta = `1` - `0`)
  
  included <- (change_in %>% filter(IOP2_ID %in% overlap_ids))$delta
  excluded <- (change_in %>% filter(!IOP2_ID %in% overlap_ids))$delta
  
  fit <- t.test(included, excluded)

  out <- data.frame(trait = i,
                    delta = fit$estimate[1] - fit$estimate[2],
                    se = fit$stderr,
                    p = fit$p.value)
  
  if(i == traits[1]){
    res <- out
  } else {
    res <- rbind(res, out)
  }
  
  row.names(res) <- NULL
}

# Adjust p-values
res$padj <- p.adjust(res$p, method='fdr')
res_full <- res

##############################################################
# All v. Muscle
for(i in traits){
  change_in <- pheno_data %>% 
    pivot_wider(id_cols='IOP2_ID',
                names_from='timepoint',
                values_from=i) %>% 
    mutate(delta = `1` - `0`)
  
  included <- (change_in %>% filter(IOP2_ID %in% muscle_ids))$delta
  excluded <- (change_in %>% filter(!IOP2_ID %in% muscle_ids))$delta
  
  fit <- t.test(included, excluded)
  
  out <- data.frame(delta = fit$estimate[1] - fit$estimate[2],
                    se = fit$stderr,
                    p = fit$p.value)
  
  if(i == traits[1]){
    res <- out
  } else {
    res <- rbind(res, out)
  }
  
  row.names(res) <- NULL
}

# Adjust p-values
res$padj <- p.adjust(res$p, method='fdr')
res_full <- cbind(res_full, res)

##############################################################
# All v. Adipose
for(i in traits){
  change_in <- pheno_data %>% 
    pivot_wider(id_cols='IOP2_ID',
                names_from='timepoint',
                values_from=i) %>% 
    mutate(delta = `1` - `0`)
  
  included <- (change_in %>% filter(IOP2_ID %in% fat_ids))$delta
  excluded <- (change_in %>% filter(!IOP2_ID %in% fat_ids))$delta
  
  fit <- t.test(included, excluded)
  
  out <- data.frame(delta = fit$estimate[1] - fit$estimate[2],
                    se = fit$stderr,
                    p = fit$p.value)
  
  if(i == traits[1]){
    res <- out
  } else {
    res <- rbind(res, out)
  }
  
  row.names(res) <- NULL
}

# Adjust p-values
res$padj <- p.adjust(res$p, method='fdr')
res_full <- cbind(res_full, res)

##############################################################
# All v. Blood
for(i in traits){
  change_in <- pheno_data %>% 
    pivot_wider(id_cols='IOP2_ID',
                names_from='timepoint',
                values_from=i) %>% 
    mutate(delta = `1` - `0`)
  
  included <- (change_in %>% filter(IOP2_ID %in% blood_ids))$delta
  excluded <- (change_in %>% filter(!IOP2_ID %in% blood_ids))$delta
  
  fit <- t.test(included, excluded)
  
  out <- data.frame(delta = fit$estimate[1] - fit$estimate[2],
                    se = fit$stderr,
                    p = fit$p.value)
  
  if(i == traits[1]){
    res <- out
  } else {
    res <- rbind(res, out)
  }
  
  row.names(res) <- NULL
}

# Adjust p-values
res$padj <- p.adjust(res$p, method='fdr')
res_full <- cbind(res_full, res)

##############################################################
# Save
write_csv(res_full, file='ST2.csv')

##############################################################