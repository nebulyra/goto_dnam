##############################################################
# SCRIPT 9: Descriptive tables
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(lme4)
library(lmerTest)

##############################################################
# Load data
load("../GOTO_Data/GOTO_targets-filtered.Rdata")
load('../GOTO_Data/Processing/GOTO_methData-filtered.Rdata')
load('../GOTO_Data/Processing/GOTO_RGset-filtered.Rdata')

# Load metabolic trait data 
pheno_data <- read_csv('../GOTO_Data/Pheno/GOTO_PhenoData.csv')

##############################################################
# Combine
# Pheno variables to numeric
pheno_data$age <- as.numeric(pheno_data$age)
pheno_data$sex <- ifelse(pheno_data$sex == 'male', 0, 1)
pheno_data$bmi <- as.numeric(pheno_data$bmi)

# Targets ID to double
targets$IOP2_ID <- as.numeric(as.character(targets$IOP2_ID))
targets$timepoint <- ifelse(targets$timepoint == 'after', 1, 0)
targets <- targets %>% select(-age, -sex, -bmi)

targets <- left_join(targets, 
                     pheno_data, by=c('IOP2_ID', 'timepoint'))

##############################################################
# Save
save(targets,
     file="../GOTO_Data/GOTO_targets-filtered.Rdata")

colData(methData) <- DataFrame(targets)
save(methData,
     file='../GOTO_Data/Processing/GOTO_methData-filtered.Rdata')

colData(RGset) <- DataFrame(targets)
save(RGset,
     file='../GOTO_Data/Processing/GOTO_RGset-filtered.Rdata')

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
# Entire GOTO population
# Table 1
for(i in traits){
  fit <- lmer(get(i) ~ timepoint + age + sex + (1|IOP2_ID),
              data = pheno_data)
  
  out <- data.frame(
    trait = i,
    delta = summary(fit)$coefficients[2,1],
    SE = summary(fit)$coefficients[2,2],
    p = summary(fit)$coefficients[2,5]
  )
  
  if(i == traits[1]){
    res <- out
  } else {
    res <- rbind(res, out)
  }
}

# Adjust p-values
res$padj <- p.adjust(res$p, method='fdr')
res_full <- res

##############################################################
# Overlap
# Table 1
for(i in traits){
  fit <- lmer(get(i) ~ timepoint + age + sex + (1|IOP2_ID),
              data = pheno_data %>% filter(IOP2_ID %in% overlap_ids))
  
  out <- data.frame(
    delta = summary(fit)$coefficients[2,1],
    SE = summary(fit)$coefficients[2,2],
    p = summary(fit)$coefficients[2,5]
  )
  
  if(i == traits[1]){
    res <- out
  } else {
    res <- rbind(res, out)
  }
}

# Adjust p-values
res$padj <- p.adjust(res$p, method='fdr')
res_full <- cbind(res_full, res)

##############################################################
# Muscle 
# Table 1
for(i in traits){
  fit <- lmer(get(i) ~ timepoint + age + sex + (1|IOP2_ID),
              data = pheno_data %>% filter(IOP2_ID %in% muscle_ids))
  
  out <- data.frame(
    delta = summary(fit)$coefficients[2,1],
    SE = summary(fit)$coefficients[2,2],
    p = summary(fit)$coefficients[2,5]
  )
  
  if(i == traits[1]){
    res <- out
  } else {
    res <- rbind(res, out)
  }
}

# Adjust p-values
res$padj <- p.adjust(res$p, method='fdr')
res_full <- cbind(res_full, res)

##############################################################
# Adipose
# Table 1
for(i in traits){
  fit <- lmer(get(i) ~ timepoint + age + sex + (1|IOP2_ID),
              data = pheno_data %>% filter(IOP2_ID %in% fat_ids))
  
  out <- data.frame(
    delta = summary(fit)$coefficients[2,1],
    SE = summary(fit)$coefficients[2,2],
    p = summary(fit)$coefficients[2,5]
  )
  
  if(i == traits[1]){
    res <- out
  } else {
    res <- rbind(res, out)
  }
}

# Adjust p-values
res$padj <- p.adjust(res$p, method='fdr')
res_full <- cbind(res_full, res)

##############################################################
# Blood
# Table 1
for(i in traits){
  fit <- lmer(get(i) ~ timepoint + age + sex + (1|IOP2_ID),
              data = pheno_data %>% filter(IOP2_ID %in% blood_ids))
  
  out <- data.frame(
    delta = summary(fit)$coefficients[2,1],
    SE = summary(fit)$coefficients[2,2],
    p = summary(fit)$coefficients[2,5]
  )
  
  if(i == traits[1]){
    res <- out
  } else {
    res <- rbind(res, out)
  }
}

# Adjust p-values
res$padj <- p.adjust(res$p, method='fdr')
res_full <- cbind(res_full, res)

##############################################################
# Save
write_csv(res_full, file='Table1.csv')

##############################################################
# Baseline characteristics
# Entire GOTO population
for(i in 9:18){
  before <- pheno_data %>% filter(timepoint==0)
  
  out <- data.frame(
    trait = colnames(pheno_data)[i],
    mean = mean(unlist(before[,i]), na.rm=T),
    sd = sd(unlist(before[,i]), na.rm=T)
  )
  
  if(i == 9){
    res <- out
  } else {
    res <- rbind(res, out)
  }
}

res_full <- res

# Overlap
for(i in 9:18){
  before <- pheno_data %>% 
    filter(timepoint==0, IOP2_ID %in% overlap_ids)
  
  out <- data.frame(
    mean = mean(unlist(before[,i]), na.rm=T),
    sd = sd(unlist(before[,i]), na.rm=T)
  )
  
  if(i == 9){
    res <- out
  } else {
    res <- rbind(res, out)
  }
}

res_full <- cbind(res_full, res)

# Muscle
for(i in 9:18){
  before <- pheno_data %>% 
    filter(timepoint==0, IOP2_ID %in% muscle_ids)
  
  out <- data.frame(
    mean = mean(unlist(before[,i]), na.rm=T),
    sd = sd(unlist(before[,i]), na.rm=T)
  )
  
  if(i == 9){
    res <- out
  } else {
    res <- rbind(res, out)
  }
}

res_full <- cbind(res_full, res)

# Adipose
for(i in 9:18){
  before <- pheno_data %>% 
    filter(timepoint==0, IOP2_ID %in% fat_ids)
  
  out <- data.frame(
    mean = mean(unlist(before[,i]), na.rm=T),
    sd = sd(unlist(before[,i]), na.rm=T)
  )
  
  if(i == 9){
    res <- out
  } else {
    res <- rbind(res, out)
  }
}

res_full <- cbind(res_full, res)

##############################################################
# Save
write_csv(res_full, file='ST1.csv')

##############################################################
