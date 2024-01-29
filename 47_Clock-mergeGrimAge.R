##############################################################
# SCRIPT 45: Predicting bAge - grimAge
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(foreach)

##############################################################
# Load betas
load('../Processing/GOTO_methData-unfiltered.Rdata')
betas <- assay(methData_unfiltered)

# Phenotype data
targets <- as.data.frame(colData(methData))

# Transform sex to Female
targets$Female <- ifelse(targets$sex == 'female', 1, 0)
xtabs(~targets$Female+targets$sex)

# Load needed CpGs
cpgs_in_GOTO <- as.vector(rownames(betas))
gold <- read.csv('../GOTO_Data/Clocks/grimAge/datMiniAnnotation3_Gold.csv')
row.names(gold)<-gold$CpG

# Determine missing
missing_cpgs <-setdiff(rownames(gold), cpgs_in_GOTO) 
missing_cpgs

# Making a df with the mean value of gold for number of 
# participants in our betas dataset
gold_missing <- gold[missing_cpgs,]
gold_missing <- as.numeric(gold_missing[,2])
gold_missing <- replicate(ncol(betas), gold_missing)
rownames(gold_missing) <- missing_cpgs

# Combining measured betas with gold dataset
betas_final <- rbind(betas, gold_missing)
betas_final2 <- betas_final[gold$CpG,]

# Save
write.csv(betas_final2, 
          file="../GOTO_Data/Clocks/grimAge/betas_imputed_with_gold_standard.csv")

##############################################################
### GrimAge
identical(colnames(betas_final2), rownames(targets))

grimage <- cbind(SampleID = colnames(betas_final2), 
                 Age = targets$age, 
                 Female = targets$Female, 
                 t(betas_final2))

# Initialize
r=19
n=dim(grimage)[1]

# Define groups
f <- rep(seq_len(ceiling(n/r)),
         each = r, length.out = n) 
ldf <- split(grimage, f = f)

# Save
for(i in unique(f)){
  tmp <- grimage[which(f %in% i),]
  write.csv(tmp, 
            paste0("../GOTO_Data/Clocks/grimAge/GrimAgeFile",
                   i, ".csv"))
}

##############################################################
# Merge
GrimAge<-foreach(i= 1:29, 
                 .combine="rbind") %do% {
  read.csv(paste0("GrimAge_output",i,".csv"),
           stringsAsFactors = F)
}

GrimAge<-GrimAge[,-1]
write.csv(GrimAge, "Grimage_output_GOTO.csv")
```

Calibrated GrimAge
```{r eval=FALSE}
GrimAge_scaled_to_LL<-foreach(i= 1:29, .combine="rbind") %do% {
  read.csv(paste0(main_dir,"/GrimAge_calculation_in_BIOS/cpgs_imputed_with_gold/GrimAge_output_scaled_to_LL/GrimAge_output",i,".csv"),
           stringsAsFactors = F)
}
GrimAge_scaled_to_LL<-GrimAge_scaled_to_LL[,-1]
write.csv(GrimAge_scaled_to_LL, paste0(main_dir,"/GrimAge_calculation_in_BIOS/cpgs_imputed_with_gold/Grimage_output_BIOS_scaled_to_LL.csv"))
```
