##############################################################
# SCRIPT 45: Predicting bAge - Bernabeu (2023)
##############################################################
# Copied from Bernabeu GitHub (cage_bage)
# Setup
rm(list=ls())
library(tidyverse)
library(survival)

##############################################################
# Load DNAm
load('../GOTO_Data/GOTO_methData-filtered.Rdata')
methData

# Betas
betas <- as.data.frame(assay(methData))
betas <- betas %>% rownames_to_column(var='cpg')

# Path to methylation table
write_csv(betas, file='../GOTO_Data/Clocks/Bern/methylationTable.csv')

methylationTable <- "../GOTO_Data/Clocks/Bern/methylationTable.csv"
methylationTable_format <- "csv"

##############################################################
# Phenotype table
# targets
targets <- as.data.frame(colData(methData))
targets <- targets %>% 
  select(sample_ID = Basename, age, sex) %>% 
  mutate(sex = ifelse(sex == 'female', 1, 0), tte = NA, death = 0)

write_csv(targets, file="../GOTO_Data/Clocks/Bern/phenotypeTable.csv")

phenotypeTable <- "../GOTO_Data/Clocks/Bern/phenotypeTable.csv"
phenotypeTable_format <- "csv"

##############################################################
# GrimAge results

grimageTable <- "../GOTO_Data/Clocks/grimAge/Grimage_output_GOTO.csv"
grimageTable_format <- "csv"

# Colnames
ageColname <- "age"
sexColname <- "sex"
tteColname <- "tte"
deathColname <- "death"

# Load in model coefficients
coefficients <- read.delim("bage_coefficients.tsv")

# Load in CpG coefficients
cpgs <- read.delim("cpg_episcore_weights.tsv")

# Load DNAm data
data <- read.csv(methylationTable, sep = ",", row.names = 1)

# Load pheno data
pheno <- read.csv(phenotypeTable, sep = ",", row.names = 1)

# Load grimage data
grim <- read.csv(grimageTable, sep = ",", row.names = 1)
row.names(grim) <- grim$SampleID

# Inspect
data[1:5, 1:5]
pheno[1:5, ]
grim[1:5, ]

colnames(data) <- gsub('X', '', colnames(data))

##############################################################
# QC and data prep for Episcore calculation
# Subset
pheno <- pheno[!is.na(pheno[tteColname]) & pheno[tteColname]>0,]
data <- data[,which(colnames(data) %in% rownames(pheno))]
grim <- grim[which(rownames(grim) %in% rownames(pheno)),]
data <- data[,rownames(pheno)]
grim <- grim[rownames(pheno),]

# Subset CpGs
coef <- data[intersect(rownames(data), cpgs$CpG_Site),]

# Scale
ids <- colnames(coef)
scaled <- apply(coef, 1, function(x) sd(x, na.rm = T)) 

coef <-  if(range(scaled)[1] == 1 & range(scaled)[2] == 1) { 
  coef
} else { 
  coef_scale <- apply(coef, 1, scale)
  coef_scale <- t(coef_scale)
  coef_scale <- as.data.frame(coef_scale)
  colnames(coef_scale) <- ids
  coef_scale
}

# Idnetify missing
coef <- if (nrow(coef)==length(unique(cpgs$CpG_Site))) { 
  message("All sites present")
  coef
} else { 
  missing_cpgs = cpgs[-which(cpgs$CpG_Site %in% rownames(coef)), 
                      c("CpG_Site", "Mean_Beta_Value")]
  message(paste(length(unique(missing_cpgs$CpG_Site)), 
                "unique sites are missing", sep = " "))
  mat = matrix(nrow=length(unique(missing_cpgs$CpG_Site)), 
               ncol = ncol(coef))
  row.names(mat) <- unique(missing_cpgs$CpG_Site)
  colnames(mat) <- colnames(coef) 
  mat[is.na(mat)] <- 1
  missing_cpgs1 <- if (length(which(duplicated(missing_cpgs$CpG_Site))) > 1) { 
    missing_cpgs[-which(duplicated(missing_cpgs$CpG_Site)),]
  } else {
    missing_cpgs
  }  
  ids = unique(row.names(mat))
  missing_cpgs1 = missing_cpgs1[match(ids,missing_cpgs1$CpG_Site),]
  mat = mat*missing_cpgs1$Mean_Beta_Value
  coef = rbind(coef,mat)
} 

# Impute mean values
na_to_mean <-function(methyl) {
  methyl[is.na(methyl)] <- mean(methyl, na.rm=T)
  return(methyl)
}
coef <- t(apply(coef,1,function(x) na_to_mean(x)))

##############################################################
# Calculate Episcores
loop <- unique(cpgs$Predictor)
out <- data.frame()
for(i in loop){ 
  tmp=coef[intersect(row.names(coef),cpgs[cpgs$Predictor %in% 
                                            i,"CpG_Site"]),]
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  if(nrow(tmp_coef) > 1) { 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_Site),]
    out[colnames(coef),i]=colSums(tmp_coef$Coefficient*tmp)
  } else {
    tmp2 = as.matrix(tmp)*tmp_coef$Coefficient 
    out[colnames(coef),i] = tmp2[,1]
  }
} 

# Save
write.table(out, 
            "../GOTO_Data/Clocks/Bern/episcore_projections.tsv", 
            sep = "\t", quote = FALSE)

# Scale GrimAge components
samples <- rownames(out)
grim_pred <- grim[samples, c("DNAmGrimAge"), drop = FALSE]
grim <- grim[samples, c("DNAmadm", "DNAmB2M", 
                        "DNAmCystatin_C", "DNAmGDF_15", 
                        "DNAmleptin", "DNAmPACKYRS", 
                        "DNAmpai_1", "DNAmTIMP_1")]

scaled_grim <- apply(grim, 2, function(x) sd(x, na.rm = T)) 
ids <- colnames(grim)

##############################################################
# bAge prediction
names(pheno)[names(pheno) == ageColname] <- "Age"
names(pheno)[names(pheno) == tteColname] <- "TTE"
names(pheno)[names(pheno) == deathColname] <- "Dead"
names(pheno)[names(pheno) == sexColname] <- "Sex"

# Merge
scores <- cbind(pheno[samples, c("TTE", "Dead", "Age", "Sex")], 
                grim, out)

# Subset
scores <- scores[, coefficients$Variable]

# Calculate
scores <- t(scores)
pred <- scores * coefficients[,"Coefficient"]
pred_pp <- colSums(pred)

# Scale 
scale_pred <- function(x, mean_pred, sd_pred, 
                       mean_test, sd_test) { 
  scaled <- mean_test + (x - mean_pred)*(sd_test/sd_pred)
  return(scaled)
}

scale_Z <- function(x, mean_pred, sd_pred) { 
  scaled <- (x - mean_pred)/sd_pred
  return(scaled)
}

mean_pred <- mean(pred_pp)
mean_test <- mean(pheno$Age) # Mean age in testing data
sd_pred <- sd(pred_pp)
sd_test <- sd(pheno$Age) # SD age in testing data
pred_pp_Z <- scale_Z(pred_pp, mean_pred, sd_pred)
pred_pp_scaled <- scale_pred(pred_pp, mean_pred, sd_pred, 
                             mean_test, sd_test)

# Make results df
pred_df <- data.frame(pred_pp_Z, pred_pp_scaled, 
                      grim_pred, pheno[samples, 
                                       c("Age", "Sex", "TTE", "Dead")])
names(pred_df) <- c("bAge", "bAge_Years", "GrimAge", 
                    "Age", "Sex", "TTE", "Dead")

# Save
write.table(data.frame("Sample" = rownames(pred_df), 
                       pred_df), 
            file = paste0("../GOTO_Data/Clocks/Bern/bage_predictions.tsv"), 
            quote = FALSE, sep = "\t", row.names = FALSE)

##############################################################