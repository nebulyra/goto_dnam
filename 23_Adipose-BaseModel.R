##############################################################
# SCRIPT 23: Adipose - Base model
##############################################################
# Setup
rm(list=ls())
library(DNAmArray)
library(SummarizedExperiment)
library(tidyverse)
library(limma)
library(sva)
library(methylGSA)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(minfi)
library(bacon)
library(lattice)
library(ggrepel)
library(irlba)
library(BiocParallel)
library(FDb.InfiniumMethylation.hg19)
library(ggfortify)
library(GenomicRanges)
library(ggrepel)
library(irlba)

##############################################################
# Load data
load("../GOTO_Data/GOTO_targets-filtered.Rdata")
load('../GOTO_Data/Processing/GOTO_methData-filtered.Rdata')

# Subset muscle
methData <- methData[ , methData$tissue == 'fat']
targets <- targets %>% filter(tissue == "fat")

##############################################################
# PC calculation
# Complete betas
complete_betas <- na.omit(assay(methData))
dim(complete_betas)

# Calculate PCs
pc <- prcomp_irlba(
  t(complete_betas), 
  n=5)

# Inspect
summary(pc)
dim(pc$x)
head(pc$x)

# Merge
targets <- as.data.frame(colData(methData))
targets <- cbind(
  targets, 
  pc$x)
colData(methData) <- DataFrame(targets)

##############################################################
# Linear mixed model
# Formula
formula <- ~timepoint + age + sex + smoking + plate + array_row + 
  PC1 + PC2 + PC3 + PC4 + PC5 

# Design matrix
design <- model.matrix(formula, 
                       data=colData(methData))

# Correlation for random effect
dupcor <- duplicateCorrelation(assay(methData),
                               design,
                               block = colData(methData)$IOP2_ID)

# Fit models
fit <- lmFit(assay(methData), design,
             block = colData(methData)$IOP2_ID,
             correlation = dupcor$consensus.correlation)

# Examine coefficients
str(fit$coefficients)

# Save
coef <- fit$coefficients[, 2]
se <- fit$stdev.unscaled[, 2] * fit$sigma
tstat <- coef / se
pval <- 2 * pt(-abs(tstat), fit$df.residual)
n <- ncol(design) + fit$df.residual

##############################################################
# Bacon
bc <- bacon::bacon(teststatistics = tstat,
                   effectsizes = coef,
                   standarderrors = se,
                   verbose = TRUE)
bc

# Estimate inflation and bias
bacon::estimates(bc)
bacon::inflation(bc)
bacon::bias(bc)

# Save bacon adjusted values
pval <- bacon::pval(bc)
tstat <- bacon::tstat(bc)

# Rerun bacon
bc <- bacon::bacon(teststatistics = tstat,
                   effectsizes = coef,
                   standarderrors = se,
                   verbose = TRUE)
bc

##############################################################
# Output
# Adjust p-values for multiple testing
padj_fdr <- p.adjust(pval, method="fdr")

# Save
limma_base <- data.frame(cpg = rownames(fit$coefficients), 
                         beta = coef, SE = se, 
                         p = pval, padj_fdr = padj_fdr, 
                         t = tstat, N = n)

# Inspect top hits
limma_base %>% arrange(padj_fdr, p) %>% 
  dplyr::select(cpg, beta, p, padj_bonf, padj_fdr, N) %>% head()

top_cpgs <- limma_base %>% filter(padj_fdr <= 0.05)

# Save
save(top_cpgs, file='../GOTO_Data/GOTO_results-top-fat.Rdata')
save(limma_base, file='../GOTO_Data/GOTO_results-full-fat.Rdata')

##############################################################
