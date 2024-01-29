##############################################################
# SCRIPT 4: DNAmArray - Genotype check (LLS)
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(GenomicRanges)
library(SummarizedExperiment)
library(omicsPrint)
library(snpStats)

##############################################################
# Load data for genotype check
# Load genotype data
lls_plink <- read.plink(bed = '../GOTO_Data/Processing/Genotypes/LLS_Offspring_Partners_Final_36_130_Overlap.bed', 
                        bim = '../GOTO_Data/Processing/Genotypes/LLS_Offspring_Partners_Final_36_130_Overlap.bim', 
                        fam = '../GOTO_Data/Processing/Genotypes/LLS_Offspring_Partners_Final_36_130_Overlap.fam')

# Load omicsPrint output - blood
load('../GOTO_Data/Processing/omicsPrint/GOTO_omicsPrint-blood.Rdata')
dnam_blood <- as.data.frame(genotype)

print(paste0("We have data on ",
             nrow(dnam_blood), " SNPs from ",
             ncol(dnam_blood), " fasted blood samples"))

# Load omicsPrint output - muscle
load('../GOTO_Data/Processing/omicsPrint/GOTO_omicsPrint-muscle.Rdata')
dnam_muscle <- as.data.frame(genotype)

print(paste0("We have data on ",
             nrow(dnam_muscle), " SNPs from ",
             ncol(dnam_muscle), " skeletal muscle samples"))

# Load omicsPrint output - fat
load('../GOTO_Data/Processing/omicsPrint/GOTO_omicsPrint-fat.Rdata')
dnam_fat <- as.data.frame(genotype)

print(paste0("We have data on ",
             nrow(dnam_fat), " SNPs from ",
             ncol(dnam_fat), " adipose tissue samples"))

# Load key
key_id_gwas <- read_csv('../GOTO_Data/Processing/Genotypes/GOTO_GWASid-IOP2id_20210211.csv')

# Format
key_id_gwas$GWASnr <- sub('LLS', '', sub('__', '_', key_id_gwas$GWASnr))

##############################################################
# Bind to phenotype data
# Load phenotype data
load('../GOTO_Data/Processing/GOTO_targets-unfiltered.Rdata')

# Blood
key_blood <- targets %>% 
  dplyr::select(sample_ID = Basename,
                IOP2_ID,
                HMU_ID, 
                DNA_labnr) %>% 
  filter(sample_ID %in% colnames(dnam_blood)) %>% 
  mutate(IOP2_ID = as.numeric(as.character(IOP2_ID)))
key_blood <- left_join(key_blood, key_id_gwas, by = 'IOP2_ID')

# Muscle
key_muscle <- targets %>% dplyr::select(sample_ID = Basename,
                                        IOP2_ID,
                                        HMU_ID, 
                                        DNA_labnr) %>% 
  filter(sample_ID %in% colnames(dnam_muscle)) %>% 
  mutate(IOP2_ID = as.numeric(as.character(IOP2_ID)))
key_muscle <- left_join(key_muscle, key_id_gwas, by = 'IOP2_ID')

# Fat
key_fat <- targets %>% dplyr::select(sample_ID = Basename,
                                     IOP2_ID,
                                     HMU_ID, 
                                     DNA_labnr) %>% 
  filter(sample_ID %in% colnames(dnam_fat)) %>% 
  mutate(IOP2_ID = as.numeric(as.character(IOP2_ID)))
key_fat <- left_join(key_fat, key_id_gwas, by = 'IOP2_ID')

##############################################################
# Annotating CpGs to SNPs
# Load annotation
cpg_snp <- read_tsv(
  "../../LLS/Shared_Data/Manifests/EPIC.hg38.commonsnp.tsv")
cpg_snp <- cpg_snp %>% 
  dplyr::select(cpg=probeID, snpID)

# Keep only CpGs that map to one SNP
cpg_snp <- cpg_snp %>% 
  group_by(cpg) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n == 1) %>% 
  dplyr::select(-n)

# Merge with omicsPrint genotypes
dnam_blood <- dnam_blood %>% rownames_to_column(var = 'cpg')
dnam_blood <- inner_join(cpg_snp, dnam_blood, by = 'cpg')

dnam_muscle <- dnam_muscle %>% rownames_to_column(var = 'cpg')
dnam_muscle <- inner_join(cpg_snp, dnam_muscle, by = 'cpg')

dnam_fat <- dnam_fat %>% rownames_to_column(var = 'cpg')
dnam_fat <- inner_join(cpg_snp, dnam_fat, by = 'cpg')

##############################################################
# Subsetting
# Format rownames
rownames(lls_plink$genotypes) <- gsub('.', '_', rownames(lls_plink$genotypes), 
                                      fixed=TRUE)

# Subset keys for LLS samples
key_blood <- key_blood %>% filter(GWASnr %in% rownames(lls_plink$genotypes))
key_muscle <- key_muscle %>% filter(GWASnr %in% rownames(lls_plink$genotypes))
key_fat <- key_fat %>% filter(GWASnr %in% rownames(lls_plink$genotypes))

# Subset omicsPrint genotype samples
dnam_blood <- dnam_blood %>% dplyr::select(cpg, snpID, key_blood$sample_ID)
dnam_muscle <- dnam_muscle %>% dplyr::select(cpg, snpID, key_muscle$sample_ID)
dnam_fat <- dnam_fat %>% dplyr::select(cpg, snpID, key_fat$sample_ID)

# Subsetting genotype chip data samples
lls_plink <- as.data.frame(t(lls_plink$genotypes))

genotype_blood <- lls_plink %>% 
  dplyr::select(key_blood$GWASnr) %>% 
  rownames_to_column(var = 'snpID') %>% 
  filter(snpID %in% dnam_blood$snpID) %>% 
  arrange(snpID)

genotype_muscle <- lls_plink %>% 
  dplyr::select(key_muscle$GWASnr) %>% 
  rownames_to_column(var = 'snpID') %>% 
  filter(snpID %in% dnam_muscle$snpID) %>% 
  arrange(snpID)

genotype_fat <- lls_plink %>% 
  dplyr::select(key_fat$GWASnr) %>% 
  rownames_to_column(var = 'snpID') %>% 
  filter(snpID %in% dnam_fat$snpID) %>% 
  arrange(snpID)

# Subsetting SNPs
dnam_blood <- dnam_blood %>% 
  filter(snpID %in% genotype_blood$snpID) %>% 
  arrange(snpID)

dnam_muscle <- dnam_muscle %>% 
  filter(snpID %in% genotype_muscle$snpID) %>% 
  arrange(snpID)

dnam_fat <- dnam_fat %>% 
  filter(snpID %in% genotype_fat$snpID) %>% 
  arrange(snpID)

# Make numeric
genotype_blood[, c(2:ncol(genotype_blood))] <- sapply(genotype_blood[, c(2:ncol(genotype_blood))], as.numeric)
genotype_muscle[, c(2:ncol(genotype_muscle))] <- sapply(genotype_muscle[, c(2:ncol(genotype_muscle))], as.numeric)
genotype_fat[, c(2:ncol(genotype_fat))] <- sapply(genotype_fat[, c(2:ncol(genotype_fat))], as.numeric)

##############################################################
# Genotype Check Analysis
# Blood - Initialize
dnam <- dnam_blood
gwas <- genotype_blood
key <- key_blood
out_blood <- data.frame()

# Table of assumed against actual match and number of SNPs
for(j in 3:ncol(dnam)){
  ind <- dnam[,j]
  for(i in 2:ncol(gwas)){
    out_blood <- rbind(out_blood, 
                       data.frame(sample_ID = (colnames(dnam)[j]),
                                  actual_match = colnames(gwas)[i], 
                                  assumed_match = key[key$sample_ID == colnames(dnam)[j],6],
                                  n_snps = data.frame(xtabs(~gwas[,i] == ind))[2,2]))
  }
}
out <- out_blood %>% 
  group_by(sample_ID) %>% 
  mutate(max = max(n_snps)) %>% 
  ungroup()

# Plot results
png('../GOTO_Data/Processing/Genotypes/Figures/GOTO_figure-genotypeCheck-blood.png', width = 1000, height = 300)
out %>% ggplot(aes(x = sample_ID, y = n_snps)) +
  geom_point(fill = ifelse(out$actual_match == out$assumed_match, 'yellow', 'grey60'),
             color = ifelse(out$n_snps == out$max, 'red', 'grey60'), shape = 21, size = 2) +
  theme_bw() + ylab('Number of matching SNPs \n') + xlab('Sample ID') + ggtitle('Comparison of genotypes for all sample-person pairs') +
  theme(axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = 'white'),
        axis.text.x.bottom = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        title = element_text(size = 20),
        axis.text.y = element_text(size = 14)) +
  scale_x_discrete(expand = c(-1.01,-1.01))
dev.off()

# Inspect
out %>% 
  filter(n_snps == max) %>% 
  arrange(max)
out %>% 
  filter(n_snps == max) %>% 
  filter(actual_match != assumed_match)

# Save
save(out, key_blood, file='../GOTO_Data/Processing/Genotypes/GOTO_genotypeCheck-LLS-blood.Rdata')

##############################################################
# Genotype Check Analysis
# Muscle - Initialize
dnam <- dnam_muscle
gwas <- genotype_muscle
key <- key_muscle
out_muscle <- data.frame()

# Table of assumed against actual match and number of SNPs
for(j in 3:ncol(dnam)){
  ind <- dnam[,j]
  for(i in 2:ncol(gwas)){
    out_muscle <- rbind(out_muscle, 
                        data.frame(sample_ID = (colnames(dnam)[j]),
                                   actual_match = colnames(gwas)[i], 
                                   assumed_match = key[key$sample_ID == colnames(dnam)[j],6],
                                   n_snps = data.frame(xtabs(~gwas[,i] == ind))[2,2]))
  }
}
out <- out_muscle %>% 
  group_by(sample_ID) %>% 
  mutate(max = max(n_snps)) %>% 
  ungroup()

# Plot
png('../GOTO_Data/Processing/Genotypes/Figures/GOTO_figure-genotypeCheck-muscle.png', width = 1000, height = 300)
out %>% ggplot(aes(x = sample_ID, y = n_snps)) +
  geom_point(fill = ifelse(out$actual_match == out$assumed_match, 'yellow', 'grey60'),
             color = ifelse(out$n_snps == out$max, 'red', 'grey60'), shape = 21, size = 2) +
  theme_bw() + ylab('Number of matching SNPs \n') + xlab('Sample ID') + ggtitle('Comparison of genotypes for all sample-person pairs') +
  theme(axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = 'white'),
        axis.text.x.bottom = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        title = element_text(size = 20),
        axis.text.y = element_text(size = 14)) +
  scale_x_discrete(expand = c(-1.01,-1.01))
dev.off()

# Inspect
out %>% filter(n_snps == max) %>% arrange(max)
out %>% filter(n_snps == max) %>% filter(actual_match != assumed_match)

# 203548970088_R08C01 is actually from IOP2 ID 60365 not 60366
# 203548980011_R03C01 is actually from IOP2 ID 60366 not 60365

# Save
key_muscle %>% filter(GWASnr %in% c('58_1_1c5', '58_1_1c5s'))
save(out, key_muscle, file='../GOTO_Data/Processing/Genotypes/GOTO_genotypeCheck-LLS-muscle.Rdata')

##############################################################
# Genotype Check Analysis
# Fat - Initialize
dnam <- dnam_fat
gwas <- genotype_fat
key <- key_fat
out_fat <- data.frame()

# Table of assumed against actual match and number of SNPs
for(j in 3:ncol(dnam)){
  ind <- dnam[,j]
  for(i in 2:ncol(gwas)){
    out_fat <- rbind(out_fat, 
                     data.frame(sample_ID = (colnames(dnam)[j]),
                                actual_match = colnames(gwas)[i], 
                                assumed_match = key[key$sample_ID == colnames(dnam)[j],6],
                                n_snps = data.frame(xtabs(~gwas[,i] == ind))[2,2]))
  }
}
out <- out_fat %>% group_by(sample_ID) %>% mutate(max = max(n_snps)) %>% ungroup()

# Plot
png('../GOTO_Data/Processing/Genotypes/Figures/GOTO_figure-genotypeCheck-fat.png', width = 1000, height = 300)
out %>% ggplot(aes(x = sample_ID, y = n_snps)) +
  geom_point(fill = ifelse(out$actual_match == out$assumed_match, 'yellow', 'grey60'),
             color = ifelse(out$n_snps == out$max, 'red', 'grey60'), shape = 21, size = 2) +
  theme_bw() + ylab('Number of matching SNPs \n') + xlab('Sample ID') + ggtitle('Comparison of genotypes for all sample-person pairs') +
  theme(axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = 'white'),
        axis.text.x.bottom = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        title = element_text(size = 20),
        axis.text.y = element_text(size = 14)) +
  scale_x_discrete(expand = c(-1.01,-1.01))
dev.off()

# Inspect
out %>% filter(n_snps == max) %>% arrange(max)
out %>% filter(n_snps == max) %>% filter(actual_match != assumed_match)

# Save
save(out, key_fat, file='../GOTO_Data/Processing/Genotypes/GOTO_genotypeCheck-LLS-fat.Rdata')

##############################################################
