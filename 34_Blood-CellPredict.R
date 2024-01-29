##############################################################
# SCRIPT 34: Blood - Predicting cell counts
##############################################################
# Setup
rm(list=ls())
library(rlang)
library(htmltools)
library(rmarkdown)
library(cli)
library(tidyverse)
library(DNAmArray)
library(lubridate)
library(Biobase)
library(MuSiC)
library(SummarizedExperiment)
library(minfi)
library(ExperimentHub)
library(FlowSorted.Blood.EPIC)
library(FlowSorted.BloodExtended.EPIC)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

##############################################################
# Load data
load("../GOTO_Data/GOTO_targets-filtered.Rdata")
load("../GOTO_Data/GOTO_methData-filtered.Rdata")
load("../GOTO_Data/Processing/GOTO_RGset-unfiltered.Rdata")

methData
RGset

# Load measured cell counts
blood_df <- read_csv("../GOTO_Data/Cell_Counts/Blood/GOTO_Cellcounts-Medication_20210401.csv")

targets <- left_join(targets, blood_df, by=c("IOP2_ID", "timepoint"))

##############################################################
# MuSiC
# Expression Set for scRNA
sc_blood <- read.table('../GOTO_Data/Cell_Counts/Blood/scRNA-blood_GSE143704.tsv', sep = "\t")

# Save cell names
cell_names <- as.data.frame(t(sc_blood[1,-1]))
colnames(cell_names) <- "Cell"

# Remove first row with cell names
sc_blood <- as.data.frame(sc_blood[-1,])
gene_names <- sc_blood[,1]

# Remove gene symbols
sc_blood <- sc_blood[,-1]

# Expression matrix
sc_blood <- data.frame(apply(sc_blood, 2, as.numeric))

# Rownames
rownames(sc_blood) <- gene_names

# Colnames
colnames(sc_blood) <- str_pad(1:ncol(sc_blood), width=3, pad="0")

# Sample labels
cell_names$Sample <- str_pad(1:ncol(sc_blood), width=3, pad="0")
rownames(cell_names) <- str_pad(1:ncol(sc_blood), width=3, pad="0")

# GOTO Expresssion Data
source('../GOTO_Data/GOTO_RNA/goto.rnaseq.functions.R')

pathIN_dat <- "../GOTO_Data/GOTO_RNA/merge.gene.counts_biopet_13052016.RData"
pathIN_cov <- "../GOTO_Data/GOTO_RNA/datasheet_RNAseq_blood_V2.csv"

filt.samp <- "tissue_blood|qc_sexswitch|qc_multdim2|qc_rep1|complete_pairs"

goto_exp <- read.gotornaseq(pathIN_dat = pathIN_dat, pathIN_cov, filt.samp = filt.samp, quiet = FALSE)

goto_exp <- goto_exp[["dat"]]

# Save counts
goto_exp <- goto_exp %>% 
  dplyr::filter(nutridrink == 0) %>% 
  dplyr::select(sampID2, intervention, 
                starts_with('ENS')) %>% 
  mutate(ID = str_c(sampID2, '_', intervention))

# Save IDs
ID_name <- goto_exp$ID
head(ID_name)

ID_df <- data.frame(Sample = ID_name)
rownames(ID_df) = ID_name

# Remove non-gene variables
goto_exp <- goto_exp %>% dplyr::select(-ID, -sampID2, -intervention)
goto_exp <- as.data.frame(t(goto_exp))
colnames(goto_exp) <- ID_name

# Map to gene name
ens2gene <- cinaR::grch38
m <- match(rownames(goto_exp), ens2gene$ensgene)
mapped.genes <- ens2gene$symbol[m]

# Subset
removed.genes <- duplicated(mapped.genes) | is.na(mapped.genes) | grepl("^MT", mapped.genes)
goto_exp <- goto_exp[!removed.genes,]
rownames(goto_exp) <- mapped.genes[!removed.genes]

goto_exp <- goto_exp[rownames(goto_exp) %in% rownames(sc_blood),]
sc_blood <- sc_blood[rownames(sc_blood) %in% rownames(goto_exp),]

# Expression Sets
C.eset <- Biobase::ExpressionSet(
  assayData = as.matrix(sc_blood), 
  phenoData = Biobase::AnnotatedDataFrame(cell_names))
C.eset

T.eset <- Biobase::ExpressionSet(assayData = as.matrix(goto_exp),
                                 phenoData = Biobase::AnnotatedDataFrame(ID_df))
T.eset

# MuSiC
deconv <- music_prop(
  bulk.eset = T.eset, 
  sc.eset = C.eset, 
  clusters = 'Cell',
  markers = NULL, 
  normalize = FALSE, 
  samples = 'Sample', 
  verbose = F)$Est.prop.weighted

# Save
save(deconv, file="../GOTO_Data/Cell_Counts/Blood/GOTO_Blood-Music.Rdata")

# Add to targets
targets <- targets %>% 
  mutate(
    ID = paste0(IOP2_ID, "_",as.numeric(timepoint)-1)
  )

# Variable names
colnames(deconv) <- paste0("cc_blood_music_", colnames(deconv))

# ID var in deconv
deconv <- as.data.frame(deconv) %>% rownames_to_column(var="ID")

# Merge
targets <- left_join(targets, deconv, by="ID")

##############################################################
# IDOL
RGset_blood <- RGset[ , RGset$tissue == 'fasted blood']

# Minfi
hub <- ExperimentHub()  
query(hub, "FlowSorted.Blood.EPIC")  
FlowSorted.Blood.EPIC <- hub[["EH1136"]]  

# Deconvolute
idol_blood <- estimateCellCounts(
  rgSet = RGset_blood, 
  referencePlatform = 'IlluminaHumanMethylationEPIC', 
  cellTypes = c("CD8T", "CD4T", "NK", "Bcell",  
                "Mono", "Neu"),
  verbose = TRUE, meanPlot = TRUE)

# Save
save(idol_blood, 
     file="../GOTO_Data/Cell_Counts/Blood/GOTO_Blood-IDOL.Rdata")

# Join with targets
idol_blood <- as.data.frame(idol_blood) %>% 
  rownames_to_column(var = 'Basename') 

# Colnames
colnames(idol_blood) <- c("Basename", "cc_blood_idol_CD8T",
                          "cc_blood_idol_CD4T", 
                          "cc_blood_idol_NK",
                          "cc_blood_idol_Bcell",
                          "cc_blood_idol_Mono",
                          "cc_blood_idol_Neu")

targets <- left_join(targets, idol_blood, by="Basename")

##############################################################
# IDOL Extended
# Load
load("../GOTO_Data/Cell_Counts/Blood/GOTO_Blood-IDOLext.Rdata")

# Merge
idol_ext <- idol_ext %>% rownames_to_column(var="Basename")

# Colnames
colnames(idol_ext) <- c("Basename", "cc_blood_ext_Bas",
                        "cc_blood_ext_Bmem", "cc_blood_ext_Bnv",
                        "cc_blood_ext_CD4mem", "cc_blood_ext_CD4nv",
                        "cc_blood_ext_CD8mem", "cc_blood_ext_CD8nv",
                        "cc_blood_ext_Eos", "cc_blood_ext_Mono",
                        "cc_blood_ext_Neu", "cc_blood_ext_NK",
                        "cc_blood_ext_Treg")

# Merge
targets <- left_join(targets, idol_ext, by="Basename")

# Add to methData and RGset
check <- targets$Basename == colnames(methData)
xtabs(~check)

check <- targets$Basename == colnames(RGset)
xtabs(~check)

check <- colnames(RGset) == colnames(methData)
xtabs(~check)

order <- colnames(methData)

targets <- targets[match(order, targets$Basename),]
rownames(targets) <- targets$Basename

colData(methData) <- DataFrame(targets)
colData(RGset) <- DataFrame(targets)

# Save
save(targets, file="../GOTO_Data/GOTO_targets-filtered.Rdata")
save(methData, file="../GOTO_Data/GOTO_methData-filtered.Rdata")

##############################################################
