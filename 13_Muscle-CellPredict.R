##############################################################
# SCRIPT 13: Muscle - Predict cell proportions with MuSiC
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(Biobase)
library(MuSiC)
library(SummarizedExperiment)
library(GEOquery)

##############################################################
# snRNA counts
# Load
load('../GOTO_Data/Cell_Counts/Muscle/snCounts.Rdata')
dim(counts)

# Replace NA with 0
counts[is.na(counts)] <- 0

# Save gene names
gene_names <- rownames(counts)

##############################################################
# Bulk RNAseq
# Load functions
source("../GOTO_Data/GOTO_RNA/goto.rnaseq.functions.R")

# Paths
pathIN_dat <- "../GOTO_Data/GOTO_RNA/merge.gene.counts_biopet_13052016.RData"
pathIN_cov <- "../GOTO_Data/GOTO_RNA/muscle_QC_covariates_filesv2.csv"

# Muscle complete pairs
filt.samp <- "tissue_muscle|qc_sexswitch|qc_multdim2|qc_rep1|complete_pairs"

# Read in
dat1 <- read.gotornaseq(pathIN_dat = pathIN_dat, 
                        pathIN_cov, 
                        filt.samp = filt.samp, 
                        type="counts",
                        quiet = FALSE)

# Format
goto_counts <- dat1[[1]][,grep('ENS',colnames(dat1[[1]]))]
goto_counts <- as.data.frame(t(goto_counts))
colnames(goto_counts) <- str_c(dat1[[1]]$sampID2,
                               "_",
                               dat1[[1]]$visitnr)

# Map ensembl ID to name
ens2gene <- cinaR::grch38
m <- match(rownames(goto_counts), ens2gene$ensgene)
mapped.genes <- ens2gene$symbol[m]

##############################################################
# Subsetting
# Remove genes
removed.genes <- duplicated(mapped.genes) | is.na(mapped.genes) | grepl("^MT", mapped.genes)
goto_counts <- goto_counts[!removed.genes,]
rownames(goto_counts) <- mapped.genes[!removed.genes]

# Subset snRNA genes
goto_counts <- goto_counts[gene_names,]
goto_counts <- goto_counts[!rowSums(is.na(goto_counts)) == 162,]
goto_counts[1:5, 1:5]

# Save goto genes
goto_genes <- rownames(goto_counts)

# Save overlapping genes
bin <- gene_names %in% goto_genes
xtabs(~bin)
new_names <- gene_names[bin]
counts <- counts[bin,]
goto_counts <- goto_counts[new_names,]

##############################################################
# Annotate to Clusters
# Load clusters
load('../GOTO_Data/Cell_Counts/Muscle/snClusters.Rdata')

cluster_key <- data.frame(
  seurat_clusters = c(0,1,2,3,4,5,6,7),
  cell_type = c('Fast Skeletal Muscle',
                'Slow Skeletal Muscle',
                'FAPs',
                'Endothelial cell',
                'T/NK/Macrophages',
                'SATs',
                'Smooth Muscle Cell',
                'Slow Skeletal Muscle')
)

nuclei_clusters <- nuclei_clusters %>% 
  mutate(
    cell_type = case_when(
      seurat_clusters == 0 ~ "Fast Skeletal Muscle",
      seurat_clusters == 1 ~ "Slow Skeletal Muscle",
      seurat_clusters == 2 ~ "FAPs",
      seurat_clusters == 3 ~ "Endothelial Cells",
      seurat_clusters == 4 ~ "T/NK/Macrophages",
      seurat_clusters == 5 ~ "SATs",
      seurat_clusters == 6 ~ "Smooth Muscle",
      seurat_clusters == 7 ~ 'Slow Skeletal Muscle'
    )
  )

nuclei_clusters <- nuclei_clusters %>% 
  dplyr::select(cell_type)

head(nuclei_clusters)

cell_types <- nuclei_clusters[colnames(counts),]
cell_types[1:5]

# Transform and label samples
meta <- as.data.frame(cell_types)
colnames(meta) <-'Cell' 
meta$Sample <- colnames(counts)
rownames(meta) <- colnames(counts)

##############################################################
# Get data ready
# Sets
counts <- as.data.frame(summary(counts))

# Remove non-zero genes
nz.gene = rownames(counts)[(rowSums(counts) != 0)]
counts <- counts[nz.gene,]

# Convert counts to matrix
dGC_to_matrix <- function(x,ncol_break = 49999){
  if(length(colnames(x))>(ncol_break+1)){
    total_cols = length(colnames(x)) 
    the_seq <- c(seq(1,total_cols,ncol_break), total_cols)  
    the_seq <- unique(the_seq) 
  }
  matrix_list <- list() 
  total_parts <- length(the_seq)-1 
  for(i in 1:total_parts){
    start_no = ifelse(i==1,1,the_seq[i]+1)  
    print(paste0(i, " is i"))
    print(paste0("start_no is", start_no))
    end_no = the_seq[i+1] 
    print(paste0("part_number:", i, ";cols-",start_no,":",end_no))
    matrix_list[[i]] <- as.matrix(x[,start_no:end_no,drop = F])
  }
  return(do.call(cbind, matrix_list)) 
}
full_mtx <- dGC_to_matrix(counts, 50000)

# Create expression sets
C.eset <- Biobase::ExpressionSet(
  assayData = full_mtx, 
  phenoData = Biobase::AnnotatedDataFrame(meta))
C.eset

T.eset <- Biobase::ExpressionSet(assayData = as.matrix(goto_counts))
T.eset

##############################################################
# Deconvolute!
deconv <- music_prop(
  bulk.eset = T.eset, 
  sc.eset = C.eset, 
  clusters = 'Cell',
  markers = NULL, 
  normalize = FALSE, 
  samples = 'Sample', 
  verbose = T)$Est.prop.weighted

# Summarize
names <- colnames(deconv)
summary(deconv)
head(deconv)
heatmap(deconv, margins=c(12,8))

# Save
save(deconv, file='../GOTO_Data/Cell_Counts/Muscle/music_deconv.Rdata')

##############################################################
# Add to data
# Load
load('../GOTO_Data/GOTO_targets-filtered.Rdata')
load('../GOTO_Data/GOTO_methData-filtered.Rdata')

# Format
targets <- targets %>% 
  mutate(join = paste0(IOP2_ID, '_', 
                       ifelse(targets$timepoint=="after", 2, 1)))


deconv <- as.data.frame(deconv) %>% 
  rownames_to_column(var="join")
colnames(deconv) <- c("join", "cc_musc_slowSM", "cc_musc_FAP", 
                      "cc_musc_fastSM", "cc_musc_SAT", 
                      "cc_musc_TNKmacro", "cc_musc_endo", 
                      "cc_musc_smooth")

# Merge
targets <- left_join(targets, deconv, by = "join")

# Save
colData(methData) <- DataFrame(targets)
save(targets, file='../GOTO_Data/GOTO_targets-filtered.Rdata')
save(methData, file="../GOTO_Data/GOTO_methData-filtered.Rdata")

##############################################################

