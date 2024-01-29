##############################################################
# SCRIPT 12: Muscle - Cluster scRNA using Seurat
##############################################################
# Setup
rm(list=ls())
library(gridExtra)
library(tidyverse)
library(Seurat)
library(plotly)
library(gplots)

# Initialize
outputDir = getwd()
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

##############################################################
# Load count data
# Make function
SeuratObjectFromFile <- function(file_path, project_id) {
  counts_table <- read.csv(file=file_path)
  rownames(counts_table) <- counts_table[,1]
  counts_table[,1] <- NULL
  seur_obj <- CreateSeuratObject(counts=counts_table, 
                                 project=project_id, 
                                 min.cells=3, 
                                 min.features=200)
  return(seur_obj)
}

# Load all datasets
HM1 <- SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098737_HM1.csv",
  project_id="HM1")
HM2 <- SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098738_HM2.csv",
  project_id="HM2")
HM3 <- SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098739_HM3.csv",
  project_id="HM3")
HM4 <- SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098740_HM4.csv",
  project_id="HM4")
HM5 <-SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098741_HM5.csv",
  project_id="HM5")
HM6 <- SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098742_HM6.csv",
  project_id="HM6")
HM7 <- SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098743_HM7.csv",
  project_id="HM7")
HM8 <- SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098744_HM8.csv",
  project_id="HM8")
HM9 <- SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098745_HM9.csv",
  project_id="HM9")
HM10 <- SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098746_HM10.csv",
  project_id="HM10")
HM11 <- SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098747_HM11.csv",
  project_id="HM11")
HM12 <- SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098748_HM12.csv",
  project_id="HM12")
HM13 <- SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098749_HM13.csv",
  project_id="HM13")
HM15 <- SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098750_HM15.csv",
  project_id="HM15")
HM21 <- SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098751_HM21.csv",
  project_id="HM21")
HM22 <- SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098752_HM22.csv",
  project_id="HM22")
HM23 <- SeuratObjectFromFile(
  file_path="../GOTO_Data/Cell_Counts/Muscle/GSE167186/out/GSM5098753_HM23.csv",
  project_id="HM23")

# Merge
object <- merge(x=HM1,
                y=c(HM2, HM3, HM4, HM5, HM6, HM7, HM8, HM9, HM10, HM11, HM12, HM13, HM15, HM21, HM22, HM23),
                add.cell.ids=c("HM1", "HM2", "HM3", "HM4", "HM5", "HM6", "HM7", "HM8", "HM9", "HM10", "HM11", "HM12", "HM13", "HM15", "HM21", "HM22", "HM23"))

##############################################################
# Quality control
# Remove cells with < 500 genes or > 10% MT  
object[["percent.mt"]] <- PercentageFeatureSet(object, 
                                               pattern = "^MT-")
cc_genes_in_dataset <- c(s.genes, g2m.genes)[c(s.genes, g2m.genes) %in% rownames(object)]
object[["percent.cc_genes"]] <- PercentageFeatureSet(object, 
                                                     features = cc_genes_in_dataset)

# Plot
VlnPlot(object, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.1) 
FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
object <- subset(object, subset = nFeature_RNA > 200 & 
                   percent.mt < 10)
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.1)
object

##############################################################
# Normalization
object <- NormalizeData(object, 
                        normalization.method = "LogNormalize", 
                        scale.factor = 10000)

# Variable genes & Scaling
object <- FindVariableFeatures(object, 
                               selection.method = "vst", 
                               nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(object), 10)
plot1 <- VariableFeaturePlot(object)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
all.genes <- rownames(object)
object <- ScaleData(object, features = all.genes)

##############################################################
# Linear Dimensional Reduction
object <- RunPCA(object, 
                 features = VariableFeatures(object = object))
print(object[["pca"]], 
      dims = 1:5, 
      nfeatures = 8)
VizDimLoadings(object, 
               dims = 1:2, 
               reduction = "pca")

# Plot
DimPlot(object, 
        reduction = "pca")
DimHeatmap(object, 
           dims = 1, 
           cells = 500, 
           balanced = TRUE)
ElbowPlot(object)

# Clustering
object <- FindNeighbors(object, 
                        dims = 1:10)
object <- FindClusters(object, 
                       resolution = 0.14)

##############################################################
# Run non-linear dimensional reduction UMAP
object <- RunUMAP(object, 
                  dims = 1:10)
DimPlot(object, 
        reduction = "umap")
saveRDS(object, file = "../GOTO_Data/Cell_Counts/Muscle/muscle.rds")

# Cluster markers
markers <- FindAllMarkers(object,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)

allMarkersPath = "../GOTO_Data/Cell_Counts/Muscle/all_markers.tsv"
write.table(x = markers,
            file = allMarkersPath,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

topMarkers <- markers %>%
  group_by(cluster) %>%
  top_n(-30, p_val_adj)
topMarkersPath = "../GOTO_Data/Cell_Counts/Muscle/topMarkers.tsv"
write.table(x = topMarkers,
            file = topMarkersPath,
            sep = "\t", quote = FALSE,
            row.names = FALSE)

# Save Counts
nuclei_clusters <- object@meta.data
save(nuclei_clusters, file="../GOTO_Data/Cell_Counts/Muscle/snClusters.Rdata")

counts <- object@assays$RNA@counts
save(counts, file='../GOTO_Data/Cell_Counts/Muscle/snCounts.Rdata')

##############################################################