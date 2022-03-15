#Traxinger, Brianna R. et al 2021
#Scripts to reproduce figures in manunscript, starting from CellRanger outputs
#One script with 1 dataset
#Using Seurat v3: Stuart*, Butler*, et al., Cell 2019 
#Raw and processed data files have been deposited in the NCBI’s Gene Expression Omnibus and are accessible through GEO series accession number PENDING.


##-- Load or Install packages, including MAST

# https://bioconductor.org/install/

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")

BiocManager::install(c("MAST", "RColorBrewer", "Seurat", "tidyverse", "circlize", 'multtest'))
library(Seurat)
library(tidyverse)
library(MAST)
library(RColorBrewer)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
install.packages('metap')
library(metap)


##-- Load the Demo PBMC dataset, required format is the unzipped file from the CellRanger workflow

LN_HSV2_CD4.data <- Read10X(data.dir = 'LN_HSV2_CD4')
LN_HSV2_Treg.data <- Read10X(data.dir = 'LN_HSV2_Treg')
LN_uninfected_CD4.data <- Read10X(data.dir = 'LN_uninfected_CD4')
LN_uninfected_Treg.data <- Read10X(data.dir = 'LN_uninfected_Treg')
VT_HSV2_CD4.data <- Read10X(data.dir = 'VT_HSV2_CD4')
VT_HSV2_Treg.data <- Read10X(data.dir = 'VT_HSV2_Treg')
VT_uninfected_CD4.data <- Read10X(data.dir = 'VT_uninfected_CD4')


#- Setup seurat object class

LN_HSV2_CD4_Seurat <- CreateSeuratObject(counts = LN_HSV2_CD4.data, 
                                         assay = 'RNA',
                                         min.cells = 3,
                                         min.features = 200,
                                         project = 'LN_HSV2_CD4')
LN_HSV2_CD4_Seurat #13132 features across 5615 samples within 1 assay  

LN_HSV2_Treg_Seurat <- CreateSeuratObject(counts = LN_HSV2_Treg.data, 
                                          assay = 'RNA',
                                          min.cells = 3,
                                          min.features = 200,
                                          project = 'LN_HSV2_Treg')
LN_HSV2_Treg_Seurat  #11856 features across 1561 samples within 1 assay 

LN_uninfected_CD4_Seurat <- CreateSeuratObject(counts = LN_uninfected_CD4.data, 
                                               assay = 'RNA',
                                               min.cells = 3,
                                               min.features = 200,
                                               project = 'LN_uninfected_CD4')
LN_uninfected_CD4_Seurat #11651 features across 2948 samples within 1 assay 


LN_uninfected_Treg_Seurat <- CreateSeuratObject(counts = LN_uninfected_Treg.data, 
                                                assay = 'RNA',
                                                min.cells = 3,
                                                min.features = 200,
                                                project = 'LN_uninfected_Treg')
LN_uninfected_Treg_Seurat #12232 features across 2184 samples within 1 assay 

VT_HSV2_CD4_Seurat <- CreateSeuratObject(counts = VT_HSV2_CD4.data, 
                                         assay = 'RNA',
                                         min.cells = 3,
                                         min.features = 200,
                                         project = 'VT_HSV2_CD4')
VT_HSV2_CD4_Seurat #13833 features across 6034 samples within 1 assay 

VT_HSV2_Treg_Seurat <- CreateSeuratObject(counts = VT_HSV2_Treg.data, 
                                          assay = 'RNA',
                                          min.cells = 3,
                                          min.features = 200,
                                          project = 'VT_HSV2_Treg')
VT_HSV2_Treg_Seurat #9689 features across 357 samples within 1 assay 
VT_uninfected_CD4_Seurat <- CreateSeuratObject(counts = VT_uninfected_CD4.data, 
                                               assay = 'RNA',
                                               min.cells = 3,
                                               min.features = 200,
                                               project = 'VT_uninfected_CD4')
VT_uninfected_CD4_Seurat #8785 features across 258 samples within 1 assay 


#merge

MergedSeurat <- merge(LN_HSV2_CD4_Seurat, y = c(LN_HSV2_Treg_Seurat, LN_uninfected_CD4_Seurat, LN_uninfected_Treg_Seurat, VT_HSV2_CD4_Seurat, VT_HSV2_Treg_Seurat, VT_uninfected_CD4_Seurat), add.cell.ids = c("LN_HSV2_CD4", "LN_HSV2_Treg", "LN_uninfected_CD4", "LN_uninfected_Treg", "VT_HSV2_CD4", "VT_HSV2_Treg", "VT_uninfected_CD4"), project = "Ms_HSV2")

table(MergedSeurat@meta.data$orig.ident)

##-- QC metrics - Percentage of mitochondrial genes
#- percent.mito variable

mito.genes <- grep(pattern = '^mt-', x = rownames(x = MergedSeurat@assays$RNA@data), value = TRUE) # 13 mitochondrial genes
percent.mito <- Matrix::colSums(MergedSeurat@assays$RNA@counts[mito.genes, ]) / Matrix::colSums(MergedSeurat@assays$RNA@counts) # Matrix package is needed


#- Use of AddMetaData to add column to object@meta.data 

MergedSeurat <- AddMetaData(object = MergedSeurat, 
                            metadata = percent.mito, 
                            col.name = 'percent.mito')


##-- QC metrics - draw a VlnPlot()
#Rna counts per sample

VlnPlot(object = MergedSeurat, 
        features = c('percent.mito', 'nCount_RNA', 'nFeature_RNA'), 
        ncol = 3, pt.size = 0.0) + stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95)
?stat_summary

#-- QC metrics - another QC metric using GenePlot(), correlation displayed on top

par(mfrow = c(1, 2))
FeatureScatter(object = MergedSeurat, feature1 = 'nCount_RNA', feature2 = 'percent.mito')
FeatureScatter(object = MergedSeurat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')

##-- Filtering (removing the outlier cells with percent.mito being higher than 10%
#set minimun RNA count to 500

MergedSeurat2 <- subset(MergedSeurat, subset = percent.mito < 0.05 & nFeature_RNA > 400 & nFeature_RNA < 5000 & nCount_RNA > 400 & nCount_RNA < 35000)

table(MergedSeurat2@meta.data$orig.ident)

#LN_HSV2_CD4       LN_HSV2_Treg  LN_uninfected_CD4 LN_uninfected_Treg 
#5404               1416               2723               1928 
#VT_HSV2_CD4       VT_HSV2_Treg  VT_uninfected_CD4 
#5443                264                241 
# -Inf and Inf should be used if you don't want a lower or upper threshold

#MergedSeurat2 #14404 features across 17419 samples

#Re-plot after filtering

par(mfrow = c(1, 2))
FeatureScatter(object = MergedSeurat2, feature1 = 'nCount_RNA', feature2 = 'percent.mito')
FeatureScatter(object = MergedSeurat2, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')

##-- Normalization, see R Markdownsheet

MergedSeurat3 <- NormalizeData(object = MergedSeurat2, 
                               normalization.method = 'LogNormalize', 
                               scale.factor = 10000)

##-- Highly variable genes, y is the dispersion/variance, cutoffs need to be set based on dataset 

MergedSeurat4 <- FindVariableFeatures(object = MergedSeurat3, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes

top10 <- head(VariableFeatures(MergedSeurat4), 10)


#identify top 20 highly variable genes

top20 <- head(VariableFeatures(MergedSeurat4), 20)

# plot variable features with and without labels

plot1 <- VariableFeaturePlot(MergedSeurat4)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
patchwork::wrap_plots(plot1, plot2, ncol = 2)


##-- Scaling the data and removing unwanted sources of variation - this function takes about 1-2 minutes, might be longer
#scaling with all samples present

MergedSeurat5 <- ScaleData(object = MergedSeurat4, 
                           model.use = 'linear',
                           vars.to.regress = c('nCount_RNA', 'percent.mito')) # nUMI is used as proxy for Cellular Detection Rate
#-- add batch if you have different batches of data. Usually regressing for number of UMI sufficient

#-- perform PCA using the scale-corrected matrix, pcs.compute is the number of principal components
MergedSeurat6 <- RunPCA(object = MergedSeurat5, 
                        pc.genes = MergedSeurat5@assays$RNA@var.features, 
                        pcs.compute = 50,
                        do.print = FALSE)

##-- PCA - PrintPCA

print(MergedSeurat6[["pca"]], dims = 1:5, nfeatures = 5) # genes.print is number of genes to print for each PC

#-- PCA - VizPCA, plot top 30 genes for PC1 and PC2, 
#-- check whether one of the PCs is driven my mitochondrial genes or cell cycle-genes

VizDimLoadings(MergedSeurat6, dims = 1:2, reduction = "pca", nfeatures = 30)

#-- PCA - PCAPlot

PCAPlot(object = MergedSeurat6, 
        dims = c(1, 2))

#-- PCElbowPlot, plots the standard deviation for each PC, look at the "elbow" of the curve (JackStraw procedure would take much longer)

ElbowPlot(object = MergedSeurat6,
          ndims = 50)


##-- tSNE, use 10 of the PCs based on the above observations, do.fast will define whether to use Barnes-Hut implementation, check for perplexity 30-50

AllSamples <- RunTSNE(object = MergedSeurat6, 
                      dims.use = 1:30,
                      perplexity = 50,
                      do.fast = TRUE,
                      seed.use = 34) 


##-- UMAP
#all samples, default clustering

AllSamples1<- SetIdent(AllSamples, value = "seurat_clusters")
UMAPall <- FindNeighbors(object = AllSamples, reduction = 'pca', dims = 1:30, k.param = 30, force.recalc = T)
UMAPall <- FindClusters(object = UMAPall, verbose = TRUE, n.start = 100)
UMAPall <- RunUMAP(UMAPall, reduction.use = "pca", dims = 1:30, seed.use = 34)
DimPlot(UMAPall, reduction = "umap", pt.size = 0.6, label = T)

# UMAP by Ident

UMAPallIDENT <- SetIdent(UMAPall, value = "orig.ident")
DimPlot(UMAPallIDENT, reduction = "umap", pt.size = 0.6)
table(Idents(UMAPall))

# of cells per cluster

table(Idents(UMAPall2))

# number of RNA counts per cluster

VlnPlot(UMAPall2, features = c("nFeature_RNA"), pt.size = 0.0)  + stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95)

##--Make new Seurat objects to re-run UMAP on specific groups

Treg <- SetIdent(AllSamples, value = "orig.ident")
Treg<- subset(Treg, idents = c("LN_uninfected_Treg", "LN_HSV2_Treg", "VT_HSV2_Treg"))

##--rescale Treg with all genes

all.genes <- rownames(Treg)
Treg2 <- ScaleData(Treg, features = all.genes, vars.to.regress = c('nCount_RNA', 'percent.mito'))

Treg2 <- SetIdent(Treg2, value = "orig.ident")


#Treg UMAP

TregUMAP<- SetIdent(Treg2, value = "seurat_clusters")
TregUMAP <- FindNeighbors(object = Treg2, reduction = 'pca', dims = 1:30, k.param = 30, force.recalc = T)
TregUMAP <- FindClusters(object = TregUMAP, verbose = TRUE, n.start = 100)
TregUMAP <- RunUMAP(TregUMAP, reduction.use = "pca", dims = 1:30, seed.use = 34)
DimPlot(TregUMAP, reduction = "umap", pt.size = 0.6, label = T)

#by Ident

TregUMAPident <- SetIdent(TregUMAP, value = "orig.ident")
DimPlot(TregUMAPident, reduction = "umap", pt.size = 0.6, cols = c("grey", "green", "blue"))

#Treg UMAP by cluster, res = 0.4

TregUMAPcluster.4 <- FindNeighbors(object = Tregcluster, reduction = 'pca', dims = 1:30, k.param = 30, force.recalc = T)
TregUMAPcluster.4 <- FindClusters(object = TregUMAPcluster.4, verbose = TRUE, n.start = 100, resolution = 0.4, force.recalc= T)
TregUMAPcluster.4 <- RunUMAP(TregUMAPcluster.4, reduction.use = "pca", dims = 1:30, seed.use = 34)
DimPlot(TregUMAPcluster.4, reduction = "umap", pt.size = 0.1, label = T)

# Number of total cells per cluster
table(Idents(TregUMAPcluster.4))

# of cells from each Treg ident total
table(TregUMAPcluster.4@meta.data$orig.ident)

# of cells from each ident per cluster
table(Idents(TregUMAPcluster.4), TregUMAPcluster.4$orig.ident)

#Violin plots

VlnPlot(Treg2, features = c("Gzmb"), group.by = "orig.ident", 
        pt.size = .5, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

VlnPlot(Treg2, features = c("Nkg7"), group.by = "orig.ident", 
        pt.size = .5, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)


##--Find all markers

#find all markers, Treg, by ident

All.markers_Treg_ident <- FindAllMarkers(object = Treg2, 
                                         only.pos = F,
                                         min.pct = 0.25, # defines minimum percent of cells expressing a given gene
                                         logfc.threshold = 0.5, # log-FC
                                         return.thresh = 0.01, # Only return markers that have a p-value < 0.01
                                         test.use = 'MAST',
                                         latent.vars = 'nCount_RNA') # nUMI as proxy for CDR

#Find all markers, Treg, by cluster
All.markers_Treg_clusters <- FindAllMarkers(object = TregUMAPcluster.4, 
                                            only.pos = F,
                                            min.pct = 0.25, # defines minimum percent of cells expressing a given gene
                                            logfc.threshold = 0.5, # log-FC
                                            return.thresh = 0.01, # Only return markers that have a p-value < 0.01
                                            test.use = 'MAST',
                                            latent.vars = 'nCount_RNA') # nUMI as proxy for CDR

write.csv(All.markers_Treg_clusters, file = "All_markersTreg_clusters_res0.4_9Mar22.csv")


##--Heatmaps

#Treg heatmap by ident, top 20
top20Allident <- All.markers_Treg_ident %>% 
  group_by(cluster) %>% 
  top_n(20, avg_log2FC)

#Change order

levels(Treg2) <- c("LN_HSV2_Treg", "VT_HSV2_Treg", "LN_uninfected_Treg")

TregHeatmapident20 <- DoHeatmap(object = Treg2, 
                              features = top20Allident$gene)

TregHeatmapident20

#Treg heat map by cluster, top 20
#top20cluster
top20Allcluster <- All.markers_Treg_clusters %>% 
  group_by(cluster) %>% 
  top_n(20, avg_log2FC)

TregHeatmapcluster <- DoHeatmap(object = TregUMAPcluster.4, 
                                features = top20Allcluster$gene)



