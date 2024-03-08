# script to annotate cell types from 20k Human PBMCs from a healthy female donor
# setwd("~/Desktop/demo/singleCell_singleR/scripts")

BiocManager::install("SingleR")


library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)


#setup directory
setwd("E:/2023-09-27-Moslehi-scrnaseq/data/GSE227718/GSM7106037_sc3.TPM")


#Reading csv file
GSM4955731_R1_N <- as.sparse(read.csv(file = "GSM7106037_sc3.TPM.csv", sep = ",",
                                      header = TRUE, row.names = 1))


# Initialize the Seurat object with the raw (non-normalized data).
GBM_4 <- CreateSeuratObject(counts = GSM4955731_R1_N, project = "GBM", min.cells = 3, min.features = 200)
GBM_4




# Input Data 10X CellRanger .HDF5 format --------------
#hdf5_obj <- Read10X_h5(filename = '../data/20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5',
#                      use.names = TRUE,
#                     unique.features = TRUE)



#pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj)

# QC and Filtering -----------
# explore QC
GBM_4$mitoPercent <- PercentageFeatureSet(GBM_4, pattern = '^MT-')
GBM_4.filtered <- subset(GBM_4, subset = nCount_RNA > 800 &
                           nFeature_RNA > 500 &
                           mitoPercent < 10)


# It is a good practice to filter out cells with non-sufficient genes identified and genes with non-sufficient expression across cells.


# pre-process standard workflow ---------------
GBM_4.filtered <- NormalizeData(object = GBM_4.filtered)
GBM_4.filtered <- FindVariableFeatures(object = GBM_4.filtered)
GBM_4.filtered <- ScaleData(object = GBM_4.filtered)
GBM_4.filtered <- RunPCA(object = GBM_4.filtered)
GBM_4.filtered <- FindNeighbors(object = GBM_4.filtered, dims = 1:20)
GBM_4.filtered <- FindClusters(object = GBM_4.filtered)
GBM_4.filtered <- RunUMAP(object = GBM_4.filtered, dims = 1:20)

# running steps above to get clusters
View(GBM_4.filtered@meta.data)

jpeg("dimplot.jpg", width = 500, height = 500)
DimPlot(GBM_4.filtered, reduction = 'umap')
dev.off()

# get reference data -----------
ref <- celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))

# expression values are log counts (log normalized counts)


# run SingleR (default mode) ---------
# default for SingleR is to perform annotation of each individual cell in the test dataset

GBM_4.assayed <- GetAssayData(GBM_4.filtered)
pred <- SingleR(test = GBM_4.assayed,
                ref = ref,
                labels = ref$label.main)

pred
pred$labels
View(as.data.frame(pred))

#GBM_4.assayed$singleR.labels <- pred$labels[match(rownames(GBM_4.filtered@meta.data), rownames(pred))]

#jpeg("sigleR_lables.jpg", width = 500, height = 500)
#DimPlot(GBM_4.assayed, reduction = 'umap', group.by = 'singleR.labels')
#dev.off()

# Annotation diagnostics ----------


# ...Based on the scores within cells -----------
pred
pred$scores


jpeg("heatmap_clusters.jpg", width = 500, height = 500)
plotScoreHeatmap(pred)
dev.off()

# ...Based on deltas across cells ----------
jpeg("distribution_clusters.jpg", width = 500, height = 500)
plotDeltaDistribution(pred)
dev.off()



# ...Comparing to unsupervised clustering ------------

comparison <- table(Assigned=pred$labels, Clusters=GBM_4.filtered$seurat_clusters)
write.csv(comparison, "comparison.csv")

jpeg("Vlnplot2.jpg", width = 500, height = 500)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))
dev.off()









