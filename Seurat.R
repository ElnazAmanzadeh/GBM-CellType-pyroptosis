Scaling the data

all.genes <- rownames(GBM_4)
GBM_4 <- ScaleData(GBM_4, features = all.genes)
#all.genes <- rownames(GBM_2)

#Perform linear dimensional reduction
GBM_4 <- RunPCA(GBM_4, features = VariableFeatures(object = GBM_4))

# Examine and visualize PCA results a few different ways
print(GBM_4[["pca"]], dims = 1:5, nfeatures = 20)

VizDimLoadings(GBM_4, dims = 1:2, reduction = "pca")


jpeg("dimplot.jpg", width = 650, height = 650)
DimPlot(GBM_4, reduction = "pca")
# 3. Close the file
dev.off()

jpeg("gimheat.jpg", width = 650, height = 650)
DimHeatmap(GBM_4, dims = 1, cells = 500, balanced = TRUE)
dev.off()

jpeg("Vlnplot.jpg", width = 650, height = 650)
DimHeatmap(GBM_4, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()


#Determine the ‘dimensionality’ of the dataset

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
GBM_4 <- JackStraw(GBM_4, num.replicate = 100)
GBM_4 <- ScoreJackStraw(GBM_4, dims = 1:15)

jpeg("Jackstraw.jpg", width = 650, height = 650)
JackStrawPlot(GBM_4, dims = 1:15)
dev.off()
#http://127.0.0.1:37799/graphics/plot_zoom_png?width=1920&height=1017
#Cluster the cells

GBM_4 <- FindNeighbors(GBM_4, dims = 1:10)
GBM_4 <- FindClusters(GBM_4, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(GBM_4), 5)


# Run non-linear dimensional reduction (UMAP/tSNE)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
GBM_4 <- RunUMAP(GBM_4, dims = 1:20)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
jpeg("dim_reduction.jpg", width = 650, height = 650)
DimPlot(GBM_4, reduction = "umap")
dev.off()

saveRDS(GBM_4, file = "GBM_2.rds")

#Finding differentially expressed features (cluster biomarkers)

# find all markers of cluster 1
cluster1.markers <- FindMarkers(GBM_4, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 300)

write.csv(cluster1.markers, 'GBM_2_markers.csv')


# find markers for every cluster compared to all remaining cells, report only the positive ones
GBM_4.markers <- FindAllMarkers(GBM_4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#GBM_4.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(GBM_4.markers, 'cluetr_markers2.csv')




cluster1.markers <- FindMarkers(GBM_4, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(GBM_4, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(GBM_4, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

jpeg("all_features.jpg", width = 650, height = 650)
FeaturePlot(GBM_4, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                                "CD8A"))
dev.off()

#top20 <- GBM_4.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
jpeg("topheat.jpg", width = 650, height = 650)
DoHeatmap(GBM_4, features = top20) + NoLegend()
dev.off()


saveRDS(GBM_4, file = "GBM_2_final.rds")
