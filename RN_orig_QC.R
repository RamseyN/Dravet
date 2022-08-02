#Following instructions from following tutorial
#https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#######
#Here, I'm just making cuttoffs in the subset argument by eye. Is there a better way?
##########3
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 6500 & percent.mt < 8)
######

#Normalizing the data
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

seurat_obj <- NormalizeData(seurat_obj)
#Finding variables features for future PCA analysis
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#Scaling the data
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

#Running a PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
DimPlot(seurat_obj, reduction = "pca")
DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)
#seems like 12 dims is wehre it starts to fall apart... (28/7/22)




ElbowPlot(seurat_obj)
#Here 10 seems like a good number of PCs to use



seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters( seurat_obj, resolution = 0.5)
head(Idents(seurat_obj), 5)


seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = 'umap')
DimPlot (seurat_obj, reduction = 'umap', group.by = 'orig.ident')
#Save this seurat obj w/ res 0.5
saveRDS(seurat_obj, file = '/Users/ramsey.najm/GitHub/Projects/RN_SCN1A_PRT43/Objects/SeuratObj_0.5')

#DimPlot(seurat_obj, reduction = 'umap', group.by = 'orig.ident')
FeaturePlot(seurat_obj, features = c('DLX1','DLX2','DLX5','LHX6'), dims = c(2,1))


#Now we want to find all markers that define each cluster given this resolution of clustering
seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, 
                                     min.pct = 0.25, logfc.threshold = 0.25)

seurat_obj.markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)

seurat_obj.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()


#highlighting cells from specific cluster?
#n <- seurat_obj@meta.data
#n$cluster0 <- n$seurat_clusters
#n[rownames(n) %in% n$seurat_clusters, 'cluster0'] <- '0'
##n
#DimPlot(n,reduction = 'umap', group.by = '0')

seurat_obj@meta.data
clus0 <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters=='0',])
clus1 <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters=='1',])
clus2 <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters=='2',])
clus3 <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters=='3',])
clus4 <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters=='4',])
clus5 <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters=='5',])
clus6 <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters=='6',])
clus7 <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters=='7',])
clus8 <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters=='8',])
clus9 <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters=='9',])
clus10 <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters=='10',])
clus11 <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters=='11',])


DimPlot (seurat_obj, reduction = 'umap', cells.highlight = clus11)


#### Quantifying % Cells/Groups

perc.cluster <- prop.table(table(Idents(seurat_obj), seurat_obj$orig.ident), margin = 2)

print(perc.cluster)
write.table(perc.cluster)
write.csv(perc.cluster,"/Users/ramsey.najm/GitHub/Projects/RN_SCN1A_PRT43/perc_cluster.csv", row.names = TRUE)
######################