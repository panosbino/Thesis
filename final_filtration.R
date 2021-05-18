library(Seurat)
library(tidyverse)
library(DoubletFinder)
library(scCATCH)


raw_data <- Read10X("~/Desktop/filtered_feature_bc_matrix")
raw_data <- CreateSeuratObject(counts = raw_data, project = "10x95", min.cells = 3, min.features = 500)
###################################################
# 
# individual_entropy <- function(p) {
#   ifelse(p == 0, 0, -log2(p)*p)
# }
# 
# entropy_data <-GetAssayData(all_data) %>% as.da() %>% mutate_all( funs((.)/sum(.))) %>% mutate_all(individual_entropy)
# 
# 
# 

###################################################



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
raw_data[["percent.mt"]] <- PercentageFeatureSet(raw_data, pattern = "^mt-")
all_data <- subset(raw_data, subset = nFeature_RNA < 10000 & percent.mt < 30)

sct_data <- SCTransform(all_data, vars.to.regress = c("nFeature_RNA","nCount_RNA")) 
all_data <- NormalizeData(all_data, normalization.method = "LogNormalize", scale.factor = 10000)


all_data <- FindVariableFeatures(all_data, selection.method = "vst", nfeatures = 3000)
all_data <- ScaleData(all_data, features = VariableFeatures(all_data))


all_data <- RunPCA(all_data, features = VariableFeatures(object = all_data))
ElbowPlot(all_data)

all_data <- FindNeighbors(all_data, dims = 1:20)
all_data <- FindClusters(all_data, resolution = 0.1)

all_data <- RunUMAP(all_data, dims = 1:20)
DimPlot(all_data, reduction = "umap", label = TRUE)
################################################################################################




scCatch.markers <- findmarkergenes(all_data, species = "Mouse", tissue = c("Brain","Blood", "Cerebellum", "Fetal brain",  "Neural tube", "Embryonic stem cell" ))
scCatch_annotations <- scCATCH(scCatch.markers, species = "Mouse", tissue = c("Brain","Blood", "Cerebellum", "Fetal brain",  "Neural tube", "Embryonic stem cell" ))



############################################################################################
# cluster.0 <- subset(all_data, idents = "0")
# cluster.0.cells <- cluster.0 %>% GetAssayData() %>% colnames()
#
# cluster.0.raw <- subset(raw_data, cells = cluster.0.cells)
#
# cluster.0.raw <- NormalizeData(cluster.0.raw, normalization.method = "LogNormalize", scale.factor = 10000)
#
#
# cluster.0.raw <- FindVariableFeatures(cluster.0.raw, selection.method = "vst", nfeatures = 3000)
# cluster.0.raw <- ScaleData(cluster.0.raw, features = VariableFeatures(all_data))
#
#
# cluster.0.raw <- RunPCA(cluster.0.raw, features = VariableFeatures(object = cluster.0.raw))
# ElbowPlot(all_data)
#
# cluster.0.raw <- FindNeighbors(cluster.0.raw, dims = 1:20)
# cluster.0.raw <- FindClusters(cluster.0.raw, resolution = 0.5)
#
# cluster.0.raw <- RunUMAP(cluster.0.raw, dims = 1:20)
# DimPlot(cluster.0.raw, reduction = "umap", label = TRUE)


#################################################################################################3
# VlnPlot(all_data, features = c("Sox2", "Fabp7"))
# VlnPlot(all_data, features = c("Htr2c", "Slc1a3"))
# #markers<- FindMarkers(pbmc, min.pct = 0.25 )
#
filtered_data <- subset(all_data, idents = c("0","4"), invert = TRUE)
#
sweep.res.data <- paramSweep_v3(filtered_data, PCs = 1:20, sct = FALSE)
sweep.data <- summarizeSweep(sweep.res.data, GT = FALSE)
bcmvn_data <- find.pK(sweep.data)
#
annotations <- filtered_data@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.16*nrow(filtered_data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#
filtered_data <- doubletFinder_v3(filtered_data, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#
# #data2 <- RunUMAP(data2, dims = 1:30 )
 DimPlot(filtered_data, reduction = "umap", label = TRUE, group.by = "DF.classifications_0.25_0.09_246")
#
#
#
no_doublets <- subset(filtered_data, DF.classifications_0.25_0.09_246 == "Singlet")
# no_doublets <- RunUMAP(no_doublets, dims = 1:50 )
# ElbowPlot(no_doublets)
# no_doublets <- FindNeighbors(no_doublets, dims = 1:12)
# no_doublets <- FindClusters(no_doublets, resolution = 0.5, algorithm = 1)
# no_doublets <- RunUMAP(no_doublets, dims = 1:12)
# DimPlot(no_doublets, reduction = "umap", label = TRUE)
#
#
# GetAssayData(no_doublets) %>% colnames() %>% write.table(file='~/Desktop/sub_bc.tsv', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
# write_rds(no_doublets, path =  "~/Desktop/no_doublets.rds")
#
# #
# # cluster_0 <- subset(no_doublets, idents = "0")
# # cluster_0 <- FindNeighbors(cluster_0, dims = 1:20)
# # cluster_0 <- FindClusters(cluster_0, resolution = 0.5)
# #
# # cluster_0 <- RunUMAP(cluster_0, dims = 1:20 )
# # DimPlot(cluster_0, reduction = "umap", label = TRUE)
#
#
#
# head(cluster2.markers, n = 5)