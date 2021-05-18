library(Seurat)
library(tidyverse)
library(scCATCH)
library(SingleR)
library(Jmisc)
library(ggrepel)
library(loomR)
#library(SeuratDisk)

no_doublets <- read_rds("~/Desktop/no_doublets.rds")
all_data <- Read10X("~/Desktop/filtered_feature_bc_matrix")
all_data <- CreateSeuratObject(counts = all_data, project = "10x95", min.cells = 3, min.features = 500)

raw_data <- subset(all_data, cells = colnames(GetAssayData(no_doublets)))

raw_data <- SCTransform(raw_data, vars.to.regress = c("nFeature_RNA","nCount_RNA"))

raw_data <- FindVariableFeatures(raw_data, selection.method = "vst", nfeatures = 3000)

raw_data <- RunPCA(raw_data, features = VariableFeatures(object = raw_data))
ElbowPlot(data, ndims = 30)

raw_data <- FindNeighbors(raw_data, dims = 1:10)
raw_data <- FindClusters(raw_data, resolution = 0.5)

raw_data <- RunUMAP(raw_data, dims = 1:10)



plot1 <- VariableFeaturePlot(raw_data)
top10 <- head(VariableFeatures(raw_data), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

DimPlot(raw_data, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.5)


data <- NormalizeData(raw_data, normalization.method = "LogNormalize", scale.factor = 10000)


data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 3000)

all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)


data <- RunPCA(data, features = VariableFeatures(object = data))
ElbowPlot(data, ndims = 30)

data <- FindNeighbors(data, dims = 1:10)
data <- FindClusters(data, resolution = 0.5)

data <- RunUMAP(data, dims = 1:10)
DimPlot(data, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.5)

#data <- subset(data, idents = "0", invert = TRUE)

mkrs <- c("Sox2", "Rbfox3", "Atp1a3", "Map2", "Fgfr3", "Aldh1l1", "Aif1", "Rgs5", "Cldn5", "Erg","Foxj1" )

DotPlot(annotated_data, features = mkrs, cols = c("lightgrey", "red") ) + RotatedAxis()

VlnPlot(annotated_data, features = "pMR641") + NoLegend()


VlnPlot(annotated_data, features = c("Sox2", "Rbfox3", "Aif1", "Erg"), ncol = 2) + NoLegend()

RidgePlot(annotated_data, features = c("Rbfox3", "Map2"))

################### Auto-annotation ###############################
scCatch.markers <- findmarkergenes(data, species = "Mouse", tissue = c("Brain","Blood", "Cerebellum", "Fetal brain",  "Neural tube", "Embryonic stem cell" ))

scCatch_annotations <- scCATCH(scCatch.markers, species = "Mouse", tissue = c("Brain","Blood", "Cerebellum", "Fetal brain",  "Neural tube", "Embryonic stem cell" ))

scCatch_annotations <- scCATCH(scCatch.markers, species = "Mouse", tissue = c("Brain","Blood", "Cerebellum", "Fetal brain",  "Neural tube", "Embryonic stem cell" ))

##################################################################



markers <- FindAllMarkers(data)




FeaturePlot(data, features = radial_glia_markers)
FeaturePlot(data, features = oligodendrocytes)
FeaturePlot(data, features = immature_neurons)

FeaturePlot(data, features = "Lum")
FeaturePlot(data, features = "pMR641", cols = c("gray", "red"))
VlnPlot(data, features = "pMR641")

VlnPlot(data, features = c("2900040C04Rik", "Ttr"))
VlnPlot(data, features = c("Kcnj13", "Folr1"))


all.genes <- rownames(GetAssayData(data)) %>% as.data.frame()

filter(markers, cluster == "1")


GetAssayData(no_doublets) %>% colnames() %>% write.table(file='~/Desktop/sub_bc2.tsv', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
save.image("~/Desktop/after_annotation.RData")



##############################################################

counts <- GetAssayData(data)

singler <- CreateSinglerObject(counts=counts,
                               project.name="SingleR_annot", # choose
                               min.genes = 200, # ignore cells with fewer than 200 transcripts
                               technology = "10x", # choose
                               species = "Mouse",
                               citation = "", # choose
                               ref.list = list(mouse.rnaseq),
                               normalize.gene.length = FALSE,        # needed for full-length platforms (e.g. smartseq)
                               variable.genes = "de",  # see vignette
                               fine.tune = FALSE, # TRUE would take very long
                               reduce.file.size = TRUE, # leave out less-often used fields 
                               do.signatures = FALSE,
                               do.main.types = TRUE,
                               numCores = SingleR.numCores)

############################################################################################33

FeaturePlot(data, features = c("Sox2", "Sox9"))

FeaturePlot(data, features = "Aif1")
VlnPlot(data, features = c("Sox2", "Sox9", "Rbfox3","Dcx", "Aif1", "Erg"))



VlnPlot(data, features = c("Tubb3","Ttr", "Pax6", "Ki67"))
VlnPlot(data, features = c("Slc1a3","Eaat"))
VlnPlot(data,features = "Stc1")
VlnPlot(data, features = c("Tubb3","Rbfox3","Map2"), ncol = 2)


mature.neuronal.markers <- c("Tubb3","Rbfox3","Map2")
immature_neurons <- c("Dcx", "Neurod1", "Tbr1", "Pax6" )
oligodendrocytes <- c("Olig1", "Olig2", "Olig3", "Cldn11", "Sox10" )
radial_glia_markers <- c("Vim", "Nes", "Pax6", "Hes1", "Hes5","Cdh2", "Sox2" )
############################################################

m <-filter(markers, cluster ==10)
#Neuronal Progenitor Cells

cluster.ids <- c("NPC1", "NPC2","NPC3","NPC4","Radial glia", "Endothelial", "Microglia", "Astrocytes", "Pericyte", "Ependymal", "Neural Stem Cell")

names(cluster.ids) <- levels(data)
annotated_data <- RenameIdents(data, cluster.ids)

DimPlot(annotated_data, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE, label.size = 5) + NoLegend()

VlnPlot(annotated_data, features = "pMR641") + NoLegend()
#############################################################

write.csv(Cells(annotated_data), file = "cellID_obs.csv", row.names = FALSE)

write.csv(Embeddings(annotated_data, reduction = "umap"), file = "cell_embeddings.csv")


loom_data <- as.loom(annotated_data, filename = "10x95.loom")

annot <- annotated_data@active.ident %>% as.data.frame()

write.csv(annot, file = "annot_clusters.csv")

clusters <- annotated_data@meta.data %>% as.data.frame()
clusters <- rownames_to_column(clusters)
clusters2 <- clusters[c("rowname","seurat_clusters")]

write.csv(clusters2, file = "clusters.csv")

################################################################
a <- GetAssayData(data) %>% as.data.frame()
b <- subset(a, rownames(a) == "Ttr") %>% t() %>% as.data.frame()
c <- filter(b, Ttr > 2)
pct <- length(c$Ttr)/length(b$Ttr)

p <- ggplot(b, aes(x=Ttr)) + 
  geom_density()
p

genes <- rownames(GetAssayData(data)) %>% as.data.frame()

d <- GetAssayData(data) %>% as.data.frame() %>% rowMeans() %>% as.data.frame() 
d_max <- max(d$.)
gns <- rownames(d)
d <- addCol(d, Ttr = 0.5)
colnames(d) <- c("exp", "Ttr")
d <- mutate(d, logexp = log2(exp))
rownames(d) <- gns
e <- subset(d, rownames(d) == "Ttr") 

p <- ggplot(d, aes(x = Ttr, y=logexp)) + 
  geom_violin(aes( fill ="red", alpha = 0.8)) +
  geom_hline(yintercept = 1.920185, linetype = "dashed") +
  geom_point()
  
  #geom_jitter(shape = 16, position = position_jitter(0.07))
  
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.01)

p + scale_y_continuous(name="Expression") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", color = "black", size = 16),
        #axis.text.y = element_text(face = "bold", color = "black", size = 14),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none"
        )
#   geom_histogram(bins = 50, alpha = 0.8) +
#   geom_density(color = "blue") +
#   geom_vline(xintercept = 3.784717)
# p


# g <- GetAssayData(data) %>% as.data.frame() %>% colMeans() %>% as.data.frame() %>% rownames_to_column()
# colnames(g) <- c("cell", "exp")
# p <- ggplot(g, aes(x=cell,y=exp)) +
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
#   #geom_histogram(bins = 50, alpha = 0.8) +
#   geom_point(aes(x = 0.5, y=exp))
# p

######################################

gene_means <- data %>%  GetAssayData() %>% as.data.frame() %>% rowMeans() %>% as.data.frame()
genes <- rownames(gene_means)
colnames(gene_means) <- c("exp")
#gene_means <- sort(gene_means$exp, decreasing = TRUE) %>% as.data.frame()
rownames(gene_means) <- genes



Ttr <- subset(gene_means, rownames(gene_means) == "Ttr") %>% as.data.frame()
gene_means <- rownames_to_column(gene_means)
colnames(gene_means) <- c("gene","exp")
#Ttr <- subset(GetAssayData(data), rownames(GetAssayData(data) %>% as.data.frame()) == "Ttr") %>% t() %>% as.data.frame()



p <- ggplot() +
  geom_violin(data = gene_means, aes(x=1, y = log2(exp)), alpha = 0.7, fill = "Purple") +
  geom_jitter(data=gene_means[which(log2(gene_means$exp)>-15),], aes(x=1, y=log(exp)), alpha = 0.5, color = "black", width = 0.3, size = 0.2)+
  geom_jitter(data = Ttr, aes(x=1, y=log2(exp)), color = "red", width = 0.3, size = 2.5)
  #geom_jitter(data = Ttr, aes(x=1, y=log2(exp)), color = "Blue", width = 0.3, size = 2.5)
p + scale_y_continuous(name="Expression", limits = c(-15,5), breaks = c(-15,-10,-5,0,6)) +
  scale_x_continuous(breaks = c(), limits = c(0.3,1.7))+
  theme(axis.title.x = element_blank(),
        axis.ticks.y = element_line(size = 1),
        axis.title.y = element_text(face = "bold", color = "black", size = 16, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        #axis.text.y = element_text(face = "bold", color = "black", size = 14),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 0.8),
        legend.position = "none"
  ) 

gene_means_Sorted <- arrange(gene_means, desc(exp))

p <- ggplot(gene_means, aes(y = exp)) +
  geom_point() +
  geom_label()

p

markers2 <- c("Acot7", "Actl6b", "Adgrl1", "Ank2", "Ankrd12", "Aplp1")

############################ Plots

VlnPlot(data, features = c("Sox2", "Rbfox3", "Aif1", "Erg"), ncol = 2)

###############################################################################
#code to calc gene pct in cells

PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}


####

PrctCellExpringGene(annotated_data, genes = c("pMR641"))
raw_data <- FindVariableFeatures(raw_data)
plot1 <- VariableFeaturePlot(raw_data)
top10 <- head(VariableFeatures(raw_data), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
