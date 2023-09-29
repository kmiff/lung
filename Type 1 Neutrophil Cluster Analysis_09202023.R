# Load Libraries
library(Seurat)
library(qs)
library(harmony)
library(cowplot)
library(ggplot2)
library(tidyverse)
library(RCurl)
library(dplyr)
library(clustree)  

#set wd
setwd("C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023")

#get seurat object 
lung <- readRDS("C:/Users/Katherine/Dropbox/KAM R Code/SCT_12122022/lung.rds")
lung

# set idents & isoalte cluster
Idents(lung)<- 'main.cluster'
levels(lung)
this_cluster <-'Type 1 Neutrophils'
lung_Neu1 <- subset(lung, idents = this_cluster)

# set default assay and normalize subcluster
DefaultAssay(lung_Neu1) <- "RNA"
Idents(lung_Neu1)<- 'orig.ident'
levels(lung_Neu1)
split_lung_Neu1 <- SplitObject(lung_Neu1 , split.by = "orig.ident")
split_lung_Neu1 <- lapply(X = split_lung_Neu1, 
                          FUN = SCTransform, 
                          method = "glmGamPoi", 
                          return.only.var.genes = FALSE,
                          vars.to.regress = c("percent.mito","S.Score","G2M.Score"))

var.features <- SelectIntegrationFeatures(object.list = split_lung_Neu1, nfeatures = 3000)                     
lung_Neu1 <- merge(x = split_lung_Neu1[[1]], y = split_lung_Neu1[2:length(split_lung_Neu1)], merge.data=TRUE)
VariableFeatures(lung_Neu1) <- var.features
lung_Neu1 <- RunPCA(lung_Neu1, verbose = FALSE)
lung_Neu1 <- RunUMAP(lung_Neu1, reduction = "pca", dims = 1:30)
lung_Neu1 <- FindNeighbors(lung_Neu1, reduction = "pca", dims = 1:30) 
lung_Neu1 <- FindClusters(lung_Neu1, resolution = c(0.1, 0.2, 0.3, 0.5, 0.6, 0.7,0.8, 1))
z1 <-DimPlot(lung_Neu1, group.by = "SCT_snn_res.0.1", label =TRUE)
z2 <-DimPlot(lung_Neu1, group.by = "SCT_snn_res.0.2", label =TRUE)
z3 <-DimPlot(lung_Neu1, group.by = "SCT_snn_res.0.3", label =TRUE)
z5 <-DimPlot(lung_Neu1, group.by = "SCT_snn_res.0.5", label =TRUE)
z6 <-DimPlot(lung_Neu1, group.by = "SCT_snn_res.0.6", label =TRUE)
z7 <-DimPlot(lung_Neu1, group.by = "SCT_snn_res.0.7", label =TRUE)
z8 <-DimPlot(lung_Neu1, group.by = "SCT_snn_res.0.8", label =TRUE)
z1.1 <-DimPlot(lung_Neu1, group.by = "SCT_snn_res.1", label =TRUE)
z1|z2|z3
z5|z6|z7
z8|z1.1
#decided to go with 0.5 resolution
lung_Neu1 <- FindClusters(lung_Neu1, resolution = 0.5)

lung_Neu1DP<-DimPlot(lung_Neu1, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
lung_Neu1DP

#check for batch effects
DimPlot(lung_Neu1, reduction = 'umap', group.by = 'orig.ident', label = FALSE)
#looks ok, so will not run harmony again

#save clusers
Idents(lung_Neu1) <- lung_Neu1$seurat_clusters
saveRDS(lung_Neu1, file = "C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/lung_Neu1.rds")

#check split umap and cell number per cluster & timpeoint
DimPlot(lung_Neu1, reduction = "umap", split.by = "orig.ident", ncol = 2, pt.size = 0.7)

# number of cells in each cluster
neu1_cluster_cell_number <- (table(Idents(lung_Neu1)))
neu1_cluster_cell_number
#results
0    1 
3924   87

# number of cells in each cluster by timepoint
neu1_cell_number_by_group<- table(Idents(lung_Neu1), lung_Neu1$orig.ident)
neu1_cell_number_by_group
# results
    N0  S12  S28  S3  T12  T28   T3
0  252   32 1598  174  101 1257  510
1    4    0   24    1    0    7   51

# now identify clusters markers
DefaultAssay(lung_Neu1) <- "RNA"
DefaultAssay(lung_Neu1)

#checked conserved markers, but none showed up for cluster 0
n1markers_cluster0 <- FindConservedMarkers(lung_Neu1,ident.1 = 0,grouping.var = 'orig.ident')
head(n1markers_cluster0)

n1markers_cluster1 <- FindConservedMarkers(lung_Neu1,ident.1 = 1,grouping.var = 'orig.ident')
head(n1markers_cluster1)

# Checking all makers 
DefaultAssay(lung_Neu1) <- "RNA"
DefaultAssay(lung_Neu1)

All_n1markers <- FindAllMarkers(object = lung_Neu1, only.pos = TRUE,logfc.threshold = 0.1, min.pct = )
library(dplyr)
All_n1markers  <- All_n1markers [ , c(6, 7, 2:4, 1, 5)]
All_n1markers<- All_n1markers %>%
  dplyr::arrange(cluster, p_val_adj)

View(All_n1markers)
write.table(All_n1markers, file = 'All_n1markers.csv', sep =',')
#checking find markers
Idents(lung_zero)<- 'seurat_clusters'
table(Idents(lung_zero))

cluster0.markers <- FindMarkers(lung_Neu1, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 10)
# results
p_val avg_log2FC pct.1 pct.2    p_val_adj
Lyst   2.157679e-69 -1.2883295 0.080 0.644 4.113183e-65
Retnlg 8.525897e-53 -1.1531717 0.081 0.586 1.625292e-48
Csf3r  2.377933e-51 -1.2807193 0.196 0.839 4.533053e-47
Mmp9   8.934636e-48 -1.1787128 0.108 0.621 1.703210e-43
Sorl1  1.042956e-46 -1.0536268 0.073 0.517 1.988186e-42
Xdh    4.777474e-43 -0.9421843 0.049 0.402 9.107299e-39
Hdc    7.447220e-43 -1.1205476 0.103 0.586 1.419664e-38
Picalm 3.263231e-42 -1.0882573 0.100 0.586 6.220697e-38
Malat1 9.201207e-42 -0.2618551 0.987 1.000 1.754026e-37
Tmcc1  1.254056e-41 -1.0588736 0.086 0.529 2.390607e-37

cluster1.markers <- FindMarkers(lung_Neu1, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 10)
#results
p_val avg_log2FC pct.1 pct.2    p_val_adj
Lyst   2.157679e-69  1.2883295 0.644 0.080 4.113183e-65
Retnlg 8.525897e-53  1.1531717 0.586 0.081 1.625292e-48
Csf3r  2.377933e-51  1.2807193 0.839 0.196 4.533053e-47
Mmp9   8.934636e-48  1.1787128 0.621 0.108 1.703210e-43
Sorl1  1.042956e-46  1.0536268 0.517 0.073 1.988186e-42
Xdh    4.777474e-43  0.9421843 0.402 0.049 9.107299e-39
Hdc    7.447220e-43  1.1205476 0.586 0.103 1.419664e-38
Picalm 3.263231e-42  1.0882573 0.586 0.100 6.220697e-38
Malat1 9.201207e-42  0.2618551 1.000 0.987 1.754026e-37
Tmcc1  1.254056e-41  1.0588736 0.529 0.086 2.390607e-37

#note: same markers showed up for both clusters using find markers,, not really sure there are two clusters
# thus decided not to go forward with analysis