# 0.0 References ----
#Used the following Tutorials to filter data & run Harmony with a Seurat Object
# http://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/docs/SeuratV3.html 
# https://www.youtube.com/watch?v=5HBzgsz8qyk
#https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html

# 1.0 OBJECT PREPARATION ----
## 1.1 Load Libraries----
library(Seurat)
library(qs)
library(harmony)
library(cowplot)
library(ggplot2)
library(tidyverse)
library(RCurl)
library(dplyr)
library(clustree)  
library(sctransform)
library(DoubletFinder)
library(RColorBrewer)
library(ggthemes)

## 1.2 Load the data---- 
#(Seurat Object from Yang before QC; already has QC metrics calculated)
lung <- qread("C:/Users/kathe/Dropbox/Mifflin Lung scRNA Seq Raw Data_April 2021/SeuratFiles_Data from Yang/combined_qc.qsave")
lung

## 1.3 SCI vs Uninjured metadata column----
#Renaming Idents to identify cell clusters
Idents(lung) <- "orig.ident"
levels (lung)
lung<- RenameIdents(lung, "N0" = "Uninjured")
lung<- RenameIdents(lung, "S12" = "Uninjured")
lung<- RenameIdents(lung, "S28" = "Uninjured")
lung<- RenameIdents(lung, "S3" = "Uninjured")
lung<- RenameIdents(lung, "T12" = "SCI")
lung<- RenameIdents(lung, "T28" = "SCI")
lung<- RenameIdents(lung, "T3" = "SCI")
# can now stash these identities as a new column in the meta data
lung <- StashIdent(lung, save.name = "injury.status")

#also added naive, sci and sham metadata colum
Idents(lung) <- "orig.ident"
levels (lung)
lung<- RenameIdents(lung, "N0" = "Naive")
lung<- RenameIdents(lung, "S12" = "Uninjured")
lung<- RenameIdents(lung, "S28" = "Uninjured")
lung<- RenameIdents(lung, "S3" = "Uninjured")
lung<- RenameIdents(lung, "T12" = "SCI")
lung<- RenameIdents(lung, "T28" = "SCI")
lung<- RenameIdents(lung, "T3" = "SCI")

# can now stash these identities as a new column in the meta data
lung <- StashIdent(lung, save.name = "naive.status")
View(lung@meta.data)

# 2.0 QUALITY CONTROL METRICS----
## 2.1Pre-Filtering Graphs ----
# set groups to visualize using Idents()
Idents(lung) <- "orig.ident"
levels (lung)
#Scatter graph
before_filter2 <-FeatureScatter(lung, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ggtitle("Before")
before_filter2
#rna graph
rna_before2 <- FeatureScatter(lung, feature1 = "nCount_RNA", feature2 = "percent.mito")+ggtitle("Before")
rna_before2
#Violin plot
Vln_before2 <- VlnPlot(lung, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol = 4)
Vln_before2 

## 2.2 Filter Data----
lung.filtered2 <- subset(lung, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito < 5)
lung.filtered2

## 2.3 Post-Filtering Graphs----
Idents(lung.filtered2) <- "orig.ident"
#Scatter graph
after_filter2 <-FeatureScatter(lung.filtered2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ggtitle("After")
before_filter2|after_filter2
#rna graph
rna_after2 <- FeatureScatter(lung.filtered2, feature1 = "nCount_RNA", feature2 = "percent.mito")+ggtitle("After")
rna_before2|rna_after2
# violin plot
after_Vln2 <- VlnPlot(lung.filtered2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol = 4)
Vln_before2|after_Vln2

##2.4 Assign Filtered object as Main object----
# ie make sure to use filtered lung in downstream analysis
lung <-lung.filtered2

## 2.5 Cell Cycle Scoring Check ----
# first Checking to see if need to account for this variance
#first normalize counts and check cell cycle scoring
#see this tutorial https://satijalab.org/seurat/articles/cell_cycle_vignette
# get genes(stored in Seurat library)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
lung <- NormalizeData(lung)
lung<- FindVariableFeatures(lung, selection.method = "vst")
lung <- ScaleData(lung)
lung <- RunPCA(lung)
lung_phaseDimP <- DimPlot(lung, reduction = 'pca', group.by = 'Phase')
lung_phaseDimP
## Do see a difference in with cells in each phase, so will need to regress out in SCTransform step

##2.6 Blood Contamination Filtering----
#looking at data PCA, saw that there was RBC contamination, so decided to score for blood markers
blood_genes_to_check = list(c("Hbb-bt","Hbb-bs", "Hba-a1", "Hba-a2", "Hba-2", "Hbb", "Hba-1", "Hbg2"))
lung <-AddModuleScore(lung, feature = blood_genes_to_check, name = "Blood.Enriched")
View(lung@meta.data)
#Visualize score distribution,
FeaturePlot(lung,features = "Blood.Enriched1", label = TRUE, repel = TRUE) +scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
#can now use this module score to regressout RBC contamination in SCTransform

## 2.7 Run DoubletFinder to Remove Data and Normalize with SCTransform----
dub_split <- SplitObject(lung, split.by = "orig.ident")
dub_split <- lung.filtered2[c("N0", "S12", "S28", "S3", "T12", "T28", "T3")]
### 2.7a Loop Through Samples to Find Doublets---- 
#recommended to run different samples individually as pK value may vary
#need to normalize and pre-cluster data first, then run doublet finder
#tutorials:
#general: https://www.youtube.com/watch?v=NqvAS4HgmrE 
#for loop: https://rpubs.com/kenneditodd/doublet_finder_example
for (i in 1:length(dub_split)){
  # print the sample we are on
  print(paste0("Lung ",i))
  
  # Pre-process seurat object with standard seurat workflow
  lung.sample <- SCTransform(dub_split[[i]], method = "glmGamPoi", 
                             return.only.var.genes = FALSE,
                             vars.to.regress = c("percent.mito","S.Score","G2M.Score","Blood.Enriched1"))
  lung.sample <- RunPCA(lung.sample, verbose = FALSE)
  lung.sample <- RunUMAP(lung.sample, dims = 1:30, verbose = FALSE)
  lung.sample <- FindNeighbors(object = lung.sample, dims = 1:30)              
  lung.sample <- FindClusters(object = lung.sample)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep_v3(lung.sample, PCs = 1:30, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- lung.sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(0.052*nrow(lung.sample@meta.data))  ## line above is assuming 5.2% doublet formation rate. this gives up the expected number of doublets
  # based off info here: https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
  # now adjusting the above number based off the prediced homotypic doublets
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  lung.sample <- doubletFinder_v3(seu = lung.sample, 
                                  PCs = 1:30, 
                                  pK = optimal.pk,
                                  pN = 0.25, 
                                  nExp = nExp.poi.adj,
                                  reuse.pANN = FALSE, 
                                  sct = TRUE)
  metadata <- lung.sample@meta.data
  colnames(metadata)[21] <- "doublet_finder" #[21] here is the metadata column with DF.classifications_0.25_0.04_XX that assigns doublet or not
  lung.sample@meta.data <- metadata 
  
  # subset and save
  lung.singlets <- subset(lung.sample, doublet_finder == "Singlet")
  dub_split[[i]] <- lung.singlets
  remove(lung.singlets)
  
}

### 2.7b Converge Singlets from Each Timepoint-----
lung.singlets <- merge(x = dub_split[[1]], y = dub_split[2:length(dub_split)])
### 2.7c Assign Variable Features----
#(see https://github.com/satijalab/seurat/issues/4145#issuecomment-786728442 )
obj.features <- SelectIntegrationFeatures(object.list = dub_split, nfeatures = 3000)
VariableFeatures(lung.singlets[["SCT"]]) <- obj.features

# 3.0 HAROMY INTEGRATION AFTER SCT NORMALIZATION----
# see tutorials/referencse (new in Seurat 5)
#https://satijalab.org/seurat/articles/integration_introduction.html
#https://satijalab.org/seurat/articles/seurat5_integration
# https://satijalab.org/seurat/reference/harmonyintegration
#already SCT noramlized each tiempoint object in doublet finder
# now just need to find variable features & integrate

## 3.1 Split Object by Timepoint----
#new method with Seurat 5
lung.singlets[["RNA"]] <- split(lung.singlets[["RNA"]], f = lung.singlets$orig.ident)
#run new PCA with new combined variable features
lung.singlets <- RunPCA(lung.singlets,npcs = 30, verbose = F)
#integrate
lung.singlets <- IntegrateLayers(object = lung.singlets, 
                                method = HarmonyIntegration,
                                normalization.method = "SCT",
                                orig.reduction = "pca", 
                                new.reduction = 'harmony',
                                verbose = FALSE)
lung.singlets <- FindNeighbors(lung.singlets, dims = 1:30, reduction = "harmony")
lung.singlets <- FindClusters(lung.singlets, resolution = 2, cluster.name = "harmony_clusters")
lung.singlets <- RunUMAP(lung.singlets, dims = 1:30, reduction = "harmony")
DimPlot(lung.singlets, group.by = c("orig.ident", "harmony_clusters"), combine = FALSE, label.size = 2)
# clusters look good and have no batch effects now
## 3.2 Merge layers----
#need to do before downstream analysis
lung.singlets[["RNA"]] <- JoinLayers(lung.singlets[["RNA"]])
Layers(lung.singlets[["RNA"]])

## 3.3 Save Object----
saveRDS(lung.singlets, file = "C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/lung.singets_11072023.rds")

# 4.0 CLUSTER IDENTIFICATION----
## 4.1 Find All Markers----
# tutorial (https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_integ_marker_identification.html)
lung <- lung.singlets
### 4.1a With RNA Assay----
#Set assay as RNA for DEG analyses
DefaultAssay(lung) <-'RNA'
#run FindALLMarkers to see cluster markers
All_markers <- FindAllMarkers(object = lung, only.pos = TRUE,logfc.threshold = 0.25, min.pct = )
# Reording output columns to be easier to read; putting cluster and gene first instead of last
# Combine markers with gene descriptions
All_markers  <- All_markers [ , c(6, 7, 2:4, 1, 5)]

# Order the rows by p-adjusted values
All_markers<- All_markers %>%
  dplyr::arrange(cluster, p_val_adj)

View(All_markers)
write.table(All_markers, file = 'All_markers_11072023.csv', sep =',')

### 4.1b Comparing to Makers from SCT Assay----
#see https://satijalab.org/seurat/articles/sctransform_v2_vignette 
lung <- PrepSCTFindMarkers(lung)
All_markers_SCT <- FindAllMarkers(object = lung, only.pos = TRUE,logfc.threshold = 0.25, min.pct = )
All_markers_SCTs  <- All_markers [ , c(6, 7, 2:4, 1, 5)]
All_marker_SCT<- All_markers_SCT %>%
  dplyr::arrange(cluster, p_val_adj)
View(All_markers_SCT)
write.table(All_markers_SCT, file = 'All_markers_SCT_11072023.csv', sep =',')
# SCT assay makers are more distinct, will use this method going forward
#have both RBC and degraded cell contamination in clusters 13, 22, 25

## 4.2 Remove Contamination with New Variable Features and Scaling---- 
#### 4.2a Isolate clusters----
Idents(lung)<- 'harmony_clusters'
levels(lung)
this_cluster <-c("0","1","10", "11", "12",  "14" ,"15", "16", "17",
                 "18", "19",
                 "2" , "20", "21", "23" ,"24", 
                 "26", "3" , "4" , "5" , "6" , "7" , "8","9")
obj <-subset(lung, idents = this_cluster)
DimPlot(obj)
### 4.2b Get New Variable Features----
obj.list <- SplitObject(obj, split.by = "orig.ident")
obj.var.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
VariableFeatures(obj[["SCT"]]) <- obj.var.features
#Remerge
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
Layers(obj[["RNA"]])

### 4.2c Re-scale and Re-cluster ----
# also regress out contamination/RBC data
obj <- ScaleData(obj,  assay = 'SCT', vars.to.regress = c("percent.mito","S.Score","G2M.Score","Blood.Enriched1"))
obj  <- RunPCA(obj  ,npcs = 30, verbose = F)
obj <- FindNeighbors(obj , dims = 1:30, reduction = "harmony")
obj <- FindClusters(obj , resolution = 2, cluster.name = "harmony_clusters_2")
obj <- RunUMAP(obj , dims = 1:30, reduction = "harmony")
DimPlot(obj , group.by = ("harmony_clusters_2"), label = TRUE, label.size = 4)

### 4.2d New Cluster Markers----
Idents(obj )<- "harmony_clusters_2"
levels(obj )
obj  <- PrepSCTFindMarkers(obj)
#since re-ran PrepSCT, make sure to put "recorrect_umi = FALSE"
#see https://satijalab.org/seurat/articles/sctransform_v2_vignette.html#identify-differential-expressed-genes-across-conditions-1
All_markers_obj <- FindAllMarkers(object = obj, recorrect_umi = FALSE, only.pos = TRUE,logfc.threshold = 0.25, min.pct = )
All_markers_obj  <- All_markers_obj[ , c(6, 7, 2:4, 1, 5)]
All_marker_obj <- All_markers_obj  %>%
  dplyr::arrange(cluster, p_val_adj)
View(All_markers_obj )
write.table(All_markers_obj , file = 'All_markers_obj_11112023.csv', sep =',')
# still have degraded cells in clusters 16 and 12, so will remove

## 4.3 Second Contamination Removal---- 
#### 4.3a Isolate clusters----
Idents(obj )<- "harmony_clusters_2"
levels(obj )
this_cluster <-c("0","1","2","3","4","5","6","7","8","9",
                 "10","11","13","14","15","17","18","19",
                 "20","21","22","23")
obj2 <-subset(obj, idents = this_cluster)
DimPlot(obj2)
### 4.3b Get New Variable Features----
obj.list2 <- SplitObject(obj2, split.by = "orig.ident")
obj.var.features2 <- SelectIntegrationFeatures(object.list = obj.list2, nfeatures = 3000)
VariableFeatures(obj2[["SCT"]]) <- obj.var.features2
#remerge
obj2[["RNA"]] <- JoinLayers(obj2[["RNA"]])
Layers(obj2[["RNA"]])

### 4.3c Re-scale and Re-cluster ----
# also regress out contamination/RBC data
obj2 <- ScaleData(obj2,  assay = 'SCT', vars.to.regress = c("percent.mito","S.Score","G2M.Score","Blood.Enriched1"))
obj2 <- RunPCA(obj2  ,npcs = 30, verbose = F)
obj2 <- FindNeighbors(obj2 , dims = 1:30, reduction = "harmony")
obj2 <- FindClusters(obj2 , resolution = 2, cluster.name = "harmony_clusters_obj2")
obj2 <- RunUMAP(obj2 , dims = 1:30, reduction = "harmony")
DimPlot(obj2 , group.by = ("harmony_clusters_obj2"), label = TRUE, label.size = 4)

### 4.3d New Cluster Markers
Idents(obj2)<- "harmony_clusters_obj2"
levels(obj2)
obj2  <- PrepSCTFindMarkers(obj2)
All_markers_obj2 <- FindAllMarkers(object = obj2, only.pos = TRUE,logfc.threshold = 0.25, min.pct = )
All_markers_obj2  <- All_markers_obj2[ , c(6, 7, 2:4, 1, 5)]
All_marker_obj2 <- All_markers_obj2  %>%
  dplyr::arrange(cluster, p_val_adj)
View(All_markers_obj2 )
write.table(All_markers_obj2 , file = 'All_markers_obj2_11122023.csv', sep =',')

## 4.4 Assign Cell Types to Cluster Idents----
obj2 <- RenameIdents(obj2, "0" = "T Lymphocytes" )
obj2 <- RenameIdents(obj2, "1" = "T Lymphocytes" )
obj2 <- RenameIdents(obj2, "2" = "Monocytes")
obj2 <- RenameIdents(obj2, "3" = "B Lymphocytes")
obj2 <- RenameIdents(obj2, "4" = "General Capillary")
obj2 <- RenameIdents(obj2, "5" = "Alveolar Macrophages")
obj2 <- RenameIdents(obj2, "6" = "Capillary Aerocytes")
obj2 <- RenameIdents(obj2, "7" = "NK Cells")
obj2 <- RenameIdents(obj2, "8" = "Neutrophils")
obj2 <- RenameIdents(obj2, "9" = "Club Cells")
obj2 <- RenameIdents(obj2, "10" = "Mesothelial")
obj2 <- RenameIdents(obj2, "11" = "Neutrophils")
obj2 <- RenameIdents(obj2, "12" = "Fibroblasts")
obj2 <- RenameIdents(obj2, "13" = "Capillary Aerocytes")
obj2 <- RenameIdents(obj2, "14" = "Interstitial Macrophages")
obj2 <- RenameIdents(obj2, "15" = "Neutrophils")
obj2 <- RenameIdents(obj2, "16" = "Monocytes")
obj2 <- RenameIdents(obj2, "17"= "General Capillary")
obj2 <- RenameIdents(obj2, "18" = "Col14a1+ Fibroblasts" )
obj2 <- RenameIdents(obj2, "19" = "B Lymphocytes")
obj2 <- RenameIdents(obj2, "20" = "Alveolar Macrophages")

### 4.4a New Labelled Polychrome DimPlot----
labelled_DP <- DimPlot(obj2, reduction = "umap", label = TRUE, label.size = 5, cols = "polychrome")
labelled_DP
#unlablled
unlabelled_DP <- DimPlot(obj2, reduction = "umap", label = FALSE, label.size = 5, cols = "polychrome")
unlabelled_DP
### 4.4b Final Cell Markers----
obj2  <- PrepSCTFindMarkers(obj2)
All_markers_Fobj2 <- FindAllMarkers(object = obj2, only.pos = TRUE,logfc.threshold = 0.25, min.pct = )
All_markers_Fobj2  <- All_markers_Fobj2[ , c(6, 7, 2:4, 1, 5)]
All_marker_obj2 <- All_markers_Fobj2  %>%
  dplyr::arrange(cluster, p_val_adj)
View(All_markers_Fobj2 )
write.table(All_markers_Fobj2 , file = 'All_markers_Fobj2_11132023.csv', sep =',')
### 4.4c Stash Cell Type Idents----
obj2[["main.cluster"]] <- Idents(obj2)

# 5.0 CELL TYPE MARKER HEATMAP----
## 5.1 Reorder Idents & Assign Cell Markers----
Idents(obj2) <- "main.cluster"
levels(obj2)
obj2$main.cluster <- factor(obj2$main.cluster, levels = c("Col14a1+ Fibroblasts","Fibroblasts","Mesothelial", 
                                                          "Capillary Aerocytes","General Capillary",
                                                          "Club Cells",
                                                          "Neutrophils",
                                                          "NK Cells", "T Lymphocytes", "B Lymphocytes",
                                                          "Alveolar Macrophages","Interstitial Macrophages","Monocytes"))
Idents(obj2) <- "main.cluster"
levels(obj2)
final.cluster.markers <- c('Col14a1', 'Mgp', 'Col3a1','Col1a2','Ogn',
                           'Inmt','Fmo2','Limch1', 'Macf1', 'Mfap4',
                           'Dcn',  'C3','Msln','Igfbp5', 'Igfbp4',
                           'Prx','Car4','Ednrb','Igfbp7','Ly6e',
                           'Plvap','Lyve1', 'Cd93','Gpihbp1','Kit',
                           'Scgb1a1','Cyp2f2','Scgb3a1','Scgb3a2','Eif3m',
                           'Retnlg','S100a9','S100a8', 'Ly6g','Ifitm1',
                           'Ccl5', 'Gzma', 'Nkg7', 'Klrb1c', 'Klra4',
                           'Trbc1','Rag1','Bcl11b','Cd8a','Rhoh',
                           'Igkc','Cd74','Ms4a1','Ighm','Bank1',
                           'Ear2', 'Ear1', 'Ftl1', 'Ctsd','Chil3',
                           'Apoe', 'C1qb','C1qc', 'C1qa', 'H2-Aa',
                           'Pou2f2', 'Csf1r', 'Cx3cr1','Clec4a3', 'Treml4')
## 5.2 Average Expression Heatmap----
all_lung <- AverageExpression(obj2, features = (final.cluster.markers), return.seurat = T, assay = 'SCT')
lung_HeatM <-DoHeatmap(all_lung, features = rev(final.cluster.markers), draw.lines = FALSE) + scale_fill_gradientn(colors = c("black", "midnightblue", "magenta4","orange", "yellow"))
lung_HeatM

## 5.3 Marker Feature Plots----
Col14FP <- FeaturePlot(obj2, features = c("Mgp")) & scale_colour_gradientn(colours = c("lightgrey", "violetred"))
Col14FP 
FibroFP <- FeaturePlot(obj2, features = c("Inmt")) & scale_colour_gradientn(colours = c("lightgrey", "violetred"))
FibroFP
MesoFP <-FeaturePlot(obj2, features = c("Dcn")) & scale_colour_gradientn(colours = c("lightgrey", "violetred"))
MesoFP
aCapFP <- FeaturePlot(obj2, features = c("Ly6a")) & scale_colour_gradientn(colours = c("lightgrey", "violetred"))
aCapFP
gcapFP <- FeaturePlot(obj2, features = c("Cd93")) & scale_colour_gradientn(colours = c("lightgrey", "violetred"))
gcapFP
clubFP <-FeaturePlot(obj2, features = c("Scgb1a1")) & scale_colour_gradientn(colours = c("lightgrey", "violetred"))
clubFP
neuFP<- FeaturePlot(obj2, features = c("S100a9")) & scale_colour_gradientn(colours = c("lightgrey", "violetred"))
neuFP
NkFP <-FeaturePlot(obj2, features = c("Ccl5")) & scale_colour_gradientn(colours = c("lightgrey", "violetred"))
NkFP
tFP <- FeaturePlot(obj2, features = c("Rag1")) & scale_colour_gradientn(colours = c("lightgrey", "violetred"))
tFP
bFP<-FeaturePlot(obj2, features = c("Igkc")) & scale_colour_gradientn(colours = c("lightgrey", "violetred"))
bFP
amFP <- FeaturePlot(obj2, features = c("Ftl1")) & scale_colour_gradientn(colours = c("lightgrey", "violetred"))
amFP
imFP<- FeaturePlot(obj2, features = c("Cd74")) & scale_colour_gradientn(colours = c("lightgrey", "violetred"))
imFP
monoFP<- FeaturePlot(obj2, features = c("Pou2f2")) & scale_colour_gradientn(colours = c("lightgrey", "violetred"))
monoFP

#save object
saveRDS(obj2, file = "C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/obj2.rds")

## 5.4 Table of Total Cell Counts----
cellcounts<- table(Idents(obj2))
cellcounts
write.table(cellcounts , file = 'cellcounts_11132023.csv', sep =',')
#counts by orig.ident
cluster_count<- table(Idents(obj2), obj2$orig.ident)    
cluster_count
write.table(cluster_count , file = 'cluster_count_11132023.csv', sep =',')

# 6.0 FREQUENCY PLOTS----
# plot clusters as proportion or percentage of total cell types at each timepoint
## 6.1 Reorder Timepoints----
obj2$orig.ident <- factor(obj2$orig.ident, levels = c("N0", "S12", "T12", "S3", "T3", "S28", "T28"))
Idents(obj2) <- "orig.ident"
## 6.2 Set Colour Palette----
colors_many <- palette.colors(n = NULL, palette = "Polychrome 36", alpha=1, recycle = FALSE)
## 6.3 Plot Frequencies by Timepoint----
all_groups_FreP <- ggplot(obj2@meta.data, aes(x=orig.ident, fill=main.cluster)) + geom_bar(position = "fill")+ scale_fill_manual(values = colors_many)
all_groups_FreP




