#Used the following Tutorials to filter data & run Harmony with a Seurat Object
# http://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/docs/SeuratV3.html 
# https://www.youtube.com/watch?v=5HBzgsz8qyk
#https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html

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

##QUALITY CONTROL METRICS
# Load the data (Seurat Object from Yang before QC; already has QC metrics calculated)
memory.limit(1e10)
lung <- qread("C:/Users/Kevin/Dropbox/Mifflin Lung scRNA Seq Raw Data_April 2021/SeuratFiles_Data from Yang/combined_qc.qsave")
lung
# set groups to visualize using Idents()
Idents(lung) <- "orig.ident"
View(lung@meta.data)
Vln_before2 <- VlnPlot(lung, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol = 4)
Vln_before2 
before_filter2 <-FeatureScatter(lung, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ggtitle("Before")
lung.filtered2 <- subset(lung, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito < 5)
lung.filtered2
# set groups to visualize using Idents()
Idents(lung.filtered2) <- "orig.ident"
after_filter2 <-FeatureScatter(lung.filtered2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ggtitle("After")
before_filter2|after_filter2
rna_before2 <- FeatureScatter(lung, feature1 = "nCount_RNA", feature2 = "percent.mito")+ggtitle("Before")
rna_before2
rna_after2 <- FeatureScatter(lung.filtered2, feature1 = "nCount_RNA", feature2 = "percent.mito")+ggtitle("After")
rna_before2|rna_after2
after_Vln2 <- VlnPlot(lung.filtered2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol = 4)
Vln_before2|after_Vln2


##STANDARD NORMALIZATION METHOD with HARMONY Bath Correction (for comparison)
# Normalizing Data
lung.filtered2 <- NormalizeData(lung.filtered2)
# Finding Variable Features
lung.filtered2 <- FindVariableFeatures(lung.filtered2, selection.method = "vst", nfeatures = 2000)
#identifying & plotting top 10 features
top10_2 <-head(VariableFeatures(lung.filtered2),10)
top10_2
vf_plot_2 <-VariableFeaturePlot(lung.filtered2) 
vf_plot_2
#Scaling the Data (removes unwanted sources of variation)
all.genes <-rownames(lung.filtered2)
lung.filtered2 <- ScaleData(lung.filtered2, features = all.genes)
#note: How to check object slots
str(lung.filtered2)
# Performing linear dimension reduction
lung.filtered2 <- RunPCA(lung.filtered2, features = VariableFeatures(lung.filtered2))
# determining dimensionally of data (way to determine what PC capture most of the variability)
ElbowPlot(lung.filtered2)
#note better to err on side of caution and pick more PC than less
# here will go with 20 despite 16 looking like ok, to be cautious
lung.filtered2 <- RunUMAP(lung.filtered2, dims = 1:20, reduction = 'pca')
# creating a dim plot to see groups and any potential batch effects
pre_harmonyDP2 <- DimPlot(lung.filtered2, reduction = 'umap', group.by = 'orig.ident')+ ggtitle("Before")
pre_harmonyDP2
# do appear to have bacth effects. Note had to change default parameters for larger data set
lung.harmony2 <- lung.filtered2 %>%
  RunHarmony(group.by.vars = 'orig.ident',kmeans_init_nstart=20, kmeans_init_iter_max=500)
# now can see Harmony reduction in slot recall fxn
lung.harmony2@reductions
# to see embeddings
lung.harmony2.embed <- Embeddings(lung.harmony2, "harmony")
lung.harmony2.embed[1:10,1:10]

#create new umap with harmony embeds
lung.harmony2 <- lung.harmony2 %>%
  RunUMAP(reduction = 'harmony',dim = 1:20)%>%
  FindNeighbors(reduction = 'harmony', dims = 1:20)
#check multiple resolutions to find best one. avoid over clustering
lung.harmony2 <- FindClusters(lung.harmony2, resolution = c(0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 1))
res0.1_2 <-DimPlot(lung.harmony2, group.by = "RNA_snn_res.0.1", label =TRUE)
res0.2_2 <-DimPlot(lung.harmony2, group.by = "RNA_snn_res.0.2", label =TRUE)
res0.3_2 <-DimPlot(lung.harmony2, group.by = "RNA_snn_res.0.3", label =TRUE)
res0.5_2 <-DimPlot(lung.harmony2, group.by = "RNA_snn_res.0.5", label =TRUE)
res0.6_2 <-DimPlot(lung.harmony2, group.by = "RNA_snn_res.0.6", label =TRUE)
res0.7_2 <-DimPlot(lung.harmony2, group.by = "RNA_snn_res.0.7", label =TRUE)
res1_2 <-DimPlot(lung.harmony2, group.by = "RNA_snn_res.1", label =TRUE)
res0.1_2|res0.2_2|res0.3_2
res0.5_2|res0.6_2|res0.7_2|res1_2
#decided to go with 0.8
lung.harmony2 <- FindClusters(lung.harmony2, resolution = c(0.8))
harmony_DP2 <- DimPlot(lung.harmony2, reduction = 'umap', group.by = 'orig.ident')+ ggtitle("After")
harmony_DP2
pre_harmonyDP2|harmony_DP2

#SCTRANSFORM NORMALIZATION & INTEGRATION
#Cell Cycle Scoring & Regression- Checking to see if need to account for this variance
#first normalize counts and check cell cycle scoring
# get genes(stored in Seurat library)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# split the dataset into a list of  seurat objects by timepoint and condition
lung.list <- SplitObject(lung, split.by = "orig.ident")
lung <- NormalizeData(lung)
lung<- FindVariableFeatures(lung, selection.method = "vst")
lung <- ScaleData(lung)
lung <- RunPCA(lung)
lung_phaseDimP <- DimPlot(lung, reduction = 'pca', group.by = 'Phase', split.by = 'Phase')
lung_phaseDimP
## Do see a difference in with cells in each phase, so need to regress out
#looking at data PCA, saw that there was RBC contamination, so decided to score for blood markers
# Adding Score for blood cells because saw RBC contamination previously
blood_genes_to_check = list(c("Hbb-bt","Hbb-bs", "Hba-a1", "Hba-a2", "Hba-2", "Hbb", "Hba-1", "Hbg2"))
lung <-AddModuleScore(lung, feature = blood_genes_to_check, name = "Blood.Enriched")
View(lung@meta.data)
#Using SCT Transform to Normalize and Integrate
## Now will run a loop to normalize data, cell cycle score,blood score and SCTransform on split data
## need to split because need to run on each sample individually
split_lung <- SplitObject(lung, split.by = "orig.ident")
split_lung <- split_lung[c("N0", "S12", "S28", "S3", "T12", "T28", "T3")]
for (i in 1:length(split_lung)) {
  split_lung[[i]] <- SCTransform(split_lung[[i]], vars.to.regress = c("percent.mito","S.Score","G2M.Score", "Blood.Enriched1"))
}
#save object
saveRDS(split_lung, file = "C:/Users/Katherine/Dropbox/Mifflin Lung scRNA Seq Raw Data_April 2021/Harmony_08172022/split_lung.rds")
# Select the most variable features to use for integration 
#(Canonical correlation analysis (CCA))
integ_features <- SelectIntegrationFeatures(object.list = split_lung,
                                            nfeatures = 3000) 
# Prepare the SCT list object for integration i.e identify anchors
split_lung <- PrepSCTIntegration(object.list = split_lung, 
                                 anchor.features = integ_features)
# Identify anchors & Find mean nearest neighours (MNN, aka best buddies)& filter
integ_anchors <- FindIntegrationAnchors(object.list = split_lung, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
#save integrated anchors
saveRDS(integ_anchors, file = "C:/Users/Kevin/Desktop/KAM R Code/integ_anchors.rds")

# Integrate across conditions
lung_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
saveRDS(lung_integrated, file = "C:/Users/Kevin/Desktop/KAM R Code/lung_integrated.rds")

#visualize newly integrated  PCA & UMAP and check for batch effects
lung_integrated <- RunPCA(object = lung_integrated)
lung_iuntegratedPCA <- PCAPlot(lung_integrated, split.by = "orig.ident")
lung_iuntegratedPCA # will help visulalize if have good overlay of all conditions in each PCA
lung_integrated <- RunUMAP(lung_integrated, 
                           dims = 1:40,
                           reduction = "pca")
lung_integrated_UMAP <- DimPlot(lung_integrated)
#Note: Want to see good integration of conditions
lung_integrated_UMAP
#Check Split Plot
lung_integrated_UMAPSplit <- DimPlot(lung_integrated, split.by = "orig.ident")
lung_integrated_UMAPSplit
saveRDS(lung_integrated, file = "C:/Users/Kevin/Dropbox/KAM R Code/lung_integrated.rds")


##COMBINING SCT & HARMONY
# I think there may  be still a S28/T28 batch effect so, ran  SCTransform again and combined with Harmony too remove
#tutorial here:https://github.com/satijalab/seurat/issues/4896
split_lung <- SplitObject(lung, split.by = "orig.ident")
split_lung <- lapply(X = split_lung, 
                     FUN = SCTransform, 
                     method = "glmGamPoi", 
                     return.only.var.genes = FALSE,
                     vars.to.regress = c("percent.mito","S.Score","G2M.Score","Blood.Enriched1"))

saveRDS(split_lung, file = "C:/Users/Kevin/Dropbox/KAM R Code/SCT_Harmony_11062022/split_lung.rds")
var.features <- SelectIntegrationFeatures(object.list = split_lung, nfeatures = 3000)                     
lung_sct <- merge(x = split_lung[[1]], y = split_lung[2:length(split_lung)], merge.data=TRUE)
VariableFeatures(lung_sct) <- var.features
lung_sct <- RunPCA(lung_sct, verbose = FALSE)
lung_sct <- RunHarmony(lung_sct, assay.use="SCT", group.by.vars = "orig.ident", plot_convergence = TRUE)
lung_sct@reductions
#convergence plot looks ok despite quick trans error, so moveforward and check if have embeddings
#didn't alter standard parameters as didn't appear to make a difference, also get convergence after 7 iterations no matter what
harmony_embeddings <- Embeddings(lung_sct, 'harmony')
harmony_embeddings[1:5, 1:5]
#visualizing if batch effects are gone, appear to be
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = lung_sct, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = lung_sct, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)
#making umap
lung_sct <- RunUMAP(lung_sct, reduction = "harmony", dims = 1:30)
lung_sct <- FindNeighbors(lung_sct, reduction = "harmony", dims = 1:30) %>% FindClusters()
lung_sctDP<-DimPlot(lung_sct, group.by = "orig.ident")
lung_sctDP
#comparing to SCT integration method without harmony
options(repr.plot.height = 5, repr.plot.width = 12)
p1I <- DimPlot(object = lung_integrated, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2I <- VlnPlot(object = lung_integrated, features = "nFeature_SCT", group.by = "orig.ident", pt.size = .1)
plot_grid(p1I,p2I)
saveRDS(lung_sct, file = "C:/Users/Kevin/Dropbox/KAM R Code/SCT_Harmony_11062022/lung_sct.rds")

#Comparing all normalization and integration UMAP plots to see best method
#looks like combined method give best separation of clusters and best integration and no batch effects
pre_harmonyDP2 + ggtitle("Before Correction")|harmony_DP2 +ggtitle("Harmony")|lung_integrated_UMAP +ggtitle("SCT Integration")|lung_sctDP+ggtitle("SCT w Harmony")


#Identifying DOUBLETS WITH DOUBLET FINDER
library(DoubletFinder)
library(tidyverse)
library(dplyr)
## Doublet Finder (https://www.youtube.com/watch?v=NqvAS4HgmrE)
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_lung <- paramSweep_v3(lung_sct, PCs = 1:30, sct = TRUE)
sweep.stats_lung <- summarizeSweep(sweep.res.list_lung, GT = FALSE)
bcmvn_lung <- find.pK(sweep.stats_lung)

View(sweep.res.list_lung)
View(sweep.stats_lung)
View(bcmvn_lung)
     
# plot to visualize results
ggplot(bcmvn_lung, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

# from this graph, can see pK value 0.04 corresponds to the highest BC metric
#now programatically storing the optical pk value using tidyverse to a variable called pK
pK <- bcmvn_lung %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 

# can check selects correct one with: 
View(pK)

# choose first value from returned list, convert to numeric value, and store
pK <- as.numeric(as.character(pK[[1]]))

# now have values for pK and pN needed to run doublet finder
#just need to determine estimate value from user guide of reagent kit
annotations <- lung_sct@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) # modeling based on annotated seurat clusters in line above                                                       
nExp_poi <- round(0.052*nrow(lung_sct@meta.data))  
## line above is assuming 5.2% doublet formation rate. this gives up the expected number of doublets
# based off info here: https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# now adjusting the above number based off the prediced homotypic doublets
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


dub_split <- SplitObject(lung_sct, split.by = "orig.ident")
dub_split <- dub_split[c("N0", "S12", "S28", "S3", "T12", "T28", "T3")]
for (i in 1:length(dub_split)) {
  annotations <- dub_split[[i]]@meta.data$seurat_clusters 
  homotypic.prop <- modelHomotypic(annotations) # modeling based on annotated seurat clusters in line above                                                       
  nExp_poi <- round(0.052*nrow(dub_split[[i]]@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  dub_split[[i]] <- doubletFinder_v3(dub_split[[i]], 
                                                     PCs = 1:30, 
                                                     pN = 0.25, 
                                                     pK = pK, 
                                                     nExp = nExp_poi.adj,
                                                     reuse.pANN = FALSE, sct = TRUE)
  
}
lung_sct_dubfinder <- merge(x = dub_split[[1]] , y = dub_split[2:length(dub_split)], merge.data=TRUE)
View(lung_sct_dubfinder@meta.data)
saveRDS(lung_sct_dubfinder, file = "C:/Users/Kevin/Dropbox/KAM R Code/SCT_Harmony_11062022/lung_sct_dubfinder.rds")

#subset merged columns to reassign idents so can create one doublet id column in meta data
Idents(lung_sct_dubfinder) <- "orig.ident"
levels(lung_sct_dubfinder)
this_subcluster <-"N0"
N0.idents <-subset(lung_sct_dubfinder, idents = this_subcluster)
Idents(N0.idents) <- "DF.classifications_0.25_0.03_343"
levels(N0.idents)
# can now stash these identities as a new column in the meta data
N0.idents <- StashIdent(N0.idents, save.name = "dub.idents")
view(N0.idents@meta.data)

Idents(lung_sct_dubfinder) <- "orig.ident"
levels(lung_sct_dubfinder)
this_subcluster <-"S12"
S12.idents <-subset(lung_sct_dubfinder, idents = this_subcluster)
levels(S12.idents)
view(S12.idents@meta.data)
Idents(S12.idents) <- "DF.classifications_0.25_0.03_104"
levels(S12.idents)
S12.idents <- StashIdent(S12.idents, save.name = "dub.idents")
view(S12.idents@meta.data)

Idents(lung_sct_dubfinder) <- "orig.ident"
levels(lung_sct_dubfinder)
this_subcluster <-"S28"
S28.idents <-subset(lung_sct_dubfinder, idents = this_subcluster)
levels(S28.idents)
view(S28.idents@meta.data)
Idents(S28.idents) <- "DF.classifications_0.25_0.03_2106"
levels(S28.idents)
S28.idents <- StashIdent(S28.idents, save.name = "dub.idents")
view(S28.idents@meta.data)

Idents(lung_sct_dubfinder) <- "orig.ident"
levels(lung_sct_dubfinder)
this_subcluster <-"S3"
S3.idents <-subset(lung_sct_dubfinder, idents = this_subcluster)
levels(S3.idents)
view(S3.idents@meta.data)
Idents(S3.idents) <- "DF.classifications_0.25_0.03_286"
levels(S3.idents)
S3.idents <- StashIdent(S3.idents, save.name = "dub.idents")
view(S3.idents@meta.data)

Idents(lung_sct_dubfinder) <- "orig.ident"
levels(lung_sct_dubfinder)
this_subcluster <-"T12"
T12.idents <-subset(lung_sct_dubfinder, idents = this_subcluster)
levels(T12.idents)
view(T12.idents@meta.data)
Idents(T12.idents) <- "DF.classifications_0.25_0.03_100"
levels(T12.idents)
T12.idents <- StashIdent(T12.idents, save.name = "dub.idents")
view(T12.idents@meta.data)

Idents(lung_sct_dubfinder) <- "orig.ident"
levels(lung_sct_dubfinder)
this_subcluster <-"T28"
T28.idents <-subset(lung_sct_dubfinder, idents = this_subcluster)
levels(T28.idents)
view(T28.idents@meta.data)
Idents(T28.idents) <- "DF.classifications_0.25_0.03_1675"
levels(T28.idents)
T28.idents <- StashIdent(T28.idents, save.name = "dub.idents")
view(T28.idents@meta.data)

Idents(lung_sct_dubfinder) <- "orig.ident"
levels(lung_sct_dubfinder)
this_subcluster <-"T3"
T3.idents <-subset(lung_sct_dubfinder, idents = this_subcluster)
levels(T3.idents)
view(T3.idents@meta.data)
Idents(T3.idents) <- "DF.classifications_0.25_0.03_926"
levels(T3.idents)
T3.idents <- StashIdent(T3.idents, save.name = "dub.idents")
view(T3.idents@meta.data)

#merge together
merged_lung_sct_dubfinder <- merge(x = N0.idents , y = c(S12.idents, S28.idents, S3.idents, T12.idents, T28.idents, T3.idents), merge.data=TRUE)

# now need to recluster data on just the singlets, so use subset function to remove doublets
lung_sct_singlets <- subset(merged_lung_sct_dubfinder, subset = dub.idents == "Singlet")
# NOw re-running SCT Transform to Normalize and Integrate since will be new data 
## Now will run a loop to normalize data, cell cycle score,blood score and SCTransform on split data
#also note this tutorial for variable features and dims: https://satijalab.org/seurat/articles/sctransform_vignette.html
# I think there may  be still a S28/T28 batch effect so, ran  SCTransform again and combined with Harmony too remove
#tutorial here:https://github.com/satijalab/seurat/issues/4896
split_lung_singlets <- SplitObject(lung_sct_singlets, split.by = "orig.ident")
split_lung_singlets <- lapply(X = split_lung_singlets, 
                     FUN = SCTransform, 
                     method = "glmGamPoi", 
                     return.only.var.genes = FALSE,
                     vars.to.regress = c("percent.mito","S.Score","G2M.Score","Blood.Enriched1"))

saveRDS(split_lung_singlets, file = "C:/Users/Kevin/Dropbox/KAM R Code/SCT_Harmony_11062022/split_lung_singlets.rds")
var.features <- SelectIntegrationFeatures(object.list = split_lung_singlets, nfeatures = 3000)                     
lung_sct_singlets <- merge(x = split_lung_singlets[[1]], y = split_lung_singlets[2:length(split_lung_singlets)], merge.data=TRUE)
VariableFeatures(lung_sct_singlets) <- var.features
lung_sct_singlets <- RunPCA(lung_sct_singlets, verbose = FALSE)
lung_sct_singlets<- RunHarmony(lung_sct_singlets, assay.use="SCT", group.by.vars = "orig.ident", plot_convergence = TRUE)
saveRDS(lung_sct_singlets, file = "C:/Users/Kevin/Dropbox/KAM R Code/SCT_Harmony_11062022/lung_scth_singlets.rds")
lung_sct_singlets@reductions
#convergence plot looks ok despite quick trans error, so moveforward and check if have embeddings
#didn't alter standard parameters as didn't appear to make a difference, also get convergence after 7 iterations no matter what
harmony_embeddings <- Embeddings(lung_sct_singlets, 'harmony')
harmony_embeddings[1:5, 1:5]

#visualizing if batch effects are gone, appear to be
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = lung_sct_singlets, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = lung_sct_singlets, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)
#making umap
lung_sct_singlets <- RunUMAP(lung_sct_singlets, reduction = "harmony", dims = 1:30)
lung_sct_singlets <- FindNeighbors(lung_sct_singlets, reduction = "harmony", dims = 1:30) %>% FindClusters()
lung_sct_singletsDP<-DimPlot(lung_sct_singlets, group.by = "orig.ident")
lung_sct_singletsDP

single_clusters <- DimPlot(lung_sct_singlets, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
single_condition <- DimPlot(lung_sct_singlets, reduction = 'umap', group.by = 'orig.ident')
single_condition|single_clusters
#save
saveRDS(lung_sct_singlets, file = "C:/Users/Kevin/Dropbox/KAM R Code/SCT_Harmony_11062022/lung_scth_singlets.rds")

#loading data for new session
lung_scth_singlets <- readRDS("C:/Users/Kevin/Dropbox/KAM R Code/SCT_Harmony_11062022/lung_scth_singlets.rds")
#visualize
single_clusters <- DimPlot(lung_scth_singlets, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
single_condition <- DimPlot(lung_scth_singlets, reduction = 'umap', group.by = 'orig.ident')
single_condition|single_clusters

# Now can also identify cell markers for harmonized, singlet cell clusters to help identify too
# Use FindConservedMarkers here because, since we have cells from all time points in each group
# want to find what is conserved over each time point to identify the cluster
# NOT what is potentially differential expressed,so don't want to use #FindAllMarkers
# Notes:
# slot depends on the type of the test used, 
# default is data slot that stores normalized data
# be sure to check DefaultAssay(ifnb_harmony) <- 'RNA'
DefaultAssay(lung_scth_singlets)
#checking markers
markers_cluster0 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 0,grouping.var = 'orig.ident')
markers_cluster1 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 1,grouping.var = 'orig.ident')
markers_cluster2 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 2,grouping.var = 'orig.ident')
markers_cluster3 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 3,grouping.var = 'orig.ident')
markers_cluster4 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 4,grouping.var = 'orig.ident')
markers_cluster5 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 5,grouping.var = 'orig.ident')
markers_cluster6 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 6,grouping.var = 'orig.ident')
markers_cluster7 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 7,grouping.var = 'orig.ident')
markers_cluster8 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 8,grouping.var = 'orig.ident')
markers_cluster9 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 9,grouping.var = 'orig.ident')
markers_cluster10 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 10,grouping.var = 'orig.ident')
markers_cluster11 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 11,grouping.var = 'orig.ident')
markers_cluster12 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 12,grouping.var = 'orig.ident')
markers_cluster13 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 13,grouping.var = 'orig.ident')
markers_cluster14 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 14,grouping.var = 'orig.ident')
markers_cluster15 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 15,grouping.var = 'orig.ident')
markers_cluster16 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 16,grouping.var = 'orig.ident')
markers_cluster17 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 17,grouping.var = 'orig.ident')
markers_cluster18 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 18,grouping.var = 'orig.ident')
markers_cluster19 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 19,grouping.var = 'orig.ident')
markers_cluster20 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 20,grouping.var = 'orig.ident')
markers_cluster21 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 21,grouping.var = 'orig.ident')
markers_cluster22 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 22,grouping.var = 'orig.ident')
markers_cluster23 <- FindConservedMarkers(lung_scth_singlets,ident.1 = 23,grouping.var = 'orig.ident')

#see maker lists and stats
head(markers_cluster0)
write.table(markers_cluster0, file = 'markers_cluster0.csv', sep =',')
head(markers_cluster1)
write.table(markers_cluster1, file = 'markers_cluster1.csv', sep =',')
head(markers_cluster2)
write.table(markers_cluster2, file = 'markers_cluster2.csv', sep =',')
head(markers_cluster3)
write.table(markers_cluster3, file = 'markers_cluster3.csv', sep =',')
head(markers_cluster4)
write.table(markers_cluster4, file = 'markers_cluster4.csv', sep =',')
head(markers_cluster5)
write.table(markers_cluster5, file = 'markers_cluster5.csv', sep =',')
head(markers_cluster6)
write.table(markers_cluster6, file = 'markers_cluster6.csv', sep =',')
head(markers_cluster7)
write.table(markers_cluster7, file = 'markers_cluster7.csv', sep =',')
head(markers_cluster8)
write.table(markers_cluster8, file = 'markers_cluster8.csv', sep =',')
head(markers_cluster9)
write.table(markers_cluster9, file = 'markers_cluster9.csv', sep =',')
head(markers_cluster10)
write.table(markers_cluster10, file = 'markers_cluster10.csv', sep =',')
head(markers_cluster11)
write.table(markers_cluster11, file = 'markers_cluster11.csv', sep =',')
head(markers_cluster12)
write.table(markers_cluster12, file = 'markers_cluster12.csv', sep =',')
head(markers_cluster13)
write.table(markers_cluster13, file = 'markers_cluster13.csv', sep =',')
head(markers_cluster14)
write.table(markers_cluster14, file = 'markers_cluster14.csv', sep =',')
head(markers_cluster15)
write.table(markers_cluster15, file = 'markers_cluster15.csv', sep =',')
head(markers_cluster16)
write.table(markers_cluster16, file = 'markers_cluster16.csv', sep =',')
head(markers_cluster17)
write.table(markers_cluster17, file = 'markers_cluster17.csv', sep =',')
head(markers_cluster18)
write.table(markers_cluster18, file = 'markers_cluster18.csv', sep =',')
head(markers_cluster19)
write.table(markers_cluster19, file = 'markers_cluster19.csv', sep =',')
head(markers_cluster20)
write.table(markers_cluster20, file = 'markers_cluster20.csv', sep =',')
head(markers_cluster21)
write.table(markers_cluster21, file = 'markers_cluster21.csv', sep =',')
head(markers_cluster22)
write.table(markers_cluster22, file = 'markers_cluster22.csv', sep =',')
head(markers_cluster23)
write.table(markers_cluster23, file = 'markers_cluster23.csv', sep =',')

#Some of the clusters didn't return conserved markers so will also run FindALLMarkers to see cluster markers
#following this tutorial (https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_integ_marker_identification.html)
All_markers <- FindAllMarkers(object = lung_scth_singlets, only.pos = TRUE,logfc.threshold = 0.25, min.pct = )
#Reording output columns to be easier to read; putting cluster and gene first instead of last
# Combine markers with gene descriptions
library(dplyr)
All_markers  <- All_markers [ , c(6, 7, 2:4, 1, 5)]

# Order the rows by p-adjusted values
All_markers<- All_markers %>%
  dplyr::arrange(cluster, p_val_adj)

View(All_markers)
write.table(All_markers, file = 'All_markers.csv', sep =',')

#noticed that there are still quite a few blood cells despite filtering out in SCT, so
# Adding Score for blood cells
blood_genes_to_check = list(c("Hbb-bt","Hbb-bs", "Hba-a1", "Hba-a2", "Hba-2", "Hbb", "Hba-1", "Hbg2"))
lung_sct_singlets <-AddModuleScore(lung_sct_singlets, feature = blood_genes_to_check, name = "Blood.Enriched")
blood_FP <- FeaturePlot(lung_scth_singlets,features = "Blood.Enriched1", pt.size = 0.1, label = TRUE ,label.color = "violetred2", label.size = 4, repel = TRUE)
blood_FP
blood_VP <- VlnPlot(lung_scth_singlets, features = "Blood.Enriched1")
blood_VP

# and excluding these cells 
lung <- subset(lung_scth_singlets, subset = Blood.Enriched1 <= 2.5)
DimPlot(lung)
 #re-checking markers to see if all RBC contamination is gone
markers_cluster0 <- FindConservedMarkers(lung,ident.1 = 0,grouping.var = 'orig.ident')
markers_cluster1 <- FindConservedMarkers(lung,ident.1 = 1,grouping.var = 'orig.ident')
markers_cluster2 <- FindConservedMarkers(lung,ident.1 = 2,grouping.var = 'orig.ident')
markers_cluster3 <- FindConservedMarkers(lung,ident.1 = 3,grouping.var = 'orig.ident')
markers_cluster4 <- FindConservedMarkers(lung,ident.1 = 4,grouping.var = 'orig.ident')
markers_cluster5 <- FindConservedMarkers(lung,ident.1 = 5,grouping.var = 'orig.ident')
markers_cluster6 <- FindConservedMarkers(lung,ident.1 = 6,grouping.var = 'orig.ident')
markers_cluster7 <- FindConservedMarkers(lung,ident.1 = 7,grouping.var = 'orig.ident')
markers_cluster8 <- FindConservedMarkers(lung,ident.1 = 8,grouping.var = 'orig.ident')
markers_cluster9 <- FindConservedMarkers(lung,ident.1 = 9,grouping.var = 'orig.ident')
markers_cluster10 <- FindConservedMarkers(lung,ident.1 = 10,grouping.var = 'orig.ident')
markers_cluster11 <- FindConservedMarkers(lung,ident.1 = 11,grouping.var = 'orig.ident')
markers_cluster12 <- FindConservedMarkers(lung,ident.1 = 12,grouping.var = 'orig.ident')
markers_cluster13 <- FindConservedMarkers(lung,ident.1 = 13,grouping.var = 'orig.ident')
markers_cluster14 <- FindConservedMarkers(lung,ident.1 = 14,grouping.var = 'orig.ident')
markers_cluster15 <- FindConservedMarkers(lung,ident.1 = 15,grouping.var = 'orig.ident')
markers_cluster16 <- FindConservedMarkers(lung,ident.1 = 16,grouping.var = 'orig.ident')
markers_cluster17 <- FindConservedMarkers(lung,ident.1 = 17,grouping.var = 'orig.ident')
markers_cluster18 <- FindConservedMarkers(lung,ident.1 = 18,grouping.var = 'orig.ident')
markers_cluster19 <- FindConservedMarkers(lung,ident.1 = 19,grouping.var = 'orig.ident')
markers_cluster20 <- FindConservedMarkers(lung,ident.1 = 20,grouping.var = 'orig.ident')
markers_cluster21 <- FindConservedMarkers(lung,ident.1 = 21,grouping.var = 'orig.ident')
markers_cluster22 <- FindConservedMarkers(lung,ident.1 = 22,grouping.var = 'orig.ident')
markers_cluster23 <- FindConservedMarkers(lung,ident.1 = 23,grouping.var = 'orig.ident')

#looking at markers still show RBC contamination, so will just subset out those clusters (12,19, 20, 21, 22) and re-run normalization and clustering
Idents(lung)<- 'seurat_clusters'
levels(lung)
this_cluster <-c('0','1','2','3','4','5','6','7','8','9','10','11','13','14','15','16','17','18','23')
lung <-subset(lung, idents = this_cluster)
saveRDS(lung, file = "C:/Users/Kevin/Dropbox/KAM R Code/lung.rds")

#re-running normalization and harmony
Idents(lung)<- 'orig.ident'
levels(lung)
split_lung <- SplitObject(lung, split.by = "orig.ident")
split_lung <- lapply(X = split_lung, 
                              FUN = SCTransform, 
                              method = "glmGamPoi", 
                              return.only.var.genes = FALSE,
                              vars.to.regress = c("percent.mito","S.Score","G2M.Score"))


var.features <- SelectIntegrationFeatures(object.list = split_lung, nfeatures = 3000)                     
lung<- merge(x = split_lung[[1]], y = split_lung[2:length(split_lung)], merge.data=TRUE)
VariableFeatures(lung) <- var.features
lung <- RunPCA(lung, verbose = FALSE)
lung<- RunHarmony(lung, assay.use="SCT", group.by.vars = "orig.ident", plot_convergence = TRUE)


#visualizing if batch effects are gone, appear to be
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = lung, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = lung, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)
#making umap
lung <- RunUMAP(lung, reduction = "harmony", dims = 1:30)
lung <- FindNeighbors(lung, reduction = "harmony", dims = 1:30) %>% FindClusters()
lungDP<-DimPlot(lung, group.by = "orig.ident")
lungDP

single_clusters <- DimPlot(lung, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
single_condition <- DimPlot(lung, reduction = 'umap', group.by = 'orig.ident')
single_condition|single_clusters

#save file of singlets
saveRDS(lung, file = "C:/Users/Kevin/Dropbox/KAM R Code/lung.rds")

# Use FindConservedMarkers here because, since we have cells from all time points in each group
# want to find what is conserved over each time point to identify the cluster
# be sure to check DefaultAssay(ifnb_harmony) <- 'RNA'
DefaultAssay(lung) <- "RNA"
DefaultAssay(lung)

#now run cluster maker identification
markers_cluster0 <- FindConservedMarkers(lung,ident.1 = 0,grouping.var = 'orig.ident')
markers_cluster1 <- FindConservedMarkers(lung,ident.1 = 1,grouping.var = 'orig.ident')
markers_cluster2 <- FindConservedMarkers(lung,ident.1 = 2,grouping.var = 'orig.ident')
markers_cluster3 <- FindConservedMarkers(lung,ident.1 = 3,grouping.var = 'orig.ident')
markers_cluster4 <- FindConservedMarkers(lung,ident.1 = 4,grouping.var = 'orig.ident')
markers_cluster5 <- FindConservedMarkers(lung,ident.1 = 5,grouping.var = 'orig.ident')
markers_cluster6 <- FindConservedMarkers(lung,ident.1 = 6,grouping.var = 'orig.ident')
markers_cluster7 <- FindConservedMarkers(lung,ident.1 = 7,grouping.var = 'orig.ident')
markers_cluster8 <- FindConservedMarkers(lung,ident.1 = 8,grouping.var = 'orig.ident')
markers_cluster9 <- FindConservedMarkers(lung,ident.1 = 9,grouping.var = 'orig.ident')
markers_cluster10 <- FindConservedMarkers(lung,ident.1 = 10,grouping.var = 'orig.ident')
markers_cluster11 <- FindConservedMarkers(lung,ident.1 = 11,grouping.var = 'orig.ident')
markers_cluster12 <- FindConservedMarkers(lung,ident.1 = 12,grouping.var = 'orig.ident')
markers_cluster13 <- FindConservedMarkers(lung,ident.1 = 13,grouping.var = 'orig.ident')
markers_cluster14 <- FindConservedMarkers(lung,ident.1 = 14,grouping.var = 'orig.ident')
markers_cluster15 <- FindConservedMarkers(lung,ident.1 = 15,grouping.var = 'orig.ident')
markers_cluster16 <- FindConservedMarkers(lung,ident.1 = 16,grouping.var = 'orig.ident')
markers_cluster17 <- FindConservedMarkers(lung,ident.1 = 17,grouping.var = 'orig.ident')
markers_cluster18 <- FindConservedMarkers(lung,ident.1 = 18,grouping.var = 'orig.ident')
markers_cluster19 <- FindConservedMarkers(lung,ident.1 = 19,grouping.var = 'orig.ident')
markers_cluster20 <- FindConservedMarkers(lung,ident.1 = 20,grouping.var = 'orig.ident')

#see maker lists and stats
head(markers_cluster0)
write.table(markers_cluster0, file = 'markers_cluster0.csv', sep =',')
head(markers_cluster1)
write.table(markers_cluster1, file = 'markers_cluster1.csv', sep =',')
head(markers_cluster2)
write.table(markers_cluster2, file = 'markers_cluster2.csv', sep =',')
head(markers_cluster3)
write.table(markers_cluster3, file = 'markers_cluster3.csv', sep =',')
head(markers_cluster4)
write.table(markers_cluster4, file = 'markers_cluster4.csv', sep =',')
head(markers_cluster5)
write.table(markers_cluster5, file = 'markers_cluster5.csv', sep =',')
head(markers_cluster6)
write.table(markers_cluster6, file = 'markers_cluster6.csv', sep =',')
head(markers_cluster7)
write.table(markers_cluster7, file = 'markers_cluster7.csv', sep =',')
head(markers_cluster8)
write.table(markers_cluster8, file = 'markers_cluster8.csv', sep =',')
head(markers_cluster9)
write.table(markers_cluster9, file = 'markers_cluster9.csv', sep =',')
head(markers_cluster10)
write.table(markers_cluster10, file = 'markers_cluster10.csv', sep =',')
head(markers_cluster11)
write.table(markers_cluster11, file = 'markers_cluster11.csv', sep =',')
head(markers_cluster12)
write.table(markers_cluster12, file = 'markers_cluster12.csv', sep =',')
head(markers_cluster13)
write.table(markers_cluster13, file = 'markers_cluster13.csv', sep =',')
head(markers_cluster14)
write.table(markers_cluster14, file = 'markers_cluster14.csv', sep =',')
head(markers_cluster15)
write.table(markers_cluster15, file = 'markers_cluster15.csv', sep =',')
head(markers_cluster16)
write.table(markers_cluster16, file = 'markers_cluster16.csv', sep =',')
head(markers_cluster17)
write.table(markers_cluster17, file = 'markers_cluster17.csv', sep =',')
head(markers_cluster18)
write.table(markers_cluster18, file = 'markers_cluster18.csv', sep =',')
head(markers_cluster19)
write.table(markers_cluster19, file = 'markers_cluster19.csv', sep =',')
head(markers_cluster20)
write.table(markers_cluster20, file = 'markers_cluster20.csv', sep =',')

#Some of the clusters didn't return conserved markers so will also run FindALLMarkers to see cluster markers
#following this tutorial (https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_integ_marker_identification.html)
All_markers <- FindAllMarkers(object = lung, only.pos = TRUE,logfc.threshold = 0.25, min.pct = )

#Reording output columns to be easier to read; putting cluster and gene first instead of last
# Combine markers with gene descriptions
library(dplyr)
All_markers  <- All_markers [ , c(6, 7, 2:4, 1, 5)]

# Order the rows by p-adjusted values
All_markers<- All_markers %>%
  dplyr::arrange(cluster, p_val_adj)

View(All_markers)
write.table(All_markers, file = 'All_markers.csv', sep =',')

#Some groups oddly didn't have conserved markers, so checking FindMarkers too
findmarkers_cluster0 <- FindMarkers(lung,ident.1 = 0, logfc.threshold = 0.1)
write.table(findmarkers_cluster0, file = 'findmarkers_cluster0.csv', sep =',')

findmarkers_cluster2 <- FindMarkers(lung,ident.1 = 2, logfc.threshold = 0.1)
write.table(findmarkers_cluster2, file = 'findmarkers_cluster2.csv', sep =',')

findmarkers_cluster19 <- FindMarkers(lung,ident.1 = 19, logfc.threshold = 0.1)
write.table(findmarkers_cluster19, file = 'findmarkers_cluster19.csv', sep =',')

findmarkers_cluster20 <- FindMarkers(lung,ident.1 = 20, logfc.threshold = 0.1)
write.table(findmarkers_cluster20, file = 'findmarkers_cluster20.csv', sep =',')

# looks like clusters 19 and 20 are perhaps singletons, as with table can see there are only 2 cells in each
Idents(lung)<- 'seurat_clusters'
table(Idents(lung))
#thus, exclusing from cluster (not re-running all analysis as only 4 cells)
this_cluster <-c('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18')
lung <-subset(lung, idents = this_cluster)
lungDimP <-DimPlot(lung, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
lungDimP
saveRDS(lung, file = "C:/Users/Kevin/Dropbox/KAM R Code/lung.rds")

# need to subset clusters 0 and 2 to see if can better identify cell types
Idents(lung)<- 'seurat_clusters'
this_cluster <-'0'
lung_zero <- subset(lung, idents = this_cluster)

Idents(lung_zero)<- 'orig.ident'
levels(lung_zero)
DefaultAssay(lung_zero) <- "RNA"

split_lung_zero <- SplitObject(lung_zero, split.by = "orig.ident")
split_lung_zero <- lapply(X = split_lung_zero, 
                     FUN = SCTransform, 
                     method = "glmGamPoi", 
                     return.only.var.genes = FALSE,
                     vars.to.regress = c("percent.mito","S.Score","G2M.Score"))


var.features <- SelectIntegrationFeatures(object.list = split_lung_zero, nfeatures = 3000)                     
lung_zero <- merge(x = split_lung_zero[[1]], y = split_lung_zero[2:length(split_lung_zero)], merge.data=TRUE)
VariableFeatures(lung_zero) <- var.features
lung_zero <- RunPCA(lung_zero, verbose = FALSE)
# will check to see if need harmony here as doesn't converge if run on this subset alone
# think already data is too similar to adjust again since it is a subcluster 
lung_zero <- RunUMAP(lung_zero, reduction = "pca", dims = 1:30)
lung_zero <- FindNeighbors(lung_zero, reduction = "pca", dims = 1:30) 
lung_zero <- FindClusters(lung_zero, resolution = c(0.1, 0.2, 0.3, 0.5, 0.6, 0.7,0.8, 1))
z1 <-DimPlot(lung_zero, group.by = "SCT_snn_res.0.1", label =TRUE)
z2 <-DimPlot(lung_zero, group.by = "SCT_snn_res.0.2", label =TRUE)
z3 <-DimPlot(lung_zero, group.by = "SCT_snn_res.0.3", label =TRUE)
z5 <-DimPlot(lung_zero, group.by = "SCT_snn_res.0.5", label =TRUE)
z6 <-DimPlot(lung_zero, group.by = "SCT_snn_res.0.6", label =TRUE)
z7 <-DimPlot(lung_zero, group.by = "SCT_snn_res.0.7", label =TRUE)
z8 <-DimPlot(lung_zero, group.by = "SCT_snn_res.0.8", label =TRUE)
z1.1 <-DimPlot(lung_zero, group.by = "SCT_snn_res.1", label =TRUE)
z1|z2|z3
z5|z6|z7
z8|z1.1
#decided to go with 0.8 resolution
lung_zero <- FindClusters(lung_zero, resolution = 0.8)

lung_zeroDP<-DimPlot(lung_zero, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
lung_zeroDP

#check for batch effects
#looks ok, so will not run harmony again
DimPlot(lung_zero, reduction = 'umap', group.by = 'orig.ident', label = FALSE)

#now identify lung_ zero clusters
DefaultAssay(lung_zero) <- "RNA"
DefaultAssay(lung_zero)

#now run cluster maker identification
Zmarkers_cluster0 <- FindConservedMarkers(lung_zero,ident.1 = 0,grouping.var = 'orig.ident')
head(Zmarkers_cluster0)
write.table(Zmarkers_cluster0, file = 'Zmarkers_cluster0.csv', sep =',')

Zmarkers_cluster1 <- FindConservedMarkers(lung_zero,ident.1 = 1,grouping.var = 'orig.ident')
head(Zmarkers_cluster1)
write.table(Zmarkers_cluster1, file = 'Zmarkers_cluster1.csv', sep =',')

Zmarkers_cluster2 <- FindConservedMarkers(lung_zero,ident.1 = 2,grouping.var = 'orig.ident')
head(Zmarkers_cluster2)
write.table(Zmarkers_cluster2, file = 'Zmarkers_cluster2.csv', sep =',')

Zmarkers_cluster3 <- FindConservedMarkers(lung_zero,ident.1 = 3,grouping.var = 'orig.ident')
head(Zmarkers_cluster3)
write.table(Zmarkers_cluster3, file = 'Zmarkers_cluster3.csv', sep =',')

Zmarkers_cluster4 <- FindConservedMarkers(lung_zero,ident.1 = 4,grouping.var = 'orig.ident')
head(Zmarkers_cluster4)
write.table(Zmarkers_cluster4, file = 'Zmarkers_cluster4.csv', sep =',')

Zmarkers_cluster5 <- FindConservedMarkers(lung_zero,ident.1 = 5,grouping.var = 'orig.ident')
head(Zmarkers_cluster5)
write.table(Zmarkers_cluster5, file = 'Zmarkers_cluster5.csv', sep =',')

Zmarkers_cluster6 <- FindConservedMarkers(lung_zero,ident.1 = 6,grouping.var = 'orig.ident')
head(Zmarkers_cluster6)
write.table(Zmarkers_cluster0, file = 'Zmarkers_cluster0.csv', sep =',')

Zmarkers_cluster0 <- FindConservedMarkers(lung_zero,ident.1 = 0,grouping.var = 'orig.ident')
head(Zmarkers_cluster0)
write.table(Zmarkers_cluster6, file = 'Zmarkers_cluster6.csv', sep =',')

Zmarkers_cluster7 <- FindConservedMarkers(lung_zero,ident.1 = 7,grouping.var = 'orig.ident')
head(Zmarkers_cluster7)
write.table(Zmarkers_cluster7, file = 'Zmarkers_cluster7.csv', sep =',')

Zmarkers_cluster8 <- FindConservedMarkers(lung_zero,ident.1 = 8,grouping.var = 'orig.ident')
head(Zmarkers_cluster8)
write.table(Zmarkers_cluster8, file = 'Zmarkers_cluster8.csv', sep =',')

# Checking all makers 
DefaultAssay(lung_zero) <- "RNA"
DefaultAssay(lung_zero)

All_Zmarkers <- FindAllMarkers(object = lung_zero, only.pos = TRUE,logfc.threshold = 0.1, min.pct = )
library(dplyr)
All_Zmarkers  <- All_Zmarkers [ , c(6, 7, 2:4, 1, 5)]
All_Zmarkers<- All_Zmarkers %>%
  dplyr::arrange(cluster, p_val_adj)

View(All_Zmarkers)
write.table(All_Zmarkers, file = 'All_Zmarkers.csv', sep =',')

Idents(lung_zero)<- 'seurat_clusters'
table(Idents(lung_zero))

cluster0.zmarkers <- FindMarkers(lung_zero, ident.1 = 0, min.pct = 0.25)
head(cluster0.zmarkers, n = 10)
#Renaming Idents to identify cell clusters
lung_zero <- RenameIdents(lung_zero, "0" = "Degraded")
lung_zero <- RenameIdents(lung_zero, "1" = "Monocytes")
lung_zero <- RenameIdents(lung_zero, "2" = "T Cells")
lung_zero <- RenameIdents(lung_zero, "3" = "T Cells")
lung_zero <- RenameIdents(lung_zero, "4" = "T Cells")
lung_zero <- RenameIdents(lung_zero, "5" = "Neutrophils")
lung_zero <- RenameIdents(lung_zero, "6" = "Monocytes")
lung_zero <- RenameIdents(lung_zero, "7" = "B Cells")
lung_zero <- RenameIdents(lung_zero, "8" = "A Cap")

# visualize newly labelled Dim Plit
lung_zero_LabelledDimP <- DimPlot(lung_zero, label = TRUE)
lung_zero_LabelledDimP

# can now stash these identities as a new column in the meta data
lung_zero <- StashIdent(lung_zero, save.name = "subcluster.label")
#then re-assign back to lung object
lung$subcluster.label <- as.character(Idents(lung))
lung$subcluster.label[Cells(lung_zero)] <- paste("Subcluster0_",Idents(lung_zero))
View(lung@meta.data)
Idents(lung) <- "subcluster.label"
lung_LabelledDimP <- DimPlot(lung, label = TRUE)
lung_LabelledDimP

# save file with new idents & save subcluster analysis object
saveRDS(lung, file = "C:/Users/Kevin/Dropbox/KAM R Code/lung.rds")
saveRDS(lung_zero, file = "C:/Users/Kevin/Dropbox/KAM R Code/lung_zero.rds")

##CLUSTER TWO
#need to subset cluster 2 to get cell types as well
Idents(lung)<- 'seurat_clusters'
this_cluster <-'2'
lung_two <- subset(lung, idents = this_cluster)

Idents(lung_two)<- 'orig.ident'
levels(lung_two)
DefaultAssay(lung_two) <- "RNA"

split_lung_two <- SplitObject(lung_two, split.by = "orig.ident")
split_lung_two <- lapply(X = split_lung_two, 
                          FUN = SCTransform, 
                          method = "glmGamPoi", 
                          return.only.var.genes = FALSE,
                          vars.to.regress = c("percent.mito","S.Score","G2M.Score"))


var.features <- SelectIntegrationFeatures(object.list = split_lung_two, nfeatures = 3000)                     
lung_two <- merge(x = split_lung_two[[1]], y = split_lung_two[2:length(split_lung_two)], merge.data=TRUE)
VariableFeatures(lung_two) <- var.features
lung_two <- RunPCA(lung_two, verbose = FALSE)
# will check to see if need harmony here as doesn't converge if run on this subset alone
# think already data is too similar to adjust again since it is a subcluster 
lung_two <- RunUMAP(lung_two, reduction = "pca", dims = 1:30)
lung_two <- FindNeighbors(lung_two, reduction = "pca", dims = 1:30) 
lung_two <- FindClusters(lung_two, resolution = c(0.1, 0.2, 0.3, 0.5, 0.6, 0.7,0.8, 1))
t1 <-DimPlot(lung_two, group.by = "SCT_snn_res.0.1", label =TRUE)
t2 <-DimPlot(lung_two, group.by = "SCT_snn_res.0.2", label =TRUE)
t3 <-DimPlot(lung_two, group.by = "SCT_snn_res.0.3", label =TRUE)
t5 <-DimPlot(lung_two, group.by = "SCT_snn_res.0.5", label =TRUE)
t6 <-DimPlot(lung_two, group.by = "SCT_snn_res.0.6", label =TRUE)
t7 <-DimPlot(lung_two, group.by = "SCT_snn_res.0.7", label =TRUE)
t8 <-DimPlot(lung_two, group.by = "SCT_snn_res.0.8", label =TRUE)
t1.1 <-DimPlot(lung_two, group.by = "SCT_snn_res.1", label =TRUE)
t1|t2|t3
t5|t6|t7
t8|t1.1
#decided to go with 0.7 resolution
lung_two <- FindClusters(lung_two, resolution = 0.7)

lung_twoDP<-DimPlot(lung_two, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
lung_twoDP

#check for batch effects
#looks ok, so will not run harmony again
DimPlot(lung_two, reduction = 'umap', group.by = 'orig.ident', label = FALSE)

#now identify lung_two clusters
DefaultAssay(lung_two) <- "RNA"
DefaultAssay(lung_two)

#now run cluster maker identification
Tmarkers_cluster0 <- FindConservedMarkers(lung_two,ident.1 = 0,grouping.var = 'orig.ident')
head(Tmarkers_cluster0)
write.table(Tmarkers_cluster0, file = 'Tmarkers_cluster0.csv', sep =',')

Tmarkers_cluster1 <- FindConservedMarkers(lung_two,ident.1 = 1,grouping.var = 'orig.ident')
head(Tmarkers_cluster1)
write.table(Tmarkers_cluster1, file = 'Tmarkers_cluster1.csv', sep =',')

Tmarkers_cluster2 <- FindConservedMarkers(lung_two,ident.1 = 2,grouping.var = 'orig.ident')
head(Tmarkers_cluster2)
write.table(Tmarkers_cluster2, file = 'Tmarkers_cluster2.csv', sep =',')

Tmarkers_cluster3 <- FindConservedMarkers(lung_two,ident.1 = 3,grouping.var = 'orig.ident')
head(Tmarkers_cluster3)
write.table(Tmarkers_cluster3, file = 'Tmarkers_cluster3.csv', sep =',')

Tmarkers_cluster4 <- FindConservedMarkers(lung_two,ident.1 = 4,grouping.var = 'orig.ident')
head(Tmarkers_cluster4)
write.table(Tmarkers_cluster4, file = 'Zmarkers_cluster4.csv', sep =',')

Tmarkers_cluster5 <- FindConservedMarkers(lung_two,ident.1 = 5,grouping.var = 'orig.ident')
head(Tmarkers_cluster5)
write.table(Tmarkers_cluster5, file = 'Zmarkers_cluster5.csv', sep =',')


# Checking all makers 
DefaultAssay(lung_two) <- "RNA"
DefaultAssay(lung_two)

All_Tmarkers <- FindAllMarkers(object = lung_two, only.pos = TRUE,logfc.threshold = 0.1)
library(dplyr)
All_Tmarkers  <- All_Tmarkers [ , c(6, 7, 2:4, 1, 5)]
All_Tmarkers<- All_Tmarkers %>%
  dplyr::arrange(cluster, p_val_adj)

View(All_Tmarkers)
write.table(All_Tmarkers, file = 'All_Tmarkers.csv', sep =',')

Idents(lung_two)<- 'seurat_clusters'
table(Idents(lung_two))

#Renaming Idents to identify cell clusters
lung_two <- RenameIdents(lung_two, "0" = "T cells")
lung_two <- RenameIdents(lung_two, "1" = "Endothelial")
lung_two <- RenameIdents(lung_two, "2" = "Degraded")
lung_two <- RenameIdents(lung_two, "3" = "Endothelial")
lung_two <- RenameIdents(lung_two, "4" = "B Cells")
lung_two <- RenameIdents(lung_two, "5" = "Monocytes")

# visualize newly labelled Dim Plot
lung_two_LabelledDimP <- DimPlot(lung_two, label = TRUE)
lung_two_LabelledDimP

# can now stash these identities as a new column in the meta data
lung_two <- StashIdent(lung_two, save.name = "subcluster.label")
#then re-assign back to lung object
Idents(lung)<- 'subcluster.label'
lung$subcluster.label <- as.character(Idents(lung))
lung$subcluster.label[Cells(lung_two)] <- paste("Subcluster2_",Idents(lung_two))
View(lung@meta.data)
Idents(lung) <- "subcluster.label"
lung_LabelledDimP <- DimPlot(lung, label = TRUE)
lung_LabelledDimP

# save file with new idents & save subcluster analysis object
saveRDS(lung, file = "C:/Users/Katherine/Dropbox/KAM R Code/lung.rds")
saveRDS(lung_two, file = "C:/Users/Katherine/Dropbox/KAM R Code/lung_two.rds")

# REMOVING DEGRADED CELLS & RECLUSTERING
# need to remove degraded cells from subcluster analysis of clusters 0 and 2. 
Idents(lung) <- "subcluster.label"
this_cluster <-c("10","Subcluster0_ Neutrophils","5","17","Subcluster0_ Monocytes","13","1","Subcluster0_ T Cells",
                 "Subcluster2_ Endothelial","3","8","9","4","7","12","15","11","Subcluster2_ T cells","6","16","14","18",
                 "Subcluster0_ B Cells","Subcluster0_ A Cap","Subcluster2_ B Cells","Subcluster2_ Monocytes")
lung <-subset(lung, idents = this_cluster)

#re-running normalization and harmony
Idents(lung)<- 'orig.ident'
levels(lung)
split_lung <- SplitObject(lung, split.by = "orig.ident")
split_lung <- lapply(X = split_lung, 
                     FUN = SCTransform, 
                     method = "glmGamPoi", 
                     return.only.var.genes = FALSE,
                     vars.to.regress = c("percent.mito","S.Score","G2M.Score"))


var.features <- SelectIntegrationFeatures(object.list = split_lung, nfeatures = 3000)                     
lung<- merge(x = split_lung[[1]], y = split_lung[2:length(split_lung)], merge.data=TRUE)
VariableFeatures(lung) <- var.features
lung <- RunPCA(lung, verbose = FALSE)
lung<- RunHarmony(lung, assay.use="SCT", group.by.vars = "orig.ident", plot_convergence = TRUE, kmeans_init_nstart=20, kmeans_init_iter_max=100)


#visualizing if batch effects are gone, appear to be
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = lung, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = lung, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)
#making umap
lung <- RunUMAP(lung, reduction = "harmony", dims = 1:30)
lung <- FindNeighbors(lung, reduction = "harmony", dims = 1:30) %>% FindClusters()
lungDP<-DimPlot(lung, group.by = "orig.ident")
lungDP

single_clusters <- DimPlot(lung, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
single_condition <- DimPlot(lung, reduction = 'umap', group.by = 'orig.ident')
single_condition|single_clusters

#save file of singlets
saveRDS(lung, file = "C:/Users/Kevin/Dropbox/KAM R Code/lung.rds")

# Use FindConservedMarkers to identify clusters
# be sure to check DefaultAssay(ifnb_harmony) <- 'RNA'
DefaultAssay(lung) <- "RNA"
DefaultAssay(lung)

#now run cluster maker identification
markers_cluster0 <- FindConservedMarkers(lung,ident.1 = 0,grouping.var = 'orig.ident')
markers_cluster1 <- FindConservedMarkers(lung,ident.1 = 1,grouping.var = 'orig.ident')
markers_cluster2 <- FindConservedMarkers(lung,ident.1 = 2,grouping.var = 'orig.ident')
markers_cluster3 <- FindConservedMarkers(lung,ident.1 = 3,grouping.var = 'orig.ident')
markers_cluster4 <- FindConservedMarkers(lung,ident.1 = 4,grouping.var = 'orig.ident')
markers_cluster5 <- FindConservedMarkers(lung,ident.1 = 5,grouping.var = 'orig.ident')
markers_cluster6 <- FindConservedMarkers(lung,ident.1 = 6,grouping.var = 'orig.ident')
markers_cluster7 <- FindConservedMarkers(lung,ident.1 = 7,grouping.var = 'orig.ident')
markers_cluster8 <- FindConservedMarkers(lung,ident.1 = 8,grouping.var = 'orig.ident')
markers_cluster9 <- FindConservedMarkers(lung,ident.1 = 9,grouping.var = 'orig.ident')
markers_cluster10 <- FindConservedMarkers(lung,ident.1 = 10,grouping.var = 'orig.ident')
markers_cluster11 <- FindConservedMarkers(lung,ident.1 = 11,grouping.var = 'orig.ident')
markers_cluster12 <- FindConservedMarkers(lung,ident.1 = 12,grouping.var = 'orig.ident')
markers_cluster13 <- FindConservedMarkers(lung,ident.1 = 13,grouping.var = 'orig.ident')
markers_cluster14 <- FindConservedMarkers(lung,ident.1 = 14,grouping.var = 'orig.ident')
markers_cluster15 <- FindConservedMarkers(lung,ident.1 = 15,grouping.var = 'orig.ident')
markers_cluster16 <- FindConservedMarkers(lung,ident.1 = 16,grouping.var = 'orig.ident')
markers_cluster17 <- FindConservedMarkers(lung,ident.1 = 17,grouping.var = 'orig.ident')
markers_cluster18 <- FindConservedMarkers(lung,ident.1 = 18,grouping.var = 'orig.ident')


#see maker lists and stats
head(markers_cluster0)
write.table(markers_cluster0, file = 'markers_cluster0.csv', sep =',')
head(markers_cluster1)
write.table(markers_cluster1, file = 'markers_cluster1.csv', sep =',')
head(markers_cluster2)
write.table(markers_cluster2, file = 'markers_cluster2.csv', sep =',')
head(markers_cluster3)
write.table(markers_cluster3, file = 'markers_cluster3.csv', sep =',')
head(markers_cluster4)
write.table(markers_cluster4, file = 'markers_cluster4.csv', sep =',')
head(markers_cluster5)
write.table(markers_cluster5, file = 'markers_cluster5.csv', sep =',')
head(markers_cluster6)
write.table(markers_cluster6, file = 'markers_cluster6.csv', sep =',')
head(markers_cluster7)
write.table(markers_cluster7, file = 'markers_cluster7.csv', sep =',')
head(markers_cluster8)
write.table(markers_cluster8, file = 'markers_cluster8.csv', sep =',')
head(markers_cluster9)
write.table(markers_cluster9, file = 'markers_cluster9.csv', sep =',')
head(markers_cluster10)
write.table(markers_cluster10, file = 'markers_cluster10.csv', sep =',')
head(markers_cluster11)
write.table(markers_cluster11, file = 'markers_cluster11.csv', sep =',')
head(markers_cluster12)
write.table(markers_cluster12, file = 'markers_cluster12.csv', sep =',')
head(markers_cluster13)
write.table(markers_cluster13, file = 'markers_cluster13.csv', sep =',')
head(markers_cluster14)
write.table(markers_cluster14, file = 'markers_cluster14.csv', sep =',')
head(markers_cluster15)
write.table(markers_cluster15, file = 'markers_cluster15.csv', sep =',')
head(markers_cluster16)
write.table(markers_cluster16, file = 'markers_cluster16.csv', sep =',')
head(markers_cluster17)
write.table(markers_cluster17, file = 'markers_cluster17.csv', sep =',')
head(markers_cluster18)
write.table(markers_cluster18, file = 'markers_cluster18.csv', sep =',')


#Some of the clusters didn't return conserved markers so will also run FindALLMarkers to see cluster markers
#following this tutorial (https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_integ_marker_identification.html)
All_markers <- FindAllMarkers(object = lung, only.pos = TRUE,logfc.threshold = 0.25, min.pct = )

#Reording output columns to be easier to read; putting cluster and gene first instead of last
# Combine markers with gene descriptions
library(dplyr)
All_markers  <- All_markers [ , c(6, 7, 2:4, 1, 5)]

# Order the rows by p-adjusted values
All_markers<- All_markers %>%
  dplyr::arrange(cluster, p_val_adj)

View(All_markers)
write.table(All_markers, file = 'All_labelledmarkers.csv', sep =',')

saveRDS(lung, file = "C:/Users/Kevin/Dropbox/KAM R Code/SCT_12122022/lung.rds")

##Cluster 0 looks like  mix of endothelial cells, so going to subcluster
Idents(lung)<- 'seurat_clusters'
this_cluster <-'0'
lung_zero <- subset(lung, idents = this_cluster)

Idents(lung_zero)<- 'orig.ident'
levels(lung_zero)
DefaultAssay(lung_zero) <- "RNA"

split_lung_zero <- SplitObject(lung_zero, split.by = "orig.ident")
split_lung_zero <- lapply(X = split_lung_zero, 
                          FUN = SCTransform, 
                          method = "glmGamPoi", 
                          return.only.var.genes = FALSE,
                          vars.to.regress = c("percent.mito","S.Score","G2M.Score"))


var.features <- SelectIntegrationFeatures(object.list = split_lung_zero, nfeatures = 3000)                     
lung_zero <- merge(x = split_lung_zero[[1]], y = split_lung_zero[2:length(split_lung_zero)], merge.data=TRUE)
VariableFeatures(lung_zero) <- var.features
lung_zero <- RunPCA(lung_zero, verbose = FALSE)
# will check to see if need harmony here as doesn't converge if run on this subset alone
# think already data is too similar to adjust again since it is a subcluster 
lung_zero <- RunUMAP(lung_zero, reduction = "pca", dims = 1:30)
lung_zero <- FindNeighbors(lung_zero, reduction = "pca", dims = 1:30) 
lung_zero <- FindClusters(lung_zero, resolution = c(0.1, 0.2, 0.3, 0.5, 0.6, 0.7,0.8, 1))
z1 <-DimPlot(lung_zero, group.by = "SCT_snn_res.0.1", label =TRUE)
z2 <-DimPlot(lung_zero, group.by = "SCT_snn_res.0.2", label =TRUE)
z3 <-DimPlot(lung_zero, group.by = "SCT_snn_res.0.3", label =TRUE)
z5 <-DimPlot(lung_zero, group.by = "SCT_snn_res.0.5", label =TRUE)
z6 <-DimPlot(lung_zero, group.by = "SCT_snn_res.0.6", label =TRUE)
z7 <-DimPlot(lung_zero, group.by = "SCT_snn_res.0.7", label =TRUE)
z8 <-DimPlot(lung_zero, group.by = "SCT_snn_res.0.8", label =TRUE)
z1.1 <-DimPlot(lung_zero, group.by = "SCT_snn_res.1", label =TRUE)
z1|z2|z3
z5|z6|z7
z8|z1.1
#decided to go with 0.8 resolution
lung_zero <- FindClusters(lung_zero, resolution = 0.2, graph.name = "SCT_snn")

lung_zeroDP<-DimPlot(lung_zero, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
lung_zeroDP

#check for batch effects
#looks ok, so will not run harmony again
DimPlot(lung_zero, reduction = 'umap', group.by = 'orig.ident', label = FALSE)

#now identify lung_ zero clusters
DefaultAssay(lung_zero) <- "RNA"
DefaultAssay(lung_zero)

#now run cluster maker identification
Zmarkers_cluster0 <- FindConservedMarkers(lung_zero,ident.1 = 0,grouping.var = 'orig.ident')
head(Zmarkers_cluster0)
write.table(Zmarkers_cluster0, file = 'Zmarkers_cluster0.csv', sep =',')

Zmarkers_cluster1 <- FindConservedMarkers(lung_zero,ident.1 = 1,grouping.var = 'orig.ident')
head(Zmarkers_cluster1)
write.table(Zmarkers_cluster1, file = 'Zmarkers_cluster1.csv', sep =',')

Zmarkers_cluster2 <- FindConservedMarkers(lung_zero,ident.1 = 2,grouping.var = 'orig.ident')
head(Zmarkers_cluster2)
write.table(Zmarkers_cluster2, file = 'Zmarkers_cluster2.csv', sep =',')

Zmarkers_cluster3 <- FindConservedMarkers(lung_zero,ident.1 = 3,grouping.var = 'orig.ident')
head(Zmarkers_cluster3)
write.table(Zmarkers_cluster3, file = 'Zmarkers_cluster3.csv', sep =',')

Zmarkers_cluster4 <- FindConservedMarkers(lung_zero,ident.1 = 4,grouping.var = 'orig.ident')
head(Zmarkers_cluster4)
write.table(Zmarkers_cluster4, file = 'Zmarkers_cluster4.csv', sep =',')

# Checking all makers 
DefaultAssay(lung_zero) <- "RNA"
DefaultAssay(lung_zero)

All_Zmarkers <- FindAllMarkers(object = lung_zero, only.pos = TRUE,logfc.threshold = 0.1, min.pct = )
All_Zmarkers  <- All_Zmarkers [ , c(6, 7, 2:4, 1, 5)]
All_Zmarkers<- All_Zmarkers %>%
  dplyr::arrange(cluster, p_val_adj)

View(All_Zmarkers)
write.table(All_Zmarkers, file = 'All_Zmarkers.csv', sep =',')

Idents(lung_zero)<- 'seurat_clusters'
table(Idents(lung_zero))

#Renaming Idents to identify cell clusters
lung_zero <- RenameIdents(lung_zero, "0" = "T Cells")
lung_zero <- RenameIdents(lung_zero, "1" = "A Cap")
lung_zero <- RenameIdents(lung_zero, "2" = "G Cap")
lung_zero <- RenameIdents(lung_zero, "3" = "B Cells")
lung_zero <- RenameIdents(lung_zero, "4" = "Monocytes")

# visualize newly labelled Dim Plit
lung_zero_LabelledDimP <- DimPlot(lung_zero, label = TRUE)
lung_zero_LabelledDimP

# can now stash these identities as a new column in the meta data
lung_zero <- StashIdent(lung_zero, save.name = "subcluster.label")
#then re-assign back to lung object
lung$subcluster.label <- as.character(Idents(lung))
lung$subcluster.label[Cells(lung_zero)] <- paste("Subcluster0_",Idents(lung_zero))
View(lung@meta.data)
Idents(lung) <- "subcluster.label"
lung_LabelledDimP <- DimPlot(lung, label = TRUE)
lung_LabelledDimP

# save file with new idents & save subcluster analysis object
saveRDS(lung, file = "C:/Users/Kevin/Dropbox/KAM R Code/SCT_12122022/lung.rds")
saveRDS(lung_zero, file = "C:/Users/Kevin/Dropbox/KAM R Code/SCT_12122022/lung_zero.rds")

##CLUSTER TWO
#need to subset cluster 2 to get cell types as well
Idents(lung)<- 'seurat_clusters'
this_cluster <-'2'
lung_two <- subset(lung, idents = this_cluster)

Idents(lung_two)<- 'orig.ident'
levels(lung_two)
DefaultAssay(lung_two) <- "RNA"

split_lung_two <- SplitObject(lung_two, split.by = "orig.ident")
split_lung_two <- lapply(X = split_lung_two, 
                         FUN = SCTransform, 
                         method = "glmGamPoi", 
                         return.only.var.genes = FALSE,
                         vars.to.regress = c("percent.mito","S.Score","G2M.Score"))


var.features <- SelectIntegrationFeatures(object.list = split_lung_two, nfeatures = 3000)                     
lung_two <- merge(x = split_lung_two[[1]], y = split_lung_two[2:length(split_lung_two)], merge.data=TRUE)
VariableFeatures(lung_two) <- var.features
lung_two <- RunPCA(lung_two, verbose = FALSE)
# will check to see if need harmony here as doesn't converge if run on this subset alone
# think already data is too similar to adjust again since it is a subcluster 
lung_two <- RunUMAP(lung_two, reduction = "pca", dims = 1:30)
lung_two <- FindNeighbors(lung_two, reduction = "pca", dims = 1:30) 
lung_two <- FindClusters(lung_two, resolution = c(0.1, 0.2, 0.3, 0.5, 0.6, 0.7,0.8, 1))
t1 <-DimPlot(lung_two, group.by = "SCT_snn_res.0.1", label =TRUE)
t2 <-DimPlot(lung_two, group.by = "SCT_snn_res.0.2", label =TRUE)
t3 <-DimPlot(lung_two, group.by = "SCT_snn_res.0.3", label =TRUE)
t5 <-DimPlot(lung_two, group.by = "SCT_snn_res.0.5", label =TRUE)
t6 <-DimPlot(lung_two, group.by = "SCT_snn_res.0.6", label =TRUE)
t7 <-DimPlot(lung_two, group.by = "SCT_snn_res.0.7", label =TRUE)
t8 <-DimPlot(lung_two, group.by = "SCT_snn_res.0.8", label =TRUE)
t1.1 <-DimPlot(lung_two, group.by = "SCT_snn_res.1", label =TRUE)
t1|t2|t3
t5|t6|t7
t8|t1.1
#decided to go with 0.3 resolution
lung_two <- FindClusters(lung_two, resolution = 0.3, graph.name = "SCT_snn")

lung_twoDP<-DimPlot(lung_two, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
lung_twoDP

#check for batch effects
#looks ok, so will not run harmony again
DimPlot(lung_two, reduction = 'umap', group.by = 'orig.ident', label = FALSE)

#now identify lung_two clusters
DefaultAssay(lung_two) <- "RNA"
DefaultAssay(lung_two)

#now run cluster maker identification
Tmarkers_cluster0 <- FindConservedMarkers(lung_two,ident.1 = 0,grouping.var = 'orig.ident')
head(Tmarkers_cluster0)
write.table(Tmarkers_cluster0, file = 'Tmarkers_cluster0.csv', sep =',')

Tmarkers_cluster1 <- FindConservedMarkers(lung_two,ident.1 = 1,grouping.var = 'orig.ident')
head(Tmarkers_cluster1)
write.table(Tmarkers_cluster1, file = 'Tmarkers_cluster1.csv', sep =',')

Tmarkers_cluster2 <- FindConservedMarkers(lung_two,ident.1 = 2,grouping.var = 'orig.ident')
head(Tmarkers_cluster2)
write.table(Tmarkers_cluster2, file = 'Tmarkers_cluster2.csv', sep =',')

Tmarkers_cluster3 <- FindConservedMarkers(lung_two,ident.1 = 3,grouping.var = 'orig.ident')
head(Tmarkers_cluster3)
write.table(Tmarkers_cluster3, file = 'Tmarkers_cluster3.csv', sep =',')

Tmarkers_cluster4 <- FindConservedMarkers(lung_two,ident.1 = 4,grouping.var = 'orig.ident')
head(Tmarkers_cluster4)
write.table(Tmarkers_cluster4, file = 'Tmarkers_cluster4.csv', sep =',')

# Checking all makers 
DefaultAssay(lung_two) <- "RNA"
DefaultAssay(lung_two)

All_Tmarkers <- FindAllMarkers(object = lung_two, only.pos = TRUE,logfc.threshold = 0.1)
All_Tmarkers  <- All_Tmarkers [ , c(6, 7, 2:4, 1, 5)]
All_Tmarkers<- All_Tmarkers %>%
  dplyr::arrange(cluster, p_val_adj)

View(All_Tmarkers)
write.table(All_Tmarkers, file = 'All_Tmarkers.csv', sep =',')

Idents(lung_two)<- 'seurat_clusters'
table(Idents(lung_two))

#Renaming Idents to identify cell clusters
lung_two <- RenameIdents(lung_two, "0" = "T cells")
lung_two <- RenameIdents(lung_two, "1" = "Neutrophils")
lung_two <- RenameIdents(lung_two, "2" = "Monocytes")
lung_two <- RenameIdents(lung_two, "3" = "A Cap")
lung_two <- RenameIdents(lung_two, "4" = "B Cells")


# visualize newly labelled Dim Plot
lung_two_LabelledDimP <- DimPlot(lung_two, label = TRUE)
lung_two_LabelledDimP

# can now stash these identities as a new column in the meta data
lung_two <- StashIdent(lung_two, save.name = "subcluster.label")
#then re-assign back to lung object
Idents(lung)<- 'subcluster.label'
lung$subcluster.label <- as.character(Idents(lung))
lung$subcluster.label[Cells(lung_two)] <- paste("Subcluster2_",Idents(lung_two))
View(lung@meta.data)
Idents(lung) <- "subcluster.label"
lung_LabelledDimP <- DimPlot(lung, label = TRUE)
lung_LabelledDimP

# save file with new idents & save subcluster analysis object
saveRDS(lung, file = "C:/Users/Kevin/Dropbox/KAM R Code/SCT_12122022/lung.rds")
saveRDS(lung_two, file = "C:/Users/Kevin/Dropbox/KAM R Code/SCT_12122022/lung_two.rds")

# using RenameIdents to recode, and then will stash idents into main data so don't lose original seurat clusters
Idents(lung) <- "subcluster.label"
lung <- RenameIdents(lung, "1" = "B Cells")
lung <- RenameIdents(lung, "3" = "Neutrophils.2")
lung <- RenameIdents(lung, "4" = "Capillary Aerocytes")
lung <- RenameIdents(lung, "5" = "Col13a1+ Fibroblasts")
lung <- RenameIdents(lung, "6" = "Col14a1+ Fibroblasts")
lung <- RenameIdents(lung, "7" = "Alveolar Macrophages")
lung <- RenameIdents(lung, "8" = "Neutrophils.1")
lung <- RenameIdents(lung, "9" = "NK Cells")
lung <- RenameIdents(lung, "10" = "Club Cells")
lung <- RenameIdents(lung, "11" = "T Cells")
lung <- RenameIdents(lung, "12" = "Pericytes")
lung <- RenameIdents(lung, "13" = "Interstitial Macrophages")
lung <- RenameIdents(lung, "14" = "Mesothelial")
lung <- RenameIdents(lung, "15" = "Platelets")
lung <- RenameIdents(lung, "16" = "Vein Endothelia")
lung <- RenameIdents(lung, "17" = "Myofibroblasts")
lung <- RenameIdents(lung, "18" = "Fibroblasts")
lung <- RenameIdents(lung, "Subcluster0_ A Cap" = "Capillary Aerocytes")
lung <- RenameIdents(lung, "Subcluster0_ G Cap" = "General Capillary")
lung <- RenameIdents(lung, "Subcluster2_ T cells" = "T Cells.1")
lung <- RenameIdents(lung, "Subcluster2_ A Cap" = "A Cap")
lung <- RenameIdents(lung, "Subcluster0_ T Cells" = "T Cells.2")
lung <- RenameIdents(lung, "Subcluster2_ Monocytes" = "Monocytes.1")
lung <- RenameIdents(lung, "Subcluster2_ B Cells" = "B Cells.1")
lung <- RenameIdents(lung, "Subcluster2_ Neutrophils" = "Neutrophils.3")
lung <- RenameIdents(lung, "Subcluster0_ Monocytes" = "Monocytes.2")
lung <- RenameIdents(lung, "Subcluster0_ B Cells" = "B Cells.2")

#Checked some clus ters with FindMarkers with both comparions and all identifiers and relabelled some clusters (see example)
markers <- FindMarkers(lung, ident.1 = "T Cells", ident.2 = c("T Cells.1"), min.pct = 0.25)
head(markers, n = 10)

lung <- RenameIdents(lung, "T Cells.1" = "Monocytes")
lung <- RenameIdents(lung, "Monocytes.1" = "Monocytes")
lung <- RenameIdents(lung, "Monocytes.2" = "Monocytes")
lung <- RenameIdents(lung, "T Cells.2" = "General Capillary")
lung <- RenameIdents(lung, "B Cells.2" = "B Cells")
lung <- RenameIdents(lung, "B Cells.1" = "B Cells")
lung <- RenameIdents(lung, "A Cap" = "Mixed Endothelia")
lung <- RenameIdents(lung, "Neutrophils.1" = "Type 1 Neutrophils")
lung <- RenameIdents(lung, "Neutrophils.2" = "Type 2 Neutrophils")
lung <- RenameIdents(lung, "Neutrophils.3" = "T Cells")

# visualize newly labelled Dim Plit
lung_DP <- DimPlot(lung, label = TRUE)
lung_DP

# can now stash these identities as a new column in the meta data
lung <- StashIdent(lung, save.name = "main.cluster")

saveRDS(lung, file = "C:/Users/Katherine/Dropbox/KAM R Code/SCT_12122022/lung.rds")

# setting the idents for next analysis and double checking levels are correct
Idents(lung) <- "main.cluster"
levels(lung)
lung$main.cluster <- factor(lung$main.cluster, levels = c("T Cells" , "Type 2 Neutrophils", "Type 1 Neutrophils", "B Cells","Monocytes",  "NK Cells", "Alveolar Macrophages",
                                                          "Interstitial Macrophages",  "Mixed Endothelia","General Capillary","Capillary Aerocytes", "Vein Endothelia","Mesothelial",
                                                          "Pericytes", "Club Cells","Col14a1+ Fibroblasts", "Col13a1+ Fibroblasts","Fibroblasts" , "Myofibroblasts", "Platelets"))
# checking now it looks split by timepoint and getting basic count/percentage data
split_DP <- DimPlot(lung, reduction = "umap", split.by = "orig.ident", ncol = 2)
split_DP

# plot clusters as proportion or percentage of total cell types at each timepoint
lung$orig.ident <- factor(lung$orig.ident, levels = c("N0", "S12", "T12", "S3", "T3", "S28", "T28"))
Idents(lung) <- "orig.ident"
all_groups_FreP <- ggplot(lung@meta.data, aes(x=orig.ident, fill=main.cluster)) + geom_bar(position = "fill")
all_groups_FreP

# plot as just counts
all_groups_CountP <-ggplot(lung@meta.data, aes(x=orig.ident, fill=main.cluster)) + geom_bar()
all_groups_CountP

# plot of all cells in each cluster
all_cells_FreP <- ggplot(lung@meta.data, aes(x=dub.idents, fill=main.cluster)) + geom_bar(position = "fill")
all_cells_FreP

# plot clusters as proportion over time
all_PropXTime <-ggplot(lung@meta.data, aes(x=main.cluster, fill=orig.ident)) + geom_bar()
all_PropXTime

# cell count table & how to export it to csv file
Idents(lung) <- "main.cluster"
cluster_count<- table(Idents(lung), lung$orig.ident)    
cluster_count
write.table(cluster_count, file = 'cluster_count.csv', sep =',')

#save object with new main idents stashed
saveRDS(lung, file = "C:/Users/Kevin/Dropbox/KAM R Code/SCT_12122022/lung.rds")


#Creating plots for final makers for clusters
final.cluster.markers <- c('Plac8','Pou2f2','Fcer1g','Ifi27l2a','Cybb',
                           'Tcf7','Cd8a','Cd4',"Cd8b1","Ms4a6b",
                           'Retnlg','S100a9',	'S100a8',	'Cxcr2',	'Ifitm1',
                           'Stfa3','Lcn2','Ltf','Ngp','Camp','Chil3',
                           'Ms4a1',	'Cd74',	'Igkc',	'Ighm',	'Cd79a',
                           'Calcrl',	'Pecam1',	'Epas1',	'Cd93',	'Ptprb',
                           'Ednrb',	'Car4',	'Emp2',	'Kdr',	'Kitl',
                           'Hhip',	'Enpp2',	'Aspn',	'Igfbp3',	'Mgp',
                           'Mfap5',	
                           'Vegfc',	'Amigo2',	'Prss23',	'Slc6a2',	'Vwf',
                           'Nrgn',	'Itga2b',	'Ppbp',	'Pf4',	'Tsc22d1',
                           'Msln',	'Rarres2',	'Igfbp5',	'Igfbp6',	'Upk3b',
                           'C1qa',	'Apoe',	'C1qb',	'C1qc',	'Ctsb',
                           'Gucy1b1',	'Postn',	'Gucy1a1',	'Pdzd2',	'Pde5a',
                           'Scgb1a1',	'Scgb3a1',	'Scgb3a2',	'Bpifa1',	'Cyp2f2',
                           'Gzma',	'Nkg7',	'Ccl5',	'Id2',	'AW112010',
                           'Ear2',	'Chil3',	'Ctsd',	'Ear1',	'Lpl',
                           'Col1a2',	'Col1a1',	'Col14a1',	'Dcn',	'Gsn',
                           'Inmt',	'Fmo2',	'Gpx3',	'Limch1',	'Aldh1a1',
                           'Pxdn',	'Npr3',	'Notch4',	'Col4a1',	'Bcam')

# to change order of groups for plotting
lung$main.cluster <- factor(lung$main.cluster, levels = c("Monocytes",
                                                          "T Cells",
                                                          "Type 1 Neutrophils",
                                                          "Type 2 Neutrophils",
                                                          "B Cells",
                                                          "General Capillary",
                                                          "Capillary Aerocytes",
                                                          "Myofibroblasts",
                                                          "Fibroblasts",
                                                          "Vein Endothelia",
                                                          "Basophils",
                                                          "Mesothelial",
                                                          "Interstitial Macrophages",
                                                          "Pericytes",
                                                          "Club Cells",
                                                          "NK Cells",
                                                          "Alveolar Macrophages",
                                                          "Col14a1+ Fibroblasts",
                                                          "Col13a1+ Fibroblasts",
                                                          "Mixed Endothelia"))
Idents(lung) <- "main.cluster"
levels(lung)

# marker dot plot
lung_DotP <-DotPlot(lung, features = rev(final.cluster.markers), assay = "RNA", cols = c("blue", "orange"), dot.scale = 6) + RotatedAxis()
lung_DotP

#to create heat map of cluster markers average expression
all_lung <- AverageExpression(lung, features = (final.cluster.markers), return.seurat = T)
lung_HeatM <-DoHeatmap(all_lung, features = rev(final.cluster.markers), draw.lines = FALSE) + scale_fill_gradientn(colors = c("black", "midnightblue", "magenta4","orange", "yellow"))
lung_HeatM

#table of total cell numbers
table(Idents(lung))







