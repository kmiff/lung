#set working directory
setwd("C:/Users/kathe/Dropbox/KAM R Code/Augur")

#open necessary libraries
library(Seurat)
library(Augur)
library(sparseMatrixStats)
library(tidyverse)
library(viridis)
library(magrittr)

#get seurat object 
lung <- readRDS("C:/Users/Katherine/Dropbox/KAM R Code/SCT_12122022/lung.rds")
lung
#fist get raw count and metadata files for analysis
DefaultAssay(lung) <- 'RNA'
Idents(lung)<- 'main.cluster'
#count data
counts_lung <- GetAssayData(lung, assay = "RNA", slot = "counts") # raw counts
write.table(counts_lung, file = 'counts_lung.csv', sep =',')
#meta data
meta_lung <-lung@meta.data
write.table(meta_lung, file = 'meta_lung.csv', sep =',')

#run Augur on default settings using counts & metadata files
#use multiclass augur function as multiple timepoints in analysis
multiclass_augur = calculate_auc(counts_lung, meta = meta_lung, cell_type_col = "main.cluster", label_col = "orig.ident")
print(multiclass_augur$AUC)
plot_umap(multiclass_augur, lung, cell_type_col = "main.cluster")
plot_lollipop(multiclass_augur)
saveRDS(multiclass_augur, file = "C:/Users/Kathe/Dropbox/KAM R Code/Augur/multiclass_augur_raw.rds")
