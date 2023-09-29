library(Seurat)
library(patchwork)
library(ggplot2)
library(igraph)

setwd("C:/Users/Katherine/Dropbox/KAM R Code/SCT_12122022/CellChat")
#get seurat object 
lung <- readRDS("C:/Users/Katherine/Dropbox/KAM R Code/SCT_12122022/lung.rds")
lung
#set idents
Idents(lung) <- "main.cluster"
levels(lung)
DefaultAssay(lung) <- "RNA"
#find DEGs
#type 1 neutrophils
#Naive vs 3d SCI (T3)
Neu1_N0vT3 <- FindMarkers(object = lung, 
                          ident.1 = WhichCells(
                            lung, 
                            cells = rownames(lung@meta.data[lung@meta.data$orig.ident=="T3",]), 
                            idents = "Type 1 Neutrophils"),
                          ident.2 = WhichCells(lung, cells = rownames(lung@meta.data[lung@meta.data$orig.ident=="N0",]), idents = "Type 1 Neutrophils"), 
                          min.pct = 0.1, logfc.threshold = 0.176, assay = "RNA", only.pos = FALSE)
write.csv(Neu1_N0vT3, "Neu1_N0vT3.csv")
#Naive (N0)  vs 12h SCI (T12)
Neu1_N0vT12 <- FindMarkers(object = lung, 
                           ident.1 = WhichCells(
                             lung, 
                             cells = rownames(lung@meta.data[lung@meta.data$orig.ident=="T12",]), 
                             idents = "Type 1 Neutrophils"),
                           ident.2 = WhichCells(lung, cells = rownames(lung@meta.data[lung@meta.data$orig.ident=="N0",]), idents = "Type 1 Neutrophils"), 
                           min.pct = 0.1, logfc.threshold = 0.176, assay = "RNA", only.pos = FALSE)
write.csv(Neu1_N0vT12, "Neu1_N0vT12.csv")
#Naive (N0) vs 28d SCI (T28)
Neu1_N0vT28 <- FindMarkers(object = lung, 
                           ident.1 = WhichCells(
                             lung, 
                             cells = rownames(lung@meta.data[lung@meta.data$orig.ident=="T28",]), 
                             idents = "Type 1 Neutrophils"),
                           ident.2 = WhichCells(lung, cells = rownames(lung@meta.data[lung@meta.data$orig.ident=="N0",]), idents = "Type 1 Neutrophils"), 
                           min.pct = 0.1, logfc.threshold = 0.176, assay = "RNA", only.pos = FALSE)
write.csv(Neu1_N0vT28, "Neu1_N0vT28.csv")
#Naive (N0) vs 3d Sham (S3)
Neu1_N0vS3 <- FindMarkers(object = lung, 
                           ident.1 = WhichCells(
                             lung, 
                             cells = rownames(lung@meta.data[lung@meta.data$orig.ident=="S3",]), 
                             idents = "Type 1 Neutrophils"),
                           ident.2 = WhichCells(lung, cells = rownames(lung@meta.data[lung@meta.data$orig.ident=="N0",]), idents = "Type 1 Neutrophils"), 
                           min.pct = 0.1, logfc.threshold = 0.176, assay = "RNA", only.pos = FALSE)
write.csv(Neu1_N0vS3, "Neu1_N0vS3.csv")
#Naive (N0) vs 28d Sham (S28)
Neu1_N0vS28 <- FindMarkers(object = lung, 
                           ident.1 = WhichCells(
                             lung, 
                             cells = rownames(lung@meta.data[lung@meta.data$orig.ident=="S28",]), 
                             idents = "Type 1 Neutrophils"),
                           ident.2 = WhichCells(lung, cells = rownames(lung@meta.data[lung@meta.data$orig.ident=="N0",]), idents = "Type 1 Neutrophils"), 
                           min.pct = 0.1, logfc.threshold = 0.176, assay = "RNA", only.pos = FALSE)
write.csv(Neu1_N0vS28, "Neu1_N0vS28.csv")
#Naive (N0) vs 12h Sham (S12)
Neu1_N0vS12 <- FindMarkers(object = lung, 
                           ident.1 = WhichCells(
                             lung, 
                             cells = rownames(lung@meta.data[lung@meta.data$orig.ident=="S12",]), 
                             idents = "Type 1 Neutrophils"),
                           ident.2 = WhichCells(lung, cells = rownames(lung@meta.data[lung@meta.data$orig.ident=="N0",]), idents = "Type 1 Neutrophils"), 
                           min.pct = 0.1, logfc.threshold = 0.176, assay = "RNA", only.pos = FALSE)
write.csv(Neu1_N0vS12, "Neu1_N0vS12.csv")
#Note
#DEGs were determined to be from these genes lists after 
#converting the p_val_adj to the negative log 10 p value adjusted
#any neg log10 p val adj greater than 2 was considered a sig DEG
#plotted as a volcano plot in prism
#these gene lists of increased (those with a positive ave_log2FC) 
#and decreased (those with a negative avg_log2FC) genes were the put into
#gProfiler (https://biit.cs.ut.ee/gprofiler/gost) for pathway analysis
#parameters used: 
#organism: mus musculus (Mouse)
#only annotated genes
# signficance threshold: benjamini-hochberg FDR at 0.05
#data sources: GO Biological process and KEGG
#top 10 pathways with highest -log1-(p.edj) were graphed in prims and reported
#to get overlapping SCI and SHAM degs, the lists of genes as determined above were put into
#https://bioinformatics.psb.ugent.be/webtools/Venn/ to create venn diagrams
# the lists of overlapping and unique genes were then put into gProfiler as above