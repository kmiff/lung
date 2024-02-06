#0.0 REFERENCES-----

#1.0 CELL CHAT OBJECT PREP----
##1.1 Load Libraries & Set Directories----
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(igraph)
library(Hmisc)
library(hash)
library(tidyverse)
library(circlize)
library(ggsci)
library(gtools)
library(ComplexHeatmap)
#set directory
setwd("C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023")

##1.2 Load Seurat Object---- 
lung <- readRDS("C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/obj2_12032023.rds")
lung

##1.3 Update Cell Chat Database----
##has missing genes, so need to update
# see tutorial: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Update-CellChatDB.html
CellChatDB <- CellChatDB.mouse 
interaction_input <- CellChatDB$interaction
complex_input <- CellChatDB$complex
cofactor_input <- CellChatDB$cofactor
geneInfo <- CellChatDB$geneInfo
write.csv(interaction_input, file = "interaction_input_CellChatDB.csv")
write.csv(complex_input, file = "complex_input_CellChatDB.csv")
write.csv(cofactor_input, file = "cofactor_input_CellChatDB.csv")
write.csv(geneInfo, file = "geneInfo_input_CellChatDB.csv")

#updated H2-Ea-ps to gene name H2-Ea in interaction file 
#updated  H2-BI to H2-Bl in interaction & geneInfo file based off of KEGG:mmu04612 pathway interaction with listed ligand KIR3DL1
#now upload updated files
options(stringsAsFactors = FALSE)
interaction_input <- read.csv(file = 'interaction_input_CellChatDB.csv', row.names = 1)
complex_input <- read.csv(file = 'complex_input_CellChatDB.csv', row.names = 1)
cofactor_input <- read.csv(file = 'cofactor_input_CellChatDB.csv', row.names = 1)
geneInfo <- read.csv(file = 'geneInfo_input_CellChatDB.csv', row.names = 1)
CellChatDB <- list()
CellChatDB$interaction <- interaction_input
CellChatDB$complex <- complex_input
CellChatDB$cofactor <- cofactor_input
CellChatDB$geneInfo <- geneInfo
#overwrite masterfile
setwd ("C:/Users/Kathe/AppData/Local/R/win-library/4.2/CellChat") # This is the folder of CellChat package downloaded from Github
CellChatDB.mouse <- CellChatDB
usethis::use_data(CellChatDB.mouse, overwrite = TRUE)

##1.4 Split Object----
split_lung <- SplitObject(lung, split.by = "orig.ident")

##1.5 Convert to Cell Chat Objects----
###1.5a Cellchat_N0----
N0 <-split_lung$N0
DefaultAssay(N0) <- 'RNA'
Idents(N0)<- 'main.cluster'
levels(N0)
data.input_N0 <- GetAssayData(N0, assay = "RNA", slot = "data") # normalized data matrix, data is the log noramlized version of counts, counts is the raw data
labels_N0 <- Idents(N0)
meta_N0 <- data.frame(group = labels_N0, row.names = names(labels_N0)) # create a dataframe of the cell labels
cellchat_N0 <- createCellChat(object = data.input_N0, meta = meta_N0, group.by = 'group')
CellChatDB <- CellChatDB.mouse
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat_N0@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat_N0 <- subsetData(cellchat_N0) # This step is necessary even if using the whole database
cellchat_N0 <- identifyOverExpressedGenes(cellchat_N0)
cellchat_N0 <- identifyOverExpressedInteractions(cellchat_N0)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat_N0 <- projectData(cellchat_N0, PPI.mouse)
cellchat_N0 <- computeCommunProb(cellchat_N0,type = "truncatedMean", trim = 0.1, raw.use = T, population.size = T)
#note:  Set population.size = TRUE if analyzing unsorted single-cell transcriptomes, 
# with the reason that abundant cell populations tend to send collectively stronger signals 
#than the rare cell populations.
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_N0 <- filterCommunication(cellchat_N0, min.cells = 10)
#return a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
slot.name = "netP"
df.net <- subsetCommunication(cellchat_N0)
#Infer the cell-cell communication at a signaling pathway level
cellchat_N0 <- computeCommunProbPathway(cellchat_N0)
#Calculate the aggregated cell-cell communication network
cellchat_N0 <- aggregateNet(cellchat_N0)
# Compute the network centrality scores
cellchat_N0 <- netAnalysis_computeCentrality(cellchat_N0, slot.name = "netP")
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
saveRDS(cellchat_N0, file ="C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023/cellchat_N0_12072023.rds")

###1.5b Cellchat_S12----
S12 <-split_lung$S12
DefaultAssay(S12) <- 'RNA'
Idents(S12)<- 'main.cluster'
levels(S12)
data.input_S12 <- GetAssayData(S12, assay = "RNA", slot = "data") # normalized data matrix
labels_S12 <- Idents(S12)
meta_S12 <- data.frame(group = labels_S12, row.names = names(labels_S12)) # create a dataframe of the cell labels
cellchat_S12 <- createCellChat(object = data.input_S12, meta = meta_S12, group.by = 'group')
# set the used database in the object
cellchat_S12@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat_S12 <- subsetData(cellchat_S12) # This step is necessary even if using the whole database
cellchat_S12 <- identifyOverExpressedGenes(cellchat_S12)
cellchat_S12 <- identifyOverExpressedInteractions(cellchat_S12)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat_S12 <- projectData(cellchat_S12, PPI.mouse)
cellchat_S12 <- computeCommunProb(cellchat_S12,type = "truncatedMean", trim = 0.1, raw.use = T, population.size = T)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_S12 <- filterCommunication(cellchat_S12, min.cells = 10)
#return a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
slot.name = "netP"
df.net <- subsetCommunication(cellchat_S12)
#Infer the cell-cell communication at a signaling pathway level
cellchat_S12 <- computeCommunProbPathway(cellchat_S12)
#Calculate the aggregated cell-cell communication network
cellchat_S12 <- aggregateNet(cellchat_S12)
# Compute the network centrality scores
cellchat_S12 <- netAnalysis_computeCentrality(cellchat_S12, slot.name = "netP")
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
saveRDS(cellchat_S12, file = "C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023/cellchat_S12_12072023.rds")

###1.5c Cellchat_T12----
T12 <-split_lung$T12
DefaultAssay(T12) <- 'RNA'
Idents(T12)<- 'main.cluster'
levels(T12)
data.input_T12 <- GetAssayData(T12, assay = "RNA", slot = "data") # normalized data matrix
labels_T12  <- Idents(T12)
meta_T12  <- data.frame(group = labels_T12 , row.names = names(labels_T12 )) # create a dataframe of the cell labels
cellchat_T12 <- createCellChat(object = data.input_T12 , meta = meta_T12 , group.by = 'group')
# set the used database in the object
cellchat_T12@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat_T12 <- subsetData(cellchat_T12) # This step is necessary even if using the whole database
cellchat_T12 <- identifyOverExpressedGenes(cellchat_T12)
cellchat_T12 <- identifyOverExpressedInteractions(cellchat_T12)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat_T12 <- projectData(cellchat_T12, PPI.mouse)
cellchat_T12<- computeCommunProb(cellchat_T12,type = "truncatedMean", trim = 0.1, raw.use = T, population.size = T)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_T12 <- filterCommunication(cellchat_T12, min.cells = 10)
#return a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
slot.name = "netP"
df.net <- subsetCommunication(cellchat_T12)
#Infer the cell-cell communication at a signaling pathway level
cellchat_T12 <- computeCommunProbPathway(cellchat_T12)
#Calculate the aggregated cell-cell communication network
cellchat_T12 <- aggregateNet(cellchat_T12)
# Compute the network centrality scores
cellchat_T12 <- netAnalysis_computeCentrality(cellchat_T12, slot.name = "netP")
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
saveRDS(cellchat_T12, file = "C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023/cellchat_T12_12072023.rds")

###1.5d Cellchat_S3----
S3 <-split_lung$S3
DefaultAssay(S3) <- 'RNA'
Idents(S3)<- 'main.cluster'
levels(S3)
data.input_S3 <- GetAssayData(S3, assay = "SCT", slot = "data") # normalized data matrix
labels_S3 <- Idents(S3)
meta_S3 <- data.frame(group = labels_S3, row.names = names(labels_S3)) # create a dataframe of the cell labels
cellchat_S3 <- createCellChat(object = data.input_S3, meta = meta_S3, group.by = 'group')
# set the used database in the object
cellchat_S3@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat_S3 <- subsetData(cellchat_S3) # This step is necessary even if using the whole database
cellchat_S3 <- identifyOverExpressedGenes(cellchat_S3)
cellchat_S3 <- identifyOverExpressedInteractions(cellchat_S3)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat_S3 <- projectData(cellchat_S3, PPI.mouse)
cellchat_S3<- computeCommunProb(cellchat_S3,type = "truncatedMean", trim = 0.1, raw.use = T, population.size = T)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_S3 <- filterCommunication(cellchat_S3, min.cells = 10)
#return a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
slot.name = "netP"
df.net <- subsetCommunication(cellchat_S3)
#Infer the cell-cell communication at a signaling pathway level
cellchat_S3 <- computeCommunProbPathway(cellchat_S3)
#Calculate the aggregated cell-cell communication network
cellchat_S3 <- aggregateNet(cellchat_S3)
# Compute the network centrality scores
cellchat_S3 <- netAnalysis_computeCentrality(cellchat_S3, slot.name = "netP")
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
saveRDS(cellchat_S3, file = "C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023/cellchat_S3_12072023.rds")

###1.5e Cellchat_T3----
T3 <-split_lung$T3
DefaultAssay(T3) <- 'RNA'
Idents(T3)<- 'main.cluster'
levels(T3)
data.input_T3 <- GetAssayData(T3, assay = "RNA", slot = "data") # normalized data matrix
labels_T3 <- Idents(T3)
meta_T3 <- data.frame(group = labels_T3, row.names = names(labels_T3)) # create a dataframe of the cell labels
cellchat_T3 <- createCellChat(object = data.input_T3, meta = meta_T3, group.by = 'group')
# set the used database in the object
cellchat_T3@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat_T3 <- subsetData(cellchat_T3) # This step is necessary even if using the whole database
cellchat_T3 <- identifyOverExpressedGenes(cellchat_T3)
cellchat_T3 <- identifyOverExpressedInteractions(cellchat_T3)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat_T3<- projectData(cellchat_T3, PPI.mouse)
cellchat_T3<- computeCommunProb(cellchat_T3,type = "truncatedMean", trim = 0.1, raw.use = T, population.size = T)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_T3 <- filterCommunication(cellchat_T3, min.cells = 10)
#return a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
slot.name = "netP"
df.net <- subsetCommunication(cellchat_T3)
#Infer the cell-cell communication at a signaling pathway level
cellchat_T3 <- computeCommunProbPathway(cellchat_T3)
#Calculate the aggregated cell-cell communication network
cellchat_T3 <- aggregateNet(cellchat_T3)
# Compute the network centrality scores
cellchat_T3 <- netAnalysis_computeCentrality(cellchat_T3, slot.name = "netP")
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
saveRDS(cellchat_T3, file = "C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023/cellchat_T3_12072023.rds")

###1.5f Cellchat_S28----
S28 <-split_lung$S28
DefaultAssay(S28) <- 'RNA'
Idents(S28)<- 'main.cluster'
levels(S28)
data.input_S28 <- GetAssayData(S28, assay = "RNA", slot = "data") # normalized data matrix
labels_S28 <- Idents(S28)
meta_S28 <- data.frame(group = labels_S28, row.names = names(labels_S28)) # create a dataframe of the cell labels
cellchat_S28 <- createCellChat(object = data.input_S28, meta = meta_S28, group.by = 'group')
# set the used database in the object
cellchat_S28@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat_S28 <- subsetData(cellchat_S28) # This step is necessary even if using the whole database
cellchat_S28 <- identifyOverExpressedGenes(cellchat_S28)
cellchat_S28 <- identifyOverExpressedInteractions(cellchat_S28)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat_S28<- projectData(cellchat_S28, PPI.mouse)
cellchat_S28<- computeCommunProb(cellchat_S28,type = "truncatedMean", trim = 0.1, raw.use = T, population.size = T)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_S28 <- filterCommunication(cellchat_S28, min.cells = 10)
#return a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
slot.name = "netP"
df.net <- subsetCommunication(cellchat_S28)
#Infer the cell-cell communication at a signaling pathway level
cellchat_S28 <- computeCommunProbPathway(cellchat_S28)
#Calculate the aggregated cell-cell communication network
cellchat_S28 <- aggregateNet(cellchat_S28)
# Compute the network centrality scores
cellchat_S28 <- netAnalysis_computeCentrality(cellchat_S28, slot.name = "netP")
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
saveRDS(cellchat_S28, file = "C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023/cellchat_S28_12072023.rds")

###1.5g Cellchat_T28----
T28 <-split_lung$T28
DefaultAssay(T28) <- 'RNA'
Idents(T28)<- 'main.cluster'
levels(T28)
data.input_T28 <- GetAssayData(T28, assay = "RNA", slot = "data") # normalized data matrix
labels_T28 <- Idents(T28)
meta_T28 <- data.frame(group = labels_T28, row.names = names(labels_T28)) # create a dataframe of the cell labels
cellchat_T28 <- createCellChat(object = data.input_T28, meta = meta_T28, group.by = 'group')
# set the used database in the object
cellchat_T28@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat_T28<- subsetData(cellchat_T28) # This step is necessary even if using the whole database
cellchat_T28 <- identifyOverExpressedGenes(cellchat_T28)
cellchat_T28 <- identifyOverExpressedInteractions(cellchat_T28)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat_T28<- projectData(cellchat_T28, PPI.mouse)
cellchat_T28<- computeCommunProb(cellchat_T28,type = "truncatedMean", trim = 0.1, raw.use = T, population.size = T)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_T28 <- filterCommunication(cellchat_T28, min.cells = 10)
#return a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
slot.name = "netP"
df.net <- subsetCommunication(cellchat_T28)
#Infer the cell-cell communication at a signaling pathway level
cellchat_T28 <- computeCommunProbPathway(cellchat_T28)
#Calculate the aggregated cell-cell communication network
cellchat_T28 <- aggregateNet(cellchat_T28)
# Compute the network centrality scores
cellchat_T28 <- netAnalysis_computeCentrality(cellchat_T28, slot.name = "netP")
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
saveRDS(cellchat_T28, file = "C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023/cellchat_T28_12072023.rds")

##1.6 Merge Data Sets----
object.list <- list(N0 = cellchat_N0, S12 = cellchat_S12, T12 = cellchat_T12, S3 = cellchat_S3,
                    T3 = cellchat_T3, S28 = cellchat_S28,T28 = cellchat_T28)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
cellchat
saveRDS(cellchat, file = "C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023/cellchat_12072023.rds")

#2.0 PLOTS FOR INCOMING & OUTGOING SIGNALS---- 
#to create plots showing interaction strength for incoming and outgoing signals
# see tutorial: 
# https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)


#3.0 HEATMAPS FOR INTERACTION STRENGTH----
##3.1a N0 vs S12----
list1 <-list(N0 = cellchat_N0, S12 = cellchat_S12)
cellchatN0vS12 <- mergeCellChat(list1, add.names = names(list1), cell.prefix = TRUE)
strengthN0vS12 <- netVisual_heatmap(cellchatN0vS12, measure = "weight")
strengthN0vS12

##3.1b N0 vs T12----
list2 <-list(N0 = cellchat_N0, T12 = cellchat_T12)
cellchatN0vT12 <- mergeCellChat(list2, add.names = names(list2), cell.prefix = TRUE)
strengthN0vT12 <- netVisual_heatmap(cellchatN0vT12, measure = "weight")
strengthN0vT12

##3.1c N0 vs S3----
list3 <-list(N0 = cellchat_N0, S3 = cellchat_S3)
cellchatN0vS3 <- mergeCellChat(list3, add.names = names(list3), cell.prefix = TRUE)
strengthN0vS3 <- netVisual_heatmap(cellchatN0vS3 , measure = "weight")
strengthN0vS3

##3.1d N0 vs T3----
list4 <-list(N0 = cellchat_N0, T3 = cellchat_T3)
cellchatN0vT3 <- mergeCellChat(list4, add.names = names(list4), cell.prefix = TRUE)
strengthN0vT3 <- netVisual_heatmap(cellchatN0vT3, measure = "weight")
strengthN0vT3

##3.1e N0 vs S28----
list5 <-list(N0 = cellchat_N0, S28 = cellchat_S28)
cellchatN0vS28 <- mergeCellChat(list5, add.names = names(list5), cell.prefix = TRUE)
strengthN0vS28 <- netVisual_heatmap(cellchatN0vS28 , measure = "weight")
strengthN0vS28

##3.1f N0 vs T28----
list6 <-list(N0 = cellchat_N0, T28 = cellchat_T28)
cellchatN0vT28 <- mergeCellChat(list6, add.names = names(list6), cell.prefix = TRUE)
strengthN0vT28 <- netVisual_heatmap(cellchatN0vT28, measure = "weight")
strengthN0vT28


#4.0 HEATMAPS FOR NUMBER OF INTERACTIONs ----
##4.1a N0 vs S12----
NumN0vS12<- netVisual_heatmap(cellchatN0vS12)
NumN0vS12

##4.1b N0 vs T12----
NumN0vT12<- netVisual_heatmap(cellchatN0vT12)
NumN0vT12

##4.1c N0 vs S3----
NumN0vS3 <- netVisual_heatmap(cellchatN0vS3)
NumN0vS3

##4.1d N0 vs T3----
NumN0vT3 <- netVisual_heatmap(cellchatN0vT3)
NumN0vT3

##4.1e N0 vs S28----
NumN0vS28<- netVisual_heatmap(cellchatN0vS28)
NumN0vS28

##4.1f N0 vs T28----
NumN0vT28 <- netVisual_heatmap(cellchatN0vT28)
NumN0vT28

#5.0 CHORD PLOTS FOR NEUTROPHILS----
#note; export these graphs as PDF to preserve connecting arrows
##5.1 Signals to neutrophils in S12----
#target: neutrophils (5)
#sources from heatmaps: neu (5)
par(mfrow = c(1,1), xpd=TRUE)
NeuChord_S12 <-netVisual_chord_gene(cellchat_S12, 
                                    sources.use = c(5), 
                                    targets.use = c(5),  
                                    show.legend = TRUE, 
                                    title.name = paste0("Signaling received by Sham Neutrophils at 12hpi ", names(cellchat_S12)))

##5.2 Signals to neutrophils in T12----
#target: neutrophils (5)
#sources from heatmap: alv fib (2)) and acap (14), neu (5)
NeuChord_T12 <-netVisual_chord_gene(cellchat_T12, 
                                    sources.use = c(5,2,14), 
                                    targets.use = c(5),  
                                    show.legend = TRUE, 
                                    title.name = paste0("Signaling received by SCI Neutrophils at 12hpi ", names(cellchat_T12)))

##5.3 Signals to neutrophils in S3----
#target: neutrophils (5)
#sources from heatmap: gcap (13) and acap (14), neu (5), fibros (1, 2)
NeuChord_S3 <-netVisual_chord_gene(cellchat_S3, 
                                    sources.use = c(1, 2, 5, 13,14), 
                                    targets.use = c(5),  
                                    show.legend = TRUE, 
                                    title.name = paste0("Signaling received by Sham Neutrophils at 3dpi ", names(cellchat_S3)))

##5.4 Signals to neutrophils in T3----
#target: neutrophils (5)
#sources from heatmap: fibros (1, 2), gcap (13) and acap (14), neu (5), alv mac (8), tlymph (12), b lymph (7)
NeuChord_T3 <-netVisual_chord_gene(cellchat_T3, 
                                    sources.use = c(1,2,5,7,8,12,13,14), 
                                    targets.use = c(5),  
                                    show.legend = TRUE, 
                                    title.name = paste0("Signaling received by SCI Neutrophils at 3dpi ", names(cellchat_T3)))

##5.5 Signals to neutrophils in S28----
#target: neutrophils (5)
#sources from heatmap: meso (3), mono (6)
NeuChord_S28 <-netVisual_chord_gene(cellchat_S28, 
                                    sources.use = c(3,6), 
                                    targets.use = c(5),  
                                    show.legend = TRUE, 
                                    title.name = paste0("Signaling received by Sham Neutrophils at 28dpi ", names(cellchat_S28)))


##5.6 Signals to neutrophils in T28----
#target: neutrophils (5)
#sources from heatmap: alv mac (8), acap (14)
NeuChord_T28 <-netVisual_chord_gene(cellchat_T28, 
                                    sources.use = c(8,14), 
                                    targets.use = c(5),  
                                    show.legend = TRUE, 
                                    title.name = paste0("Signaling received by SCI Neutrophils at 28dpi ", names(cellchat_T28)))

##6.0 BAR PLOTS FOR TOTAL INTERACTION STRENGTH AND NUMBER----
##6.1a N0 vs S12----
gg1 <- compareInteractions(cellchatN0vS12, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchatN0vS12, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

##6.1b N0 vs T12----
gg3 <- compareInteractions(cellchatN0vT12, show.legend = F, group = c(1,2))
gg4 <- compareInteractions(cellchatN0vT12, show.legend = F, group = c(1,2), measure = "weight")
gg3 + gg4

##6.1c N0 vs S3----
gg5 <- compareInteractions(cellchatN0vS3, show.legend = F, group = c(1,2))
gg6 <- compareInteractions(cellchatN0vS3, show.legend = F, group = c(1,2), measure = "weight")
gg5 + gg6

##6.1d N0 vs T3----
gg7 <- compareInteractions(cellchatN0vT3, show.legend = F, group = c(1,2))
gg8 <- compareInteractions(cellchatN0vT3, show.legend = F, group = c(1,2), measure = "weight")
gg7 + gg8

##6.1e N0 vs S28----
gg9 <- compareInteractions(cellchatN0vS28, show.legend = F, group = c(1,2))
gg10 <- compareInteractions(cellchatN0vS28, show.legend = F, group = c(1,2), measure = "weight")
gg9 + gg10

##6.1f N0 vs T28----
gg11 <- compareInteractions(cellchatN0vT28, show.legend = F, group = c(1,2))
gg12 <- compareInteractions(cellchatN0vT28, show.legend = F, group = c(1,2), measure = "weight")
gg11 + gg12

#7.0 MATCHED CHORD PLOTS FOR NEUTROPHILS----
#note; export these graphs as PDF to preserve connecting arrows
#redone to match custom chord plot inputs
##7.1 To Open saved Cell Chat Objects----
#individual objects
setwd("C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023")
cellchat_N0 <- readRDS("C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023/cellchat_N0_12072023.rds")
cellchat_S3 <- readRDS("C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023/cellchat_S3_12072023.rds")
cellchat_S12 <- readRDS("C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023/cellchat_S12_12072023.rds")
cellchat_S28 <- readRDS("C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023/cellchat_S28_12072023.rds")
cellchat_T3 <- readRDS("C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023/cellchat_T3_12072023.rds")
cellchat_T12 <- readRDS("C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023/cellchat_T12_12072023.rds")
cellchat_T28 <- readRDS("C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023/cellchat_T28_12072023.rds")
#merged object
cellchat <- readRDS("C:/Users/kathe/Dropbox/KAM R Code/Code_Edits_08212023/Lung_11062023/Lung_MITO20_11282023/Cell Chat_12072023/cellchat_12072023.rds")

##7.2 Signals to neutrophils in S12----
#target: neutrophils (5)
#sources from heatmaps: neu (5), t lymph (12), b lymph (7), acap (14), gcap(13)
par(mfrow = c(1,1), xpd=TRUE)
NeuChord_S12_match <-netVisual_chord_gene(cellchat_S12, 
                                          sources.use = c(5, 7, 12, 13, 14), 
                                          targets.use = c(5),  
                                          show.legend = TRUE, 
                                          title.name = paste0("Signaling received by Sham Neutrophils at 12hpi ", names(cellchat_S12)))

##7.3 Signals to neutrophils in T12----
#target: neutrophils (5)
#sources from heatmap: neu (5), t lymph (12), b lymph (7), acap (14), gcap(13)
NeuChord_T12_match <-netVisual_chord_gene(cellchat_T12, 
                                          sources.use = c(5, 7, 12, 13, 14), 
                                          targets.use = c(5),  
                                          show.legend = TRUE, 
                                          title.name = paste0("Signaling received by SCI Neutrophils at 12hpi ", names(cellchat_T12)))

##7.4 Signals to neutrophils in S3----
#target: neutrophils (5)
#sources from heatmap: neu (5), t lymph (12), b lymph (7), acap (14), gcap(13), alv mac (8), mono (6) 
NeuChord_S3_match <-netVisual_chord_gene(cellchat_S3, 
                                         sources.use = c(5,6,7,8,12,13,14), 
                                         targets.use = c(5),  
                                         show.legend = TRUE, 
                                         title.name = paste0("Signaling received by Sham Neutrophils at 3dpi ", names(cellchat_S3)))

##7.5 Signals to neutrophils in T3----
#target: neutrophils (5)
#sources from heatmap: neu (5), t lymph (12), b lymph (7), acap (14), gcap(13), alv mac (8), mono (6)
NeuChord_T3_match <-netVisual_chord_gene(cellchat_T3, 
                                         sources.use = c(5,6,7,8,12,13,14), 
                                         targets.use = c(5),  
                                         show.legend = TRUE, 
                                         title.name = paste0("Signaling received by SCI Neutrophils at 3dpi ", names(cellchat_T3)))

##7.6 Signals to neutrophils in S28----
#target: neutrophils (5)
#sources from heatmap: neu (5), t lymph (12), b lymph (7), acap (14), alv mac (8)
NeuChord_S28_match <-netVisual_chord_gene(cellchat_S28, 
                                          sources.use = c(5,7,8,12,14), 
                                          targets.use = c(5),  
                                          show.legend = TRUE, 
                                          title.name = paste0("Signaling received by Sham Neutrophils at 28dpi ", names(cellchat_S28)))


##7.7 Signals to neutrophils in T28----
#target: neutrophils (5)
#sources from heatmap: neu (5), t lymph (12), b lymph (7), acap (14), alv mac (8)
NeuChord_T28_match <-netVisual_chord_gene(cellchat_T28, 
                                          sources.use = c(5,7,8,12,14),
                                          targets.use = c(5),  
                                          show.legend = TRUE, 
                                          title.name = paste0("Signaling received by SCI Neutrophils at 28dpi ", names(cellchat_T28)))
