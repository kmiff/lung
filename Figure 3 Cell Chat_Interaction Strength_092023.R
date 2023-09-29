#Load Libraries & Set Directories
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(igraph)

setwd("C:/Users/Kathe/Dropbox/KAM R Code/SCT_12122022/CellChat")

#get seurat object 
lung <- readRDS("C:/Users/Katherine/Dropbox/KAM R Code/SCT_12122022/lung.rds")
lung

#updating database because of missing genes
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
setwd("C:/Users/Kathe/Documents/R/win-library/4.2/CellChat")
setwd ("C:/Users/Kathe/AppData/Local/R/win-library/4.2/CellChat") # This is the folder of CellChat package downloaded from Github
CellChatDB.mouse <- CellChatDB
usethis::use_data(CellChatDB.mouse, overwrite = TRUE)

#split lung, convert to cell chat object & run cell chat on EACH OBJECT before merging together
split_lung <- SplitObject(lung, split.by = "orig.ident")
N0 <-split_lung$N0
DefaultAssay(N0) <- 'RNA'
Idents(N0)<- 'main.cluster'
levels(N0)
data.input_N0 <- GetAssayData(N0, assay = "RNA", slot = "data") # normalized data matrix
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
saveRDS(cellchat_N0, file = "C:/Users/Kathe/Dropbox/KAM R Code/SCT_12122022/CellChat/cellchat_N0.rds")

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
saveRDS(cellchat_S12, file = "C:/Users/Kathe/Dropbox/KAM R Code/SCT_12122022/CellChat/cellchat_S12.rds")

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
saveRDS(cellchat_T12, file = "C:/Users/Kathe/Dropbox/KAM R Code/SCT_12122022/CellChat/cellchat_T12.rds")

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
saveRDS(cellchat_S3, file = "C:/Users/Kathe/Dropbox/KAM R Code/SCT_12122022/CellChat/cellchat_S3.rds")

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
saveRDS(cellchat_T3, file = "C:/Users/Kathe/Dropbox/KAM R Code/SCT_12122022/CellChat/cellchat_T3.rds")

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
saveRDS(cellchat_S28, file = "C:/Users/Kathe/Dropbox/KAM R Code/SCT_12122022/CellChat/cellchat_S28.rds")

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
saveRDS(cellchat_T28, file = "C:/Users/Kathe/Dropbox/KAM R Code/SCT_12122022/CellChat/cellchat_T28.rds")

# since S3 has one extra cell group (mixed endothelial), need to 
#lift up cellchat objects to this one and then merge into one cell object
## See: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.html#load-cellchat-object-of-each-dataset

# Define the cell labels to lift up
group.new = levels(cellchat_S3@idents)
cellchat_N0 <- liftCellChat(cellchat_N0, group.new)
cellchat_S12 <- liftCellChat(cellchat_S12, group.new)
cellchat_T12 <- liftCellChat(cellchat_T12, group.new)
cellchat_T3 <- liftCellChat(cellchat_T3, group.new)
cellchat_S28 <- liftCellChat(cellchat_S28, group.new)
cellchat_T28 <- liftCellChat(cellchat_T28, group.new)

#Merge Data Sets
object.list <- list(N0 = cellchat_N0, S12 = cellchat_S12, T12 = cellchat_T12, S3 = cellchat_S3,
                    T3 = cellchat_T3, S28 = cellchat_S28,T28 = cellchat_T28)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
cellchat
saveRDS(cellchat, file = "C:/Users/Kathe/Dropbox/KAM R Code/SCT_12122022/CellChat/cellchat.rds")
#Merge Data Sets

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
