#------------------------------------------------All conditions_12222023
#load packages
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
setwd("C:/Users/kathe/Dropbox/KAM R Code/FB_Oligo_CellChat_032024")

#--------preprocessing data
#load our lung data
combine.combined <- load("C:/Users/kathe/Dropbox/KAM R Code/FB_Oligo_CellChat_032024/combined.rda")
levels(combine.combined)
combine.combined[["main.cluster"]] <- Idents(combine.combined)
Idents(combine.combined)<- "main.cluster"
levels(combine.combined)
split <- SplitObject(combine.combined, split.by = "orig.ident")
#preprocessing data from each condition
meta_list<-list()
count_list<-list()
for(i in 1:6){
  object <-split[[i]]
  DefaultAssay(object) <- 'RNA'
  Idents(object)<- 'main.cluster'
  #count data
  counts_object <- GetAssayData(object, assay = "RNA", slot = "counts") # raw counts
  #normalize counts for cell chat
  counts_object <-normalizeData(counts_object)
  meta_object <-object@meta.data
  #not all groups have Mixed Endo pop, so need to exclude mismatched levels
  meta_object$main.cluster <- droplevels(meta_object$main.cluster, exclude = setdiff(levels(meta_object$main.cluster),unique(meta_object$main.cluster)))
  meta_object$main.cluster.f<-gsub(" ", "-", meta_object$main.cluster)
  meta_list[[i]]<-meta_object
  count_list[[i]]<-counts_object
}



#----------------------run functions
id <- 0
df_list <- c() # the list to hold the CellChat object
# run CellChat for all datasets
{
  #ctrl.00d 
  net.ctrl.00d  <- use_CellChat(count_list[[1]], meta_list[[1]],"main.cluster.f")
  id <- id + 1
  df_list[[id]] <- 'net.ctrl.00d'
  #plx.07d
  net.plx.07d <- use_CellChat(count_list[[2]], meta_list[[2]],"main.cluster.f")
  id <- id + 1
  df_list[[id]] <- 'net.plx.07d'
  #plx.00d
  net.plx.00d <- use_CellChat(count_list[[3]], meta_list[[3]],"main.cluster.f")
  id <- id + 1
  df_list[[id]] <- 'net.plx.00d'
  #ctrl.28d
  net.ctrl.28d <- use_CellChat(count_list[[4]], meta_list[[4]],"main.cluster.f")
  id <- id + 1
  df_list[[id]] <- 'net.ctrl.28d'
  #ctrl.07d
  net.ctrl.07d <- use_CellChat(count_list[[5]], meta_list[[5]],"main.cluster.f")
  id <- id + 1
  df_list[[id]] <- 'net.ctrl.07d'
  #plx.28d
  net.plx.28d <- use_CellChat(count_list[[6]], meta_list[[6]],"main.cluster.f")
  id <- id + 1
  df_list[[id]] <- 'net.plx.28d'

}

# convert the CCIs into the format to plot chord diagrams
dir<-"/fs/ess/PCON0022/guoqi/Lung_sci_katherine/Results/Cell_plot_1222"
for (i in 1:length(df_list)) {
  condition <- df_list[[i]] # the experiment condition
  prefix <- paste0(gsub(pattern = '^net.', replacement = '', df_list[[i]]))
  save_data(dir,get(df_list[i][[1]]), prefix)
}

# plot diagrams in batch
mat.list <- Sys.glob("/fs/ess/PCON0022/guoqi/Lung_sci_katherine/Results/Cell_plot_1222/*_mat.csv") # the list of adjacency matrix files
cat.list <- Sys.glob("/fs/ess/PCON0022/guoqi/Lung_sci_katherine/Results/Cell_plot_1222/*_cat.csv") # the list of the gene category files, i.e., cell types
for (i in 1:length(mat.list)) {
  print(mat.list[i])
  a<-genChord(mat = mat.list[i], cat = cat.list[i],
              prefix = strsplit(x = mat.list[i], split = '_')[[1]][1])
}

#---------function
#-------------the function to get L-R results from cellchat package
use_CellChat <- function(counts, meta,celltype_column) {
  # counts : the gene expression matrix
  # meta   : the meta data matched with the gene expression matrix
  # celltype_column.  : cell type column in the metadata
  counts <- as.matrix(counts) # transform data frame into matrix
  cellchat <- createCellChat(object = counts, meta = meta, group.by = celltype_column) # create a CellChat object
  levels(cellchat@idents) # check the cell types
  CellChatDB <- CellChatDB.mouse # set the database according to the organism
  dplyr::glimpse(CellChatDB$interaction) # show the structure of the database
  unique(CellChatDB$interaction$annotation) # check the annotation of interactions
  CellChatDB.use <- CellChatDB # use all databases
  cellchat@DB <- CellChatDB.use # set the used database in the object
  cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.mouse)
  cellchat <- computeCommunProb(cellchat, raw.use = T, nboot = 1, population.size = T,trim = 0.1,type = "truncatedMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat) # get the pairs of cell types or ligand receptors in the interactions
  cell_type_pair <- paste0(df.net$source, '_', df.net$target) # build the pairs of interaction cell types
  mean_exp<-c()
  mean_exp <- sapply(1:nrow(df.net), function(x) {
    ct1 <- df.net$source[x] # source cell typeQ
    ct2 <- df.net$target[x] # target cell type
    cells1 <- rownames(meta[meta[[celltype_column]] == ct1,]) # cells belonging to the source cell type
    cells2 <- rownames(meta[meta[[celltype_column]] == ct2,]) # cells belonging to the target cell type
    #character "_" split;then;character including "-", do not change anything
    gene_pair1 <- df.net$ligand[x] # genes corresponding to the ligand
    if(grepl("_",gene_pair1)){
      #split
      gene_pair1<-strsplit(gene_pair1, split = '_')[[1]]
      #do not change anything if "-" existing
      genes1_unchange<-gene_pair1[grepl("-",gene_pair1)]
      genes1_lower<-capitalize(tolower(gene_pair1[!grepl("-",gene_pair1)]))
      genes1<-c(genes1_unchange,genes1_lower)
    }
    else{
      genes1=gene_pair1
    }
    
    gene_pair2 <- df.net$receptor[x] # genes corresponding to the receptor
    if(grepl("_",gene_pair2)){
      #split
      gene_pair2<-strsplit(gene_pair2, split = '_')[[1]]
      #do not change anything if "-" existing
      genes2_unchange<-gene_pair2[grepl("-",gene_pair2)]
      genes2_lower<-capitalize(tolower(gene_pair2[!grepl("-",gene_pair2)]))
      genes2<-c(genes2_unchange,genes2_lower)
    }
    else{
      genes2=gene_pair2
    }
    if (length(genes2) > 1 & genes2[2] == 'R2') {
      stop("I'm an error")
      genes2[2] <- gsub(pattern = '1', replacement = '2', genes2[1])
    }
    cat("Ligand genes: ", genes1, "\n")
    cat("Receptor genes: ", genes2, "\n\n")
    mean1 <- mean(counts[genes1, cells1]) # the mean expression value of the ligand within the source cell type
    mean2 <- mean(counts[genes2, cells2]) # the mean expression value of the receptor within the source cell type
    avg <- (mean1 + mean2) / 2 # the mean expression value of the ligands and receptors
    return(avg)
  })
  df.net <- cbind(cell_type_pair, mean_exp, df.net) # add the mean expression and cell type pairs
  return(df.net)
}
#------------the function to save the adjacency matrix and L-R categories
save_data <- function(dir,net.df, prefix) {
  # dir.    : directory for the results
  # net.df  : the data frame to save the CCIs
  # prefix  : the prefix of the files to output
  setwd(dir)
  node_v1 <- vector() # the first node
  node_v2 <- vector() # the second node
  weight_v <- vector() # the edge weight
  
  for (i in 1:nrow(net.df)) {
    cci <- strsplit(net.df$interaction_name_2[i], split = ' ')[[1]] # split the interaction name
    ct.pair <- strsplit(net.df$cell_type_pair[i], split = '_')[[1]] # get the cell type pair
    ct1 <- ct.pair[1]
    ct2 <- ct.pair[2]
    prefix1 <- ct1
    prefix2 <- ct2
    node1 <- paste0(prefix1, '_', cci[1])
    id2 <- cci[length(cci)]
    weight <- net.df$prob[i]
    if (length(grep('\\+', id2)) > 0) {
      array <- strsplit(id2, split = '\\+')[[1]]
      node2 <- paste0(prefix2, '_', substr(array[1], 2, nchar(array[1])))
      node3 <- paste0(prefix2, '_', substr(array[2], 1, nchar(array[2]) - 1))
      node_v1 <- c(node_v1, node1, node1)
      node_v2 <- c(node_v2, node2, node3)
      weight_v <- c(weight_v, weight, weight)
    } else {
      node2 <- paste0(prefix2, '_', id2)
      node_v1 <- c(node_v1, node1)
      node_v2 <- c(node_v2, node2)
      weight_v <- c(weight_v, weight)
    }
  }
  
  if (length(node_v1) <= 0) {
    return()
  }
  
  g_df <- data.frame(node_v1, node_v2, weight_v)
  g <- graph.data.frame(g_df, directed = T)
  E(g)$weight <- g_df[[3]]
  adj <- get.adjacency(g, attr = 'weight')
  write.csv(x = as.data.frame(as.matrix(adj)), file = paste0(prefix, '_mat.csv'))
  
  Genes <- names(V(g))
  arrays <- strsplit(Genes, split = '_')
  ID <- sapply(arrays, function(x) {
    return(x[1])
  })
  category <- data.frame(Genes = Genes, ID = ID)
  write.csv(x = category, file = paste0(prefix, '_cat.csv'))
}
#------------the function to save the adjacency matrix and L-R categories
save_data <- function(dir,net.df, prefix) {
  # dir.    : directory for the results
  # net.df  : the data frame to save the CCIs
  # prefix  : the prefix of the files to output
  setwd(dir)
  node_v1 <- vector() # the first node
  node_v2 <- vector() # the second node
  weight_v <- vector() # the edge weight
  
  for (i in 1:nrow(net.df)) {
    cci <- strsplit(net.df$interaction_name_2[i], split = ' ')[[1]] # split the interaction name
    ct.pair <- strsplit(net.df$cell_type_pair[i], split = '_')[[1]] # get the cell type pair
    ct1 <- ct.pair[1]
    ct2 <- ct.pair[2]
    prefix1 <- ct1
    prefix2 <- ct2
    node1 <- paste0(prefix1, '_', cci[1])
    id2 <- cci[length(cci)]
    weight <- net.df$prob[i]
    if (length(grep('\\+', id2)) > 0) {
      array <- strsplit(id2, split = '\\+')[[1]]
      node2 <- paste0(prefix2, '_', substr(array[1], 2, nchar(array[1])))
      node3 <- paste0(prefix2, '_', substr(array[2], 1, nchar(array[2]) - 1))
      node_v1 <- c(node_v1, node1, node1)
      node_v2 <- c(node_v2, node2, node3)
      weight_v <- c(weight_v, weight, weight)
    } else {
      node2 <- paste0(prefix2, '_', id2)
      node_v1 <- c(node_v1, node1)
      node_v2 <- c(node_v2, node2)
      weight_v <- c(weight_v, weight)
    }
  }
  
  if (length(node_v1) <= 0) {
    return()
  }
  
  g_df <- data.frame(node_v1, node_v2, weight_v)
  g <- graph.data.frame(g_df, directed = T)
  E(g)$weight <- g_df[[3]]
  adj <- get.adjacency(g, attr = 'weight')
  write.csv(x = as.data.frame(as.matrix(adj)), file = paste0(prefix, '_mat.csv'))
  
  Genes <- names(V(g))
  arrays <- strsplit(Genes, split = '_')
  ID <- sapply(arrays, function(x) {
    return(x[1])
  })
  category <- data.frame(Genes = Genes, ID = ID)
  write.csv(x = category, file = paste0(prefix, '_cat.csv'))
}
#-------------the function to plot a chord diagram for a single datasete
genChord <- function(mat, cat, prefix) {
  # mat    : the adjacency matrix for the CCIs
  # cat   : the category for each gene, i.e., cell type
  # prefix : the prefix of the figure and legend files
  library(tidyverse)
  library(circlize)
  library(ggsci)
  library(igraph)
  library(gtools)
  library(ComplexHeatmap)
  # read data
  graph_adj <- read.csv(mat, row.names = 1)
  graph_module <- read.csv(cat, row.names = 1)
  g <- graph.adjacency(as.matrix(graph_adj), weighted = T)
  graph_module <- graph_module %>%
    mutate(color = as_factor(ID))%>%  
    mutate(color = fct_recode(color ,
                             # "#5A5156FF" = "Leptomeningeal-cells",
                             # "#E4E1E3FF" = "Erythroid-cells",
                              "#F6222EFF" = "Oligodendrocyte-lineage",
                              #"#FE00FAFF" = "Neutrophils",
                              #"#16FF32FF" = "Intermediate-progenitors",
                              #"#3283FEFF" = "B-cells",
                             # "#FEAF16FF" = "Ependymal-cells",
                             # "#B00068FF" = "T-cells",
                             # "#1CFFCEFF" = "Astrocytes",
                             # "#90AD1CFF" = "Monocytes",
                             # "#2ED9FFFF" = "Endothelial-cells",
                             # "#DEA0FDFF" = "MDMs",
                              "#AA0DFEFF" = "Microglia",
                             #just want oligo microglia interactions
    ))
  
  # filter graph by weights
  g <- graph.adjacency(as.matrix(graph_adj), weighted = T)
  
  raw_edges <-  as.data.frame(cbind(get.edgelist(g), E(g)$weight)) %>%
    mutate(
      V1 = gsub('\\.', '-', V1),
      V2 = gsub('\\.', '-', V2),
      V3 = as.numeric(V3),
      V4 = 1
    )
  edges <- raw_edges %>%
    arrange(V3)
  
  # Normalize the weight score to quantiles
  # quartiles_weight <- quantcut(edges$V3, 4)
  # levels(quartiles_weight) <- c(1:4)
  # edges$V3 <- as.integer(quartiles_weight)
  
  
  nodes <-  unique(c(edges$V1, edges$V2))
  
  sectors <- sort(unique(c(raw_edges$V1, raw_edges$V2)))
  
  # diagram colors
  grid_col <- graph_module %>%
    dplyr::filter(Genes %in% nodes) %>% 
    dplyr::select(Genes, color) %>%
    mutate(color = as.character(color)) %>%
    deframe()
  
  #https://vuetifyjs.com/en/styles/colors/#material-colors
  col_fun = colorRamp2(range(edges$V3), c("#FFFDE7", "#013220"))
  
  
  circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(-0.15,0.2))
  circos.initialize(sectors, xlim = c(0, 1))
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.border = NA)
  
  # we go back to the first track and customize sector labels
  circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.name = get.cell.meta.data("sector.index")
      xlim = get.cell.meta.data("xlim")
      this_node_text_color <- graph_module %>%
        dplyr::filter(Genes == sector.name) %>%
        pull(color) %>%
        as.character()
      
      circos.rect(
        xlim[1],
        0,
        xlim[2],
        1,
        col = this_node_text_color,
        border = NA
      )
      
      circos.text(
        mean(xlim),
        2,
        CELL_META$sector.index,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        col= this_node_text_color
      )
    },
    bg.border = NA
  ) 
  
  for (i in seq_len(nrow(edges))) {
    link <- edges[i,]
    circos.link(link[[1]],
                c(0, 1),
                link[[2]],
                c(0, 1),
                col = col_fun(link[[3]]),
                border = NA)
  }
  
  
  # plot legend
  #lgd <- Legend(title ="Score", col_fun = col_fun)
  #grid.draw(lgd)
  #
  #png(
  #  paste(paste0(prefix, "_legend.png"), sep = ""),
  #  width = 1000,
  #  height = 1000,
  #  res = 300
  #)
  #grid.draw(lgd)
  #dev.off()
  
}




# N0 <-split_lung$N0
# DefaultAssay(N0) <- 'RNA'
# Idents(N0)<- 'main.cluster'
# levels(N0)
# #count data
# counts_N0 <- GetAssayData(N0, assay = "RNA", slot = "counts") # raw counts
# #write.table(counts_N0, file = 'counts_N0.csv', sep =',')
# #normalize counts for cell chat
# counts_N0 <-normalizeData(counts_N0)
# #meta data
# meta_N0 <-N0@meta.data
# meta_N0$main.cluster <- droplevels(meta_N0$main.cluster, exclude = setdiff(levels(meta_N0$main.cluster),unique(meta_N0$main.cluster)))
# #switch blank " " to "-"
# meta_N0$main.cluster.f<-gsub(" ", "-", meta_N0$main.cluster)
# #S12
# S12 <-split_lung$S12
# DefaultAssay(S12) <- 'RNA'
# Idents(S12)<- 'main.cluster'
# #count data
# counts_S12 <- GetAssayData(S12, assay = "RNA", slot = "counts") # raw counts
# #normalize counts for cell chat
# counts_S12 <-normalizeData(counts_S12)
# meta_S12 <-S12@meta.data
# #not all groups have Mixed Endo pop, so need to exclude mismatched levels
# meta_S12$main.cluster <- droplevels(meta_S12$main.cluster, exclude = setdiff(levels(meta_S12$main.cluster),unique(meta_S12$main.cluster)))
# meta_S12$main.cluster.f<-gsub(" ", "-", meta_S12$main.cluster)
# #S18