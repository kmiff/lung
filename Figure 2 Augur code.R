#set working directory
setwd("C:/Users/kathe/Dropbox/KAM R Code/Augur")

#open necessary libraries
library(Seurat)
library(Augur)
library(sparseMatrixStats)
library(tidyverse)
library(viridis)
library(magrittr)

#open Seurat Object
lung <- readRDS("C:/Users/Kathe/Dropbox/KAM R Code/SCT_12122022/lung.rds")
lung
View(lung@meta.data)
#run Augur on default settings
augur = calculate_auc(lung, meta, cell_type_col = "main.cluster", label_col = "orig.ident")
augur$AUC
# A tibble: 17 x 2
cell_type                  auc
<fct>                    <dbl>
  1 Type 1 Neutrophils       0.672
2 Col14a1+ Fibroblasts     0.660
3 Type 2 Neutrophils       0.659
4 Col13a1+ Fibroblasts     0.657
5 Mesothelial              0.645
6 Vein Endothelia          0.645
7 Capillary Aerocytes      0.644
8 General Capillary        0.642
9 Monocytes                0.637
10 NK Cells                 0.617
11 Pericytes                0.615
12 Alveolar Macrophages     0.608
13 Basophils                0.604
14 T Cells                  0.601
15 Interstitial Macrophages 0.601
16 B Cells                  0.599
17 Club Cells               0.594


multiclass_augur = calculate_auc(lung, cell_type_col = "main.cluster", label_col = "orig.ident")
print(multiclass_augur$AUC)
plot_umap(multiclass_augur, lung, cell_type_col = "main.cluster")
plot_lollipop(multiclass_augur)
saveRDS(multiclass_augur, file = "C:/Users/Kathe/Dropbox/KAM R Code/Augur/multiclass_augur.rds")

# Begin by specifying the comparisons of interest
comparisons = list(
  c("N0", "S12"),
  c("N0", "S3"),
  c("N0", "S28"), 
  c("N0", "T12"), 
  c("N0" ,"T3"),
  c("N0" ,"T28"), 
  c("S12" ,"T12"), 
  c("S3" ,"T3"),
  c("S28" ,"T28"))
# Create a list to store Augur results
pairwise_results = list()
# Run Augur on each comparison
for (comparison in comparisons) {
  # Subset the Seurat object to these two experimental conditions only
  Idents(lung) = lung$orig.ident
  sub = subset(lung, idents = comparison)
  # Run Augur
  augur = calculate_auc(sub,cell_type_col = "main.cluster", label_col = "orig.ident" )
  # Store the results
  comparison_name = paste(comparison, collapse = ":")
  pairwise_results[[comparison_name]] = augur
}
pairwise_aucs = pairwise_results %>%
  map(extract2, "AUC") %>%
  bind_rows(.id = "comparison")
print(pairwise_aucs)
saveRDS(pairwise_results, file = "C:/Users/Kathe/Dropbox/KAM R Code/Augur/pairwise_results.rds")
#individual plots
lollipop1 <-plot_lollipop(pairwise_results[["N0:S12"]])
print(lollipop1 + ggtitle("N0vS12"))
lollipop2 <-plot_lollipop(pairwise_results[["N0:S3"]])
print(lollipop2 + ggtitle("N0vS3"))
lollipop3 <-plot_lollipop(pairwise_results[["N0:S28"]])
print(lollipop3 + ggtitle("N0vS28"))
lollipop4 <-plot_lollipop(pairwise_results[["N0:T12"]])
print(lollipop4 + ggtitle("N0vT12"))
lollipop5 <-plot_lollipop(pairwise_results[["N0:T3"]])
print(lollipop5 + ggtitle("N0vT3"))
lollipop6 <-plot_lollipop(pairwise_results[["N0:T28"]])
print(lollipop6 + ggtitle("N0vT28"))
lollipop7 <-plot_lollipop(pairwise_results[["S12:T12"]])
print(lollipop7 + ggtitle("S12vT12"))
lollipop8 <-plot_lollipop(pairwise_results[["S3:T3"]])
print(lollipop8 + ggtitle("S3vT3"))
lollipop9 <-plot_lollipop(pairwise_results[["S28:T28"]])
print(lollipop9 + ggtitle("S28vT28"))