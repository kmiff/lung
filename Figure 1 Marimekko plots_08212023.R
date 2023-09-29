library(Seurat)
library(tidyverse)
library(readxl)
library(viridis)
library(ggrepel)
library(ggpubr)
library(MASS)
library(RColorBrewer)
library(ggmosaic)
library(ggplot2)
library(ggthemes)
#open seurat object
lung <- readRDS("C:/Users/kathe/Dropbox/KAM R Code/SCT_12122022/lung.rds")

#establish functions needed (based off https://github.com/rsankowski/sankowski_et_al_murine_CAMs/blob/main/R/functions.R paper code)

#set Colour Palette
colors_many <- palette.colors(n = NULL, palette = "Polychrome 36", alpha=1, recycle = FALSE)

#Marimekko plot without stats
mosaicGG2 <- function(data, X, FILL, colors = colors_many, rect_col = 'white', line_width = 0.25) {
  require(dplyr)
  require(reshape2)
  #require(ggthemes)
  # Proportions in raw data
  DF <- as.data.frame.matrix(table(data[[X]], data[[FILL]]))
  DF$groupSum <- rowSums(DF)
  DF$xmax <- cumsum(DF$groupSum)
  DF$xmin <- DF$xmax - DF$groupSum
  DF$X <- row.names(DF)
  DF$groupSum <- NULL
  DF_melted <- melt(DF, id = c("X", "xmin", "xmax"), variable.name = "FILL")
  DF_melted <- DF_melted %>%
    group_by(X) %>%
    mutate(ymax = cumsum(value/sum(value)),
           ymin = ymax - value/sum(value))
  
  # Chi-sq test
  results <- chisq.test(table(data[[FILL]], data[[X]])) # fill and then x
  resid <- reshape2::melt(results$residuals)
  names(resid) <- c("FILL", "X", "residual")
  
  # Merge data
  DF_all <- merge(DF_melted, resid)
  
  # Positions for labels
  DF_all$xposn <- DF_all$xmin + (DF_all$xmax - DF_all$xmin)/2
  index <- DF_all$xmax == max(DF_all$xmax)
  #DF_all$yposn <- DF_all$ymin[index] + (DF_all$ymax[index] - DF_all$ymin[index])/2
  yposn = 0
  # Plot
  g <- ggplot(DF_all, aes(ymin = ymin,  ymax = ymax, xmin = xmin,
                          xmax = xmax, fill = FILL)) +
    geom_rect(col = rect_col, lwd = line_width) +
    geom_text(aes(x = xposn, label = X),
              y = 1, size = 3, angle = 90, hjust = 1, show.legend = FALSE,check_overlap = T) +
    geom_text(aes(x = max(xmax),  y = yposn, label = FILL),
              size = 3, hjust = 1, show.legend = FALSE,check_overlap = T) +
    scale_fill_manual(FILL, values = colors) +
    scale_x_continuous(X, expand = c(0,0)) +
    scale_y_continuous("Proportion", expand = c(0,0)) +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(g)
}

# hypergeometric tests (based on https://github.com/rsankowski/sankowski_et_al_murine_CAMs/blob/main/R/functions.R paper code )
# Ranges of P‐values for significantly enriched regions per cluster, based on hypergeometric tests.
# Benjamini–Hochberg correction was applied to correct for multiple testing (*P < 0.05; **P < 0.01; ***P < 0.001).
hyper_test_n <- function(data = df, var1 = "Cluster", var2 = "Region") {
  require(tidyverse) 
  require(broom)
  
  .df <- data.frame()
  for (i in unique(data[[var2]])) {
    data2 <- data
    data2[[var2]] <- factor(ifelse(data2[[var2]] == i, i, paste0("non_",i)), levels = c(i, paste0("non_",i)))
    clusters <- as_tibble(table(data2[[`var1`]]), .name_repair = 'unique')
    colnames(clusters) <- c(var1, 'cluster_size')
    vars <- as_tibble(table(data2[,var1], data2[,var2]), .name_repair = 'unique')
    colnames(vars) <- c(var1, var2, "freq_var2")
    vars_wide <- spread(vars, var2, freq_var2)
    
    vars_df <- vars_wide %>%
      left_join(clusters)
    
    
    #hypergeometric test
    #option a
    test_df<- data.frame(q=vars_df[,i], 
                         m=sum(vars_df[,i]), 
                         n=sum(vars_df[,paste0("non_",i)]),
                         k=vars_df[,4])
    
    colnames(test_df)[1] <- "q"
    
    p_hyper <- apply(test_df, MARGIN = 1, function(x) 1-phyper(max(0,x[[1]]-1), x[[2]], x[[3]], x[[4]])) #probability to get q or more successes in populaton
    
    test_df$p_hyper <- p_hyper
    test_df$Cluster <- vars_df[[`var1`]]
    test_df$enrichment_var <- i
    .df <- .df %>%
      bind_rows(test_df[,c("q","m","n","cluster_size","p_hyper","Cluster","enrichment_var")])
  }
  
  
  .df$padj <- p.adjust(.df$p_hyper, method="BH")
  .df$Significance <- ifelse(.df$padj<0.05 & .df$padj>0.01, '*',
                             ifelse(.df$padj<0.01 & .df$padj>0.001, '**',
                                    ifelse(.df$padj<0.001, '***','n.s.')))
  
  return(.df)
}


#extract metadata
metadata <- lung@meta.data

#Fig 1 Marimekko barplot & hypergeometric tests
mosaicGG2(metadata, "orig.ident", "main.cluster") 

hyper_test_n(metadata, "orig.ident", "main.cluster") %>%
  write_csv("marimekko_clusters_condition_all_cells.csv")

