library(dendextend)
library(circlize)
library(stringdist)
library(dplyr)
library(ggplot2)
library(viridis)
#rm(list = ls())

plot_dendrogram <- function(file_path, clone_id, num_colors) {
  #a <- load('C:/Users/sunny/Downloads/22_06_01_pdac-visium-v1.rda')
  #a <- load('C:/Users/sunny/Downloads/changeo_annotated.rda')
  a <- load(file_path)
  ##igh <- res[res$locus == 'IGH',]
  #igk_df <- changeo[changeo$locus == "IGK",]
  #igk_6013_df <- igk_df[igk_df$clone_id == "IGK_7199",]
  cid_df <- changeo[changeo$clone_id == clone_id,]
  # Take counts and cdr3 to collapse rows
  #igk_6013_subset <- igk_6013_df[, c('cdr3', 'count')]
  cid_subset <- cid_df[, c('cdr3', 'count')]
  #igk_6013_collapsed <- aggregate(.~cdr3,data=igk_6013_subset, FUN=sum )
  cid_collapsed <- aggregate(.~cdr3,data=cid_subset, FUN=sum )
  
  #d <- data.frame(matrix(ncol=(nrow(igk_6013_collapsed) + 1), nrow=nrow(igk_6013_collapsed)))
  d <- data.frame(matrix(ncol=(nrow(cid_collapsed) + 1), nrow=nrow(cid_collapsed)))
  ##cols <- igk_6013_df$cdr3_aa
  #cols <- igk_6013_collapsed$cdr3
  cols <- cid_collapsed$cdr3
  ##colnames(d) <- c('cdr3_aa', cols)
  ##d$cdr3_aa <- igk_6013_df$cdr3_aa
  
  colnames(d) <- c('cdr3', cols)
  #d$cdr3 <- igk_6013_collapsed$cdr3
  d$cdr3 <- cid_collapsed$cdr3
  for (i in 2:ncol(d)) {
    d[, i] <- stringdist(cols[i - 1], d[, 1], method='lv')
  }
  # Load distance matrix
  #d <- dist('dist matrix')
  
  dist <- as.dist(d[, 2:ncol(d)])
  
  #igk$cdr3[igk$barcode == 'AATTGCAGCAATCGAC']
  # Create distance matrix using hierarchical clustering
  hc <- as.dendrogram(hclust(dist))
  
  # Perform labeling of counts, branches, etc.
  hc <- color_branches(hc)
  hc <- hc %>% set('branches_lwd', 3)
  
  #counts = igk_6013_collapsed$count
  counts = cid_collapsed$count
#  counts = c(1, 30, 50, 75, 100, 120, 150, 170, 200, 210)
  
  #breaks = seq(min(counts), max(counts), length.out=11)
  breaks = seq(min(counts), max(counts), length.out=num_colors+1)
  #colors <- viridis(10)
  colors <- viridis(num_colors)
  labels <- hc %>% labels
  label_color <- c()
  label_ct <- 1
  for (label in labels) {
    color <- NULL
    break_ct <- 2
    while (is.null(color)) {
      if (counts[label_ct] <= breaks[break_ct]) {
        color <- colors[break_ct - 1]
      }
      break_ct <- break_ct + 1
    }
    label_color <- append(label_color, color)
    label_ct <- label_ct + 1
  }
  cid_collapsed$color <- label_color
  
  hc <- hc %>% set('leaves_pch', 19) %>% set('leaves_cex', 2) %>% set('leaves_col', label_color)
  # Circular dendrogram
  #hcircle <- circlize_dendrogram(hc)
  
  #pdf('changeo_dendrogram.pdf')
  circlize_dendrogram(hc)
  return(list('dend' = hc, 'cid_collapsed' = cid_collapsed))
}


#op <- plot_dendrogram('C:/Users/sunny/Downloads/changeo_annotated.rda', 'IGK_7199', 10)
