### Heatmap with updated labels
# Meredith Wetzel

rm(list = ls())

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(tidyverse)
library(DT)
library(lme4)
library(lmerTest)
library(corrplot)
library(ComplexHeatmap)
library(viridis)
library(RColorBrewer)

setwd("/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/scripts/")

base_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/"
script_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/scripts/"
data_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/data/"

source(paste0(script_dir,"IMCpipeline_functions_packages.R")) 

output<-readRDS(paste0(script_dir,'output_w15.RDS'))

#fully annotated/gated object from DS_IMCgating_plots script
aggr_flowframe <- readRDS("aggr_flowframe_withClusterCelltypes_Final_Updated.RDS")

colnames(aggr_flowframe)

#create aggregated table, grouped by cluster annotation
expr_mean_cluster_tbl <- data.frame(sample_id = aggr_flowframe[,"cluster_annot"], aggr_flowframe) %>%
  group_by(sample_id) %>%  summarize_all(funs(mean))

rownames(expr_mean_cluster_tbl) <- expr_mean_cluster_tbl$sample_id
dfAnnot_cluster <- rownames(expr_mean_cluster_tbl)

#remove the sample_id column 
expr_mean_cluster_tbl <- expr_mean_cluster_tbl[,3:40]
rownames(expr_mean_cluster_tbl) <- dfAnnot_cluster

clust <- dfAnnot_cluster
names(clust) <- c("lightgreen","purple","cyan","forestgreen","brown","lightblue","pink","red","orange","blue")

color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
expr_mean_cluster_tbl_fix <- expr_mean_cluster_tbl[clust,] #order the df by the colors above
rownames(expr_mean_cluster_tbl_fix) <- clust

clust_colors <- setNames(names(clust), clust)
#create another vector with white for the text labels to be on
clust2 <- dfAnnot_cluster
names(clust2) <- c('white','white','white','white','white','white','white','white','white','white')
clust_colors2 <- setNames(names(clust2), clust2)

# Create the annotation with the color mapping
# counts of cells per cluster for barplot in heatmap annot
annot_count <- table(aggr_flowframe$cluster_annot)
cell_props <- (annot_count / sum(annot_count)) * 100
annot_count <- annot_count[dfAnnot_cluster] #so order matches clust/coloring/heatmap
cell_props <- cell_props[dfAnnot_cluster] #so order matches clust/coloring/heatmap

HA2 <- HeatmapAnnotation(
  Cluster_Name = anno_text(as.character(clust), gp = gpar(fontsize = 10)),
  " " = anno_simple(clust, col = clust_colors), "% of Cells" = anno_barplot(
    x = as.numeric(cell_props),
    gp = gpar(fill = 'grey'),
    bar_width = 0.8,
    border = TRUE
  )
)

# Generate the heatmap
Heatmap(
  t(apply(expr_mean_cluster_tbl_fix, 2, scale)),
  top_annotation = HA2,
  clustering_distance_rows = 'pearson',
  clustering_distance_columns = 'pearson',
  col = color,
  name = "cluster expression"
)

pdf("../figPlots/heatmaps/ClusterCelltypes_final_scaled_heatmap.pdf", height = 8, width = 8)
Heatmap(
  t(apply(expr_mean_cluster_tbl_fix, 2, scale)),
  top_annotation = HA2,
  clustering_distance_rows = 'pearson',
  clustering_distance_columns = 'pearson',
  col = color,
  name = "cluster expression"
)
dev.off()
