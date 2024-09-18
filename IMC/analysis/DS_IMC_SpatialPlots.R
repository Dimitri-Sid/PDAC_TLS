#### Final Annotation Spatial Plots
# Meredith Wetzel
# August 14, 20204

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
library(xlsx)
set.seed(1234)
sessionInfo()

setwd("/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/scripts/")

base_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/"
script_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/scripts/"
data_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/data/"

source(paste0(script_dir,"IMCpipeline_functions_packages.R")) 

output<-readRDS(paste0(script_dir,'output_w15.RDS'))

aggr_flowframe <- readRDS("aggr_flowframe_withClusterCelltypes_Final_Updated.RDS")

# Spatial plots with overall cell types

unique(aggr_flowframe$cluster_annot)

aggr_flowframe$annot_color <- "black"
aggr_flowframe$annot_color[aggr_flowframe$cluster_annot == 'DC'] <- 'forestgreen'
aggr_flowframe$annot_color[aggr_flowframe$cluster_annot == 'B_GC'] <- 'lightgreen'
aggr_flowframe$annot_color[aggr_flowframe$cluster_annot == 'Bcell'] <- 'cyan'
aggr_flowframe$annot_color[aggr_flowframe$cluster_annot == 'Myof'] <- 'lightblue'
aggr_flowframe$annot_color[aggr_flowframe$cluster_annot == 'Tcell'] <- 'orange'
aggr_flowframe$annot_color[aggr_flowframe$cluster_annot == 'B_GC_prolif'] <- 'purple'
aggr_flowframe$annot_color[aggr_flowframe$cluster_annot == 'Th'] <- 'blue'
aggr_flowframe$annot_color[aggr_flowframe$cluster_annot == 'T_prolif'] <- 'red'
aggr_flowframe$annot_color[aggr_flowframe$cluster_annot == 'NK'] <- 'pink'
aggr_flowframe$annot_color[aggr_flowframe$cluster_annot == 'Epi'] <- 'brown'

#saveRDS(aggr_flowframe, file = "aggr_flowframe_withClusterCelltypes_Final_Colors.RDS")

unique(aggr_flowframe$annot_color) # make sure there's no black/unlabeled cells

myColors <- aggr_flowframe$annot_color
names(myColors) <- aggr_flowframe$cluster_annot
colScale <- scale_colour_manual(name = "cluster_annot",values = myColors)

# visualize plot

ggplot(aggr_flowframe, aes(x=X_position, y=Y_position, color=cluster_annot)) + geom_point(size=.01) +facet_wrap(~ROI_ID) + 
  guides(colour = guide_legend(override.aes = list(size=3))) +ggtitle("Cluster Annotations by ROI_ID") + 
  colScale + theme(legend.title = element_text(size = 5), legend.text = element_text(size = 4)) + theme_classic()

# Save plot as pdf - updated dimensions, point size, and legend text size

pdf("../figPlots/spatial_plots/updatedCellType_spatialPlots_noLines.pdf", height = 12, width = 12)
ggplot(aggr_flowframe, aes(x=X_position, y=Y_position, color=cluster_annot)) + geom_point(size=.01) +facet_wrap(~ROI_ID) + 
  guides(colour = guide_legend(override.aes = list(size=3))) +ggtitle("Cluster Annotations by ROI_ID") + 
  colScale + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 9)) + theme_classic()
dev.off()

