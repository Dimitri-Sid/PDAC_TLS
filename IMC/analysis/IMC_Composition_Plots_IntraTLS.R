## Composition Plots - IntraTLS

# Meredith Wetzel - Modified from Dimitri
# Updated: August 16, 2024

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
library(viridis)
library(ggpubr)
library(rlang)
set.seed(1234)
sessionInfo()

setwd("/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/scripts/")

base_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/"
script_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/scripts/"
data_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/data/"

source(paste0(script_dir,"IMCpipeline_functions_packages.R")) 
output<-readRDS(paste0(script_dir,'output_w15.RDS'))
aggr_flowframe <- readRDS("aggr_flowframe_withClusterCelltypes_Final_Updated.RDS")

plot_cluster_composition <- function(ggdf, group, grouplabel){
  gp <- ggplot(ggdf, aes(x = group, y = proportion, fill=cluster, order=cluster)) +
    geom_bar(stat = "identity", position="fill", width=0.8) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=5),
          axis.text.y = element_text(color="black"),
          axis.ticks = element_line(size=0.25),
          strip.text=element_text(size=6),
          strip.background = element_rect(fill=NA, color=NA),
          legend.key.size = unit(.75,'lines'),
          legend.text = element_text(size=6),#8
          legend.key = element_rect(fill="white"),
    ) +
    ylab("% of Total Cells")+
    xlab(grouplabel)+
    scale_fill_manual(values = colorassigned,
                      breaks = clusterlevels,
                      labels = clusterlevels)+
    scale_y_continuous(expand = c(0,0))+
    guides(fill=guide_legend(ncol=1))
  return(gp)
}

# add new labels with location and sample
aggr_flowframe$loc_sample <- paste0(aggr_flowframe$Location, "_", aggr_flowframe$sample_id)
aggr_flowframe$loc_sample

#################### BREAK UP DATA INTO SUBGROUPINGS #########################

aggr_IntraTLS <- aggr_flowframe[aggr_flowframe$Location == "Intra",]

aggr_IntraTLS_R <- aggr_IntraTLS[aggr_IntraTLS$Response == "R",]
aggr_IntraTLS_NR <- aggr_IntraTLS[aggr_IntraTLS$Response == "NR",]

#################### Prep Data #########################

## Set up levels
samplevels<- output$meta_data$sample_id
tissuelevels<- output$meta_data$Location
responselevels<- output$meta_data$Response
loc_sample_levels <- paste0(output$meta_data$Location, "_", output$meta_data$sample_id)

output$meta_data$loc_sample <- paste0(output$meta_data$Location, "_", output$meta_data$sample_id)

# Set up the data frame for proportional plotting
counts_table <- table(aggr_IntraTLS$cluster_annot, aggr_IntraTLS$sample_id)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)
ggdf <- melt(data.frame(cluster = rownames(props),props, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion",
             variable.name = "sample_id")

ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdf$Location <- factor(output$meta_data$Location[match(ggdf$sample_id,output$meta_data$sample_id)], levels=unique(tissuelevels))
ggdf$Response <- factor(output$meta_data$Response[match(ggdf$sample_id,output$meta_data$sample_id)], levels=c("NR","R"))
ggdf$Arm <- factor(output$meta_data$Arm[match(ggdf$sample_id,output$meta_data$sample_id)], levels=c("A","B", "C"))
ggdf$loc_sample <- factor(output$meta_data$loc_sample[match(ggdf$sample_id,output$meta_data$sample_id)], levels=sort(loc_sample_levels))

# Set up Colors
clusterlevels <- unique(aggr_IntraTLS$cluster_annot)
colorassigned<-kovesi.diverging_rainbow_bgymr_45_85_c67(length(clusterlevels))
clusternames<-clusterlevels
names(colorassigned)<-clusternames
colorassigned[["DC"]] <- "forestgreen"
colorassigned[["B_GC"]] <- "lightgreen"
colorassigned[["Bcell"]] <- "cyan"
colorassigned[["Myof"]] <- "lightblue"
colorassigned[["Tcell"]] <- "orange"
colorassigned[["B_GC_prolif"]] <- "purple"
colorassigned[["Th"]] <- "blue"
colorassigned[["Epi"]] <- "brown"
colorassigned[["T_prolif"]] <- "red"
colorassigned[["NK"]] <- "pink"

pdf("../figPlots/composition/CompositionPlots_clusterCelltype_IntraTLS.pdf", height = 3, width = 4)
plot_cluster_composition(ggdf, group = ggdf$sample_id, grouplabel = "sample_id")
plot_cluster_composition(ggdf, group = ggdf$Arm, grouplabel = "Arm")
plot_cluster_composition(ggdf, group = ggdf$Location, grouplabel = "Location")
plot_cluster_composition(ggdf, group = ggdf$Response, grouplabel = "Response")
plot_cluster_composition(ggdf, group = ggdf$loc_sample, grouplabel = "Location + SampleID")
dev.off()

saveRDS(ggdf, "../figPlots/composition/IntraTLS_props_ggdf.RDS")
