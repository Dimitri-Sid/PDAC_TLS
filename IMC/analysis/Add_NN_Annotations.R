
# Add NN Annotations to flowframe object
# Meredith Wetzel

# Calculate the nearest neighbors
rm(list = ls())

library(projectR)
library(Seurat)
library(CoGAPS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(pROC)
library(FlowSOM)
library(rstatix)
library(tidyverse)
library(ggpubr)
library(DT)
library(reshape)

setwd("/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/scripts/")

base_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/"
script_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/scripts/"
data_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/data/"

source(paste0(script_dir,"IMCpipeline_functions_packages.R")) 

output<-readRDS(paste0(script_dir,'output_w15.RDS'))

aggr_flowframe <- readRDS("aggr_flowframe_withClusterCelltypes_Final_Updated.RDS")

names(aggr_flowframe)[names(aggr_flowframe) == 'Number_Neighbors'] <- 'N_N'
names(aggr_flowframe)[names(aggr_flowframe) == 'neighbour_4_CellId1'] <- 'NN_1'
names(aggr_flowframe)[names(aggr_flowframe) == 'neighbour_4_CellId2'] <- 'NN_2'
names(aggr_flowframe)[names(aggr_flowframe) == 'neighbour_4_CellId3'] <- 'NN_3'
names(aggr_flowframe)[names(aggr_flowframe) == 'neighbour_4_CellId4'] <- 'NN_4'
names(aggr_flowframe)[names(aggr_flowframe) == 'neighbour_4_CellId5'] <- 'NN_5'
names(aggr_flowframe)[names(aggr_flowframe) == 'neighbour_4_CellId6'] <- 'NN_6'
names(aggr_flowframe)[names(aggr_flowframe) == 'neighbour_4_CellId7'] <- 'NN_7'
names(aggr_flowframe)[names(aggr_flowframe) == 'neighbour_4_CellId8'] <- 'NN_8'
names(aggr_flowframe)[names(aggr_flowframe) == 'neighbour_4_CellId9'] <- 'NN_9'
names(aggr_flowframe)[names(aggr_flowframe) == 'neighbour_4_CellId10'] <- 'NN_10'

rownames(aggr_flowframe)<-paste(aggr_flowframe$sample_id,aggr_flowframe$CellId,sep="_")

#Purpose: concatenate sample ID with cell ID for a unique identifier
aggr_flowframe$NN_1[aggr_flowframe$N_N > 0] <- paste(aggr_flowframe$sample_id[aggr_flowframe$N_N > 0],aggr_flowframe$NN_1[aggr_flowframe$N_N > 0],sep="_")
aggr_flowframe$NN_2[aggr_flowframe$N_N > 1] <- paste(aggr_flowframe$sample_id[aggr_flowframe$N_N > 1],aggr_flowframe$NN_2[aggr_flowframe$N_N > 1],sep="_")
aggr_flowframe$NN_3[aggr_flowframe$N_N > 2] <- paste(aggr_flowframe$sample_id[aggr_flowframe$N_N > 2],aggr_flowframe$NN_3[aggr_flowframe$N_N > 2],sep="_")
aggr_flowframe$NN_4[aggr_flowframe$N_N > 3] <- paste(aggr_flowframe$sample_id[aggr_flowframe$N_N > 3],aggr_flowframe$NN_4[aggr_flowframe$N_N > 3],sep="_")
aggr_flowframe$NN_5[aggr_flowframe$N_N > 4] <- paste(aggr_flowframe$sample_id[aggr_flowframe$N_N > 4],aggr_flowframe$NN_5[aggr_flowframe$N_N > 4],sep="_")
aggr_flowframe$NN_6[aggr_flowframe$N_N > 5] <- paste(aggr_flowframe$sample_id[aggr_flowframe$N_N > 5],aggr_flowframe$NN_6[aggr_flowframe$N_N > 5],sep="_")
aggr_flowframe$NN_7[aggr_flowframe$N_N > 6] <- paste(aggr_flowframe$sample_id[aggr_flowframe$N_N > 6],aggr_flowframe$NN_7[aggr_flowframe$N_N > 6],sep="_")
aggr_flowframe$NN_8[aggr_flowframe$N_N > 7] <- paste(aggr_flowframe$sample_id[aggr_flowframe$N_N > 7],aggr_flowframe$NN_8[aggr_flowframe$N_N > 7],sep="_")
aggr_flowframe$NN_9[aggr_flowframe$N_N > 8] <- paste(aggr_flowframe$sample_id[aggr_flowframe$N_N > 8],aggr_flowframe$NN_9[aggr_flowframe$N_N > 8],sep="_")
aggr_flowframe$NN_10[aggr_flowframe$N_N > 9] <- paste(aggr_flowframe$sample_id[aggr_flowframe$N_N > 9],aggr_flowframe$NN_10[aggr_flowframe$N_N > 9],sep="_")

#Purpose: match NN1 ID with the rowname to populate a column with each cell type (obj type). Look up neighboring cell's info

aggr_flowframe$N1type <- aggr_flowframe$cluster_annot[match(aggr_flowframe$NN_1,rownames(aggr_flowframe))]
aggr_flowframe$N2type <- aggr_flowframe$cluster_annot[match(aggr_flowframe$NN_2,rownames(aggr_flowframe))]
aggr_flowframe$N3type <- aggr_flowframe$cluster_annot[match(aggr_flowframe$NN_3,rownames(aggr_flowframe))]
aggr_flowframe$N4type <- aggr_flowframe$cluster_annot[match(aggr_flowframe$NN_4,rownames(aggr_flowframe))]
aggr_flowframe$N5type <- aggr_flowframe$cluster_annot[match(aggr_flowframe$NN_5,rownames(aggr_flowframe))]
aggr_flowframe$N6type <- aggr_flowframe$cluster_annot[match(aggr_flowframe$NN_6,rownames(aggr_flowframe))]
aggr_flowframe$N7type <- aggr_flowframe$cluster_annot[match(aggr_flowframe$NN_7,rownames(aggr_flowframe))]
aggr_flowframe$N8type <- aggr_flowframe$cluster_annot[match(aggr_flowframe$NN_8,rownames(aggr_flowframe))]
aggr_flowframe$N9type <- aggr_flowframe$cluster_annot[match(aggr_flowframe$NN_9,rownames(aggr_flowframe))]
aggr_flowframe$N10type <- aggr_flowframe$cluster_annot[match(aggr_flowframe$NN_10,rownames(aggr_flowframe))]

#save the aggr_flowframe with the neighbors:

saveRDS(aggr_flowframe, "aggr_flowframe_withClusterCelltypes_Final_Updated_NNannotated.RDS")

