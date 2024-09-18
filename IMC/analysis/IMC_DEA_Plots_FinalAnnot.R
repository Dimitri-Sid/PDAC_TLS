## DEA Plots
# Meredith Wetzel
# August 14, 2024

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

setwd("/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/scripts/")

base_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/"
script_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/scripts/"
data_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/data/"

source(paste0(script_dir,"IMCpipeline_functions_packages.R")) 

output<-readRDS(paste0(script_dir,'output_w15.RDS'))

#fully annotated/gated object from DS_IMCgating_plots script
aggr_flowframe <- readRDS("aggr_flowframe_withClusterCelltypes_Final_Updated.RDS")

######################### READ FUNCTION ##############################

plot_DE <- function(genes, df, xval, method, comparisons, facet_variable, file_name = "DE_plots.pdf", height = 10, width = 10) {
  plot_list <- list()
  for (gene in genes) {
    p <- ggplot(df, aes(x = !!sym(xval), y = !!sym(gene), fill = !!sym(xval))) +
      geom_boxplot() +
      scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
      theme(plot.title = element_text(size = 11)) +
      xlab("") + 
      facet_wrap(as.formula(paste("~", facet_variable))) +
      stat_compare_means(method = method, comparisons = comparisons, paired = FALSE) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
    plot_list[[gene]] <- p
  }
  pdf(file_name, height = height, width = width)
  for (gene in genes) {
    print(plot_list[[gene]])
  }
  dev.off()
}

##################### SUBSET DATA ##################################

aggr_IntraTLS <- aggr_flowframe[aggr_flowframe$Location == "Intra",]
aggr_PeriTLS <- aggr_flowframe[aggr_flowframe$Location == "Peri",]
aggr_allTLS <- aggr_flowframe[aggr_flowframe$Location == "Intra" | aggr_flowframe$Location == "Peri",]
aggr_LN <- aggr_flowframe[aggr_flowframe$Location == "LN",]

##################### RUN FUNCTION ##################################

genes <- colnames(aggr_flowframe[,2:39])
genes

comparisons <- list(c("R","NR"))

###################### Intra TLS ####################################

plot_DE(genes = genes, df = aggr_IntraTLS, xval = "Response", method = "wilcox.test", comparisons, facet_variable = "cluster_annot", 
        file_name = "../figPlots/DEA/IntraTLS_DEA_WC_Response_Boxplot.pdf", height = 10, width = 10)

###################### Peri TLS ####################################

plot_DE(genes = genes, df = aggr_PeriTLS, xval = "Response", method = "wilcox.test", comparisons, facet_variable = "cluster_annot", 
        file_name = "../figPlots/DEA/PeriTLS_DEA_WC_Response_Boxplot.pdf", height = 10, width = 10)

###################### all TLS ####################################

plot_DE(genes = genes, df = aggr_allTLS, xval = "Response", method = "wilcox.test", comparisons, facet_variable = "cluster_annot", 
        file_name = "../figPlots/DEA/allTLS_DEA_WC_Response_Boxplot.pdf", height = 10, width = 10)

###################### LN ####################################

plot_DE(genes = genes, df = aggr_LN, xval = "Response", method = "wilcox.test", comparisons, facet_variable = "cluster_annot", 
        file_name = "../figPlots/DEA/LN_DEA_WC_Response_Boxplot.pdf", height = 10, width = 10)
