## DEA Statistics
# Meredith Wetzel
# August 19, 2024

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
library(xlsx)
set.seed(1234)

setwd("/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/scripts/")

base_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/"
script_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/scripts/"
data_dir<-"/Users/meredithwetzel/Library/CloudStorage/OneDrive-JohnsHopkins/DS_IMC_TLS/data/"

source(paste0(script_dir,"IMCpipeline_functions_packages.R")) 
output<-readRDS(paste0(script_dir,'output_w15.RDS'))
aggr_flowframe <- readRDS("aggr_flowframe_withClusterCelltypes_Final_Updated.RDS")

######################### READ FUNCTION ##############################

get_WC_stats <- function(genes, df, group, comparison) {
  stat_list <- list()
  for (gene in genes) {
    formula <- as.formula(paste(gene, "~", comparison))
    stats <- df %>%
      group_by(!!sym(group)) %>%
      wilcox_test(formula) %>%
      adjust_pvalue(method = "fdr") %>%
      add_significance("p.adj")
    stats_df <- as.data.frame(stats)
    stat_list[[gene]] <- stats_df
  }
  return(stat_list)  # Returns the list of stats dfs
}

##################### SUBSET DATA ##################################

aggr_IntraTLS <- aggr_flowframe[aggr_flowframe$Location == "Intra",]
aggr_PeriTLS <- aggr_flowframe[aggr_flowframe$Location == "Peri",]
aggr_allTLS <- aggr_flowframe[aggr_flowframe$Location == "Intra" | aggr_flowframe$Location == "Peri",]
aggr_LN <- aggr_flowframe[aggr_flowframe$Location == "LN",]

##################### RUN FUNCTION ##################################

genes <- colnames(aggr_flowframe[,2:39])
genes

###################### Intra TLS ####################################

Intra_genes <- genes[-c(22,23,28,29)]
# genes excluded from analysis:
genes[c(22,23,28,29)]

IntraTLS_stats <- get_WC_stats(genes = Intra_genes, df = aggr_IntraTLS, group = 'cluster_annot', comparison = 'Response')

###################### Peri TLS ####################################

Peri_genes <- genes[-c(22,23,25,28,29,36,37,38)]
# genes excluded from analysis:
genes[c(22,23,25,28,29,36,37,38)]

aggr_PeriTLS <- aggr_PeriTLS[!(aggr_PeriTLS$cluster_annot == "B_GC_prolif"),]
PeriTLS_stats <- get_WC_stats(genes = Peri_genes, df = aggr_PeriTLS, group = 'cluster_annot', comparison = 'Response')

###################### all TLS ####################################

allTLS_genes <- genes[-c(22,23,28,29)]
# genes excluded from analysis:
genes[c(22,23,28,29)]

allTLS_stats <- get_WC_stats(genes = allTLS_genes, df = aggr_allTLS, group = 'cluster_annot', comparison = 'Response')

###################### LN ####################################

LN_stats <- get_WC_stats(genes = genes, df = aggr_LN, group = 'cluster_annot', comparison = 'Response')

################# Extract DFs of Results ###########################

IntraTLS_stats_df <- do.call("rbind", IntraTLS_stats)
PeriTLS_stats_df <- do.call("rbind", PeriTLS_stats)
allTLS_stats_df <- do.call("rbind", allTLS_stats)
LN_stats_df <- do.call("rbind", LN_stats)

##### Save DFs of Outputs #####

write.xlsx(IntraTLS_stats_df, file="../figPlots/DEA/IMC_DEA_WC_stats_Padj.xlsx", sheetName="IntraTLS_Stats", append = TRUE, row.names=FALSE, col.names = TRUE)
write.xlsx(PeriTLS_stats_df, file="../figPlots/DEA/IMC_DEA_WC_stats_Padj.xlsx", sheetName="PeriTLS_Stats", append = TRUE, row.names=FALSE, col.names = TRUE)
write.xlsx(allTLS_stats_df, file="../figPlots/DEA/IMC_DEA_WC_stats_Padj.xlsx", sheetName="allTLS_Stats", append = TRUE, row.names=FALSE, col.names = TRUE)
write.xlsx(LN_stats_df, file="../figPlots/DEA/IMC_DEA_WC_stats_Padj.xlsx", sheetName="LN_Stats", append = TRUE, row.names=FALSE, col.names = TRUE)

