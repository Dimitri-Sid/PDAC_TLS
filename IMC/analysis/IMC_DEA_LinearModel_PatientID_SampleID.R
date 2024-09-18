## DEA with Daniel's Linear Model - IMC Data (Final Annotations)
# 1.PatientID only 
# 2.sample_id only 
# 3.PatientID + sample_id <<<<<<<<<<<
# 4.without PatientID or sample_id

# Meredith Wetzel
# August 16, 2024

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

#################### BREAK UP DATA INTO SUBGROUPINGS #########################

aggr_IntraTLS <- aggr_flowframe[aggr_flowframe$Location == "Intra",]
aggr_PeriTLS <- aggr_flowframe[aggr_flowframe$Location == "Peri",]
aggr_allTLS <- aggr_flowframe[aggr_flowframe$Location == "Intra" | aggr_flowframe$Location == "Peri",]
aggr_LN <- aggr_flowframe[aggr_flowframe$Location == "LN",]

########################## READ IN FUNCTION #########################

# function to run model on all genes
fit_gene_models <- function(genes, data) {
  res <- list()
  # Loop through each gene
  for (gene in genes) {
    formula <- as.formula(paste(gene, "~ Response * cluster_annot + (1|PatientID) + (1|sample_id)"))
    model <- lmer(formula, data = data) 
    res[[gene]] <- model #(just model)
    #res[[gene]] <- summary(model) #save to list obj (summaries)
  }
  return(res)  # Returns the list of models/summaries
} 

# extract genes
genes <- colnames(aggr_flowframe[,2:39])
genes

################### Run Model on FULL DATASET ################################

full_res <- fit_gene_models(genes = genes, data = aggr_flowframe)

################### Run Model on Intra TLS ################################

#IntraTLS_res <- fit_gene_models(genes = genes, data = aggr_IntraTLS)

# genes taken out of the analysis: 22,23,28,29
genes[c(22,23,28,29)]

IntraTLS_res_1to21 <- fit_gene_models(genes = genes[1:21], data = aggr_IntraTLS)
IntraTLS_res_24to27 <- fit_gene_models(genes = genes[24:27], data = aggr_IntraTLS)
IntraTLS_res_30to38 <- fit_gene_models(genes = genes[30:38], data = aggr_IntraTLS)
IntraTLS_res <- c(IntraTLS_res_1to21,IntraTLS_res_24to27,IntraTLS_res_30to38)

################### Run Model on Peri TLS ################################

#PeriTLS_res <- fit_gene_models(genes = genes, data = aggr_PeriTLS)

#genes taken out of the analysis: 22,23,28,29,36,37,38
genes[c(22,23,28,29,36,37,38)]

PeriTLS_res_1to21 <- fit_gene_models(genes = genes[1:21], data = aggr_PeriTLS)
PeriTLS_res_24to27 <- fit_gene_models(genes = genes[24:27], data = aggr_PeriTLS)
PeriTLS_res_30to35 <- fit_gene_models(genes = genes[30:35], data = aggr_PeriTLS)
PeriTLS_res <- c(PeriTLS_res_1to21,PeriTLS_res_24to27,PeriTLS_res_30to35)

################### Run Model on all TLS ################################

#allTLS_res <- fit_gene_models(genes = genes, data = aggr_allTLS)

# genes taken out of the analysis: 22,23,28,29
genes[c(22,23,28,29)]

allTLS_res_1to21 <- fit_gene_models(genes = genes[1:21], data = aggr_allTLS)
allTLS_res_24to27 <- fit_gene_models(genes = genes[24:27], data = aggr_allTLS)
allTLS_res_30to38 <- fit_gene_models(genes = genes[30:38], data = aggr_allTLS)
allTLS_res <- c(allTLS_res_1to21,allTLS_res_24to27,allTLS_res_30to38)

################### Run Model on LN ################################

LN_res <- fit_gene_models(genes = genes, data = aggr_LN)

################## Extract P Values ####################################

# define function to loop through model function results to make a list of df

get_model_coef <- function(model_res) {
  res <- list()
  res_genes <- names(model_res)
  for (gene in res_genes) {
    df <- as.data.frame(coef(summary(model_res[[gene]])))
    res[[gene]] <- df
  }
  return(res)  # Return the list of dfs/pvalues
} 

################## Run Function to Create DFs ####################################

##### Whole Dataset #####

full_res_coef <- get_model_coef(full_res)
# merge all coef df together into one
full_res_coef_df <- do.call("rbind", full_res_coef)
# use rownames to create columns for the gene/expression and comparison from the model
res <- str_split_fixed(rownames(full_res_coef_df), pattern = "\\.", n = 2)
colnames(res) <- c("Gene","Comparison")
# merge string columns to coef df for ease of filtering
full_res_coef_df <- cbind(full_res_coef_df,res)
# split strings again to extract cluster annotations - only for the ones of interest (response&cluster)
res2 <- str_split_fixed(full_res_coef_df$Comparison, pattern = "ResponseR:cluster_annot", n = 2)
colnames(res2) <- c("Extra","cluster_annot")
full_res_coef_df <- cbind(full_res_coef_df,res2)
full_res_coef_df <- full_res_coef_df[,!(names(full_res_coef_df) %in% "Extra")] 

##### Intra TLS Subset #####

IntraTLS_res_coef <- get_model_coef(IntraTLS_res)
IntraTLS_res_coef_df <- do.call("rbind", IntraTLS_res_coef)
res <- str_split_fixed(rownames(IntraTLS_res_coef_df), pattern = "\\.", n = 2)
colnames(res) <- c("Gene","Comparison")
IntraTLS_res_coef_df <- cbind(IntraTLS_res_coef_df,res)
IntraTLS_res_coef_df$Location <- "Intra"
res2 <- str_split_fixed(IntraTLS_res_coef_df$Comparison, pattern = "ResponseR:cluster_annot", n = 2)
colnames(res2) <- c("Extra","cluster_annot")
IntraTLS_res_coef_df <- cbind(IntraTLS_res_coef_df,res2)
IntraTLS_res_coef_df <- IntraTLS_res_coef_df[,!(names(IntraTLS_res_coef_df) %in% "Extra")]

##### Peri TLS Subset #####

PeriTLS_res_coef <- get_model_coef(PeriTLS_res)
PeriTLS_res_coef_df <- do.call("rbind", PeriTLS_res_coef)
res <- str_split_fixed(rownames(PeriTLS_res_coef_df), pattern = "\\.", n = 2)
colnames(res) <- c("Gene","Comparison")
PeriTLS_res_coef_df <- cbind(PeriTLS_res_coef_df,res)
PeriTLS_res_coef_df$Location <- "Peri"
res2 <- str_split_fixed(PeriTLS_res_coef_df$Comparison, pattern = "ResponseR:cluster_annot", n = 2)
colnames(res2) <- c("Extra","cluster_annot")
PeriTLS_res_coef_df <- cbind(PeriTLS_res_coef_df,res2)
PeriTLS_res_coef_df <- PeriTLS_res_coef_df[,!(names(PeriTLS_res_coef_df) %in% "Extra")] 

##### All TLS Subset #####

allTLS_res_coef <- get_model_coef(allTLS_res)
allTLS_res_coef_df <- do.call("rbind", allTLS_res_coef)
res <- str_split_fixed(rownames(allTLS_res_coef_df), pattern = "\\.", n = 2)
colnames(res) <- c("Gene","Comparison")
allTLS_res_coef_df <- cbind(allTLS_res_coef_df,res)
res2 <- str_split_fixed(allTLS_res_coef_df$Comparison, pattern = "ResponseR:cluster_annot", n = 2)
colnames(res2) <- c("Extra","cluster_annot")
allTLS_res_coef_df <- cbind(allTLS_res_coef_df,res2)
allTLS_res_coef_df <- allTLS_res_coef_df[,!(names(allTLS_res_coef_df) %in% "Extra")] 

##### LN Subset #####

LN_res_coef <- get_model_coef(LN_res)
LN_res_coef_df <- do.call("rbind", LN_res_coef)
res <- str_split_fixed(rownames(LN_res_coef_df), pattern = "\\.", n = 2)
colnames(res) <- c("Gene","Comparison")
LN_res_coef_df <- cbind(LN_res_coef_df,res)
LN_res_coef_df$Location <- "LN"
res2 <- str_split_fixed(LN_res_coef_df$Comparison, pattern = "ResponseR:cluster_annot", n = 2)
colnames(res2) <- c("Extra","cluster_annot")
LN_res_coef_df <- cbind(LN_res_coef_df,res2)
LN_res_coef_df <- LN_res_coef_df[,!(names(LN_res_coef_df) %in% "Extra")] 

##### Save DFs of Outputs #####

write.xlsx(full_res_coef_df, file="../figPlots/DEA/IMC_DEA_LMER_stats_PatientID_SampleID.xlsx", sheetName="WholeDataset_Stats", row.names=TRUE, col.names = TRUE)
write.xlsx(IntraTLS_res_coef_df, file="../figPlots/DEA/IMC_DEA_LMER_stats_PatientID_SampleID.xlsx", sheetName="IntraTLS_Stats", append = TRUE, row.names=TRUE, col.names = TRUE)
write.xlsx(PeriTLS_res_coef_df, file="../figPlots/DEA/IMC_DEA_LMER_stats_PatientID_SampleID.xlsx", sheetName="PeriTLS_Stats", append = TRUE, row.names=TRUE, col.names = TRUE)
write.xlsx(allTLS_res_coef_df, file="../figPlots/DEA/IMC_DEA_LMER_stats_PatientID_SampleID.xlsx", sheetName="allTLS_Stats", append = TRUE, row.names=TRUE, col.names = TRUE)
write.xlsx(LN_res_coef_df, file="../figPlots/DEA/IMC_DEA_LMER_stats_PatientID_SampleID.xlsx", sheetName="LN_Stats", append = TRUE, row.names=TRUE, col.names = TRUE)

table(aggr_flowframe$cluster_annot, aggr_flowframe$Response)

# show which cell types were dropped in each group:

unique(aggr_flowframe$cluster_annot)[!(unique(aggr_flowframe$cluster_annot) %in% unique(full_res_coef_df$cluster_annot))]
unique(aggr_flowframe$cluster_annot)[!(unique(aggr_flowframe$cluster_annot) %in% unique(IntraTLS_res_coef_df$cluster_annot))]
unique(aggr_flowframe$cluster_annot)[!(unique(aggr_flowframe$cluster_annot) %in% unique(PeriTLS_res_coef_df$cluster_annot))]
unique(aggr_flowframe$cluster_annot)[!(unique(aggr_flowframe$cluster_annot) %in% unique(allTLS_res_coef_df$cluster_annot))]
unique(aggr_flowframe$cluster_annot)[!(unique(aggr_flowframe$cluster_annot) %in% unique(LN_res_coef_df$cluster_annot))]
