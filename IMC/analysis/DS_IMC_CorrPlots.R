#### Correlation Plots
# Meredith Wetzel
# August 13, 20204

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

aggr_flowframe <- readRDS("aggr_flowframe_withClusterCelltypes_Final.RDS")

#################### BREAK UP DATA INTO SUBGROUPINGS #########################

aggr_IntraTLS <- aggr_flowframe[aggr_flowframe$Location == "Intra",]
aggr_PeriTLS <- aggr_flowframe[aggr_flowframe$Location == "Peri",]
aggr_allTLS <- aggr_flowframe[aggr_flowframe$Location == "Intra" | aggr_flowframe$Location == "Peri",]
aggr_LN <- aggr_flowframe[aggr_flowframe$Location == "LN",]

aggr_IntraTLS_R <- aggr_IntraTLS[aggr_IntraTLS$Response == "R",]
aggr_IntraTLS_NR <- aggr_IntraTLS[aggr_IntraTLS$Response == "NR",]

aggr_PeriTLS_R <- aggr_PeriTLS[aggr_PeriTLS$Response == "R",]
aggr_PeriTLS_NR <- aggr_PeriTLS[aggr_PeriTLS$Response == "NR",]

aggr_allTLS_R <- aggr_allTLS[aggr_allTLS$Response == "R",]
aggr_allTLS_NR <- aggr_allTLS[aggr_allTLS$Response == "NR",]

aggr_LN_R <- aggr_LN[aggr_LN$Response == "R",]
aggr_LN_NR <- aggr_LN[aggr_LN$Response == "NR",]

#################### Correlations #########################

# using default of Pearson for correlation calc

aggr_IntraTLS_R_corr <- cor(aggr_IntraTLS_R[2:39], use = 'complete.obs') #subsetting the flowframe to just the numeric gene expression cols
IntraTLS_R_corrStats <- as.data.frame(aggr_IntraTLS_R_corr)
aggr_IntraTLS_NR_corr <- cor(aggr_IntraTLS_NR[2:39], use = 'complete.obs')
IntraTLS_NR_corrStats <- as.data.frame(aggr_IntraTLS_NR_corr)

aggr_PeriTLS_R_corr <- cor(aggr_PeriTLS_R[2:39], use = 'complete.obs')
PeriTLS_R_corrStats <- as.data.frame(aggr_PeriTLS_R_corr)
aggr_PeriTLS_NR_corr <- cor(aggr_PeriTLS_NR[2:39], use = 'complete.obs')
PeriTLS_NR_corrStats <- as.data.frame(aggr_PeriTLS_NR_corr)

aggr_allTLS_R_corr <- cor(aggr_allTLS_R[2:39], use = 'complete.obs')
allTLS_R_corrStats <- as.data.frame(aggr_allTLS_R_corr)
aggr_allTLS_NR_corr <- cor(aggr_allTLS_NR[2:39], use = 'complete.obs')
allTLS_NR_corrStats <- as.data.frame(aggr_allTLS_NR_corr)

aggr_LN_R_corr <- cor(aggr_LN_R[2:39], use = 'complete.obs')
LN_R_corrStats <- as.data.frame(aggr_LN_R_corr)
aggr_LN_NR_corr <- cor(aggr_LN_NR[2:39], use = 'complete.obs')
LN_NR_corrStats <- as.data.frame(aggr_LN_NR_corr)

# x = NA

pdf("../figPlots/correlations/IMC_correlation_plots.pdf", height = 8, width = 8)
corrplot(aggr_IntraTLS_R_corr, method = 'circle', order = 'alphabet', title = "Intra TLS Responder", na.label = "x",mar=c(0,0,2,0))
corrplot(aggr_IntraTLS_NR_corr, method = 'circle', order = 'alphabet', title = "Intra TLS Non-Responder", na.label = "x",mar=c(0,0,2,0))
corrplot(aggr_PeriTLS_R_corr, method = 'circle', order = 'alphabet', title = "Peri TLS Responder", na.label = "x",mar=c(0,0,2,0))
corrplot(aggr_PeriTLS_NR_corr, method = 'circle', order = 'alphabet', title = "Peri TLS Non-Responder", na.label = "x",mar=c(0,0,2,0))
corrplot(aggr_allTLS_R_corr, method = 'circle', order = 'alphabet', title = "Intra + Peri TLS Responder", na.label = "x",mar=c(0,0,2,0))
corrplot(aggr_allTLS_NR_corr, method = 'circle', order = 'alphabet', title = "Intra + Peri TLS Non-Responder", na.label = "x",mar=c(0,0,2,0))
corrplot(aggr_LN_R_corr, method = 'circle', order = 'alphabet', title = "LN Responder", na.label = "x",mar=c(0,0,2,0))
corrplot(aggr_LN_NR_corr, method = 'circle', order = 'alphabet', title = "LN Non-Responder", na.label = "x",mar=c(0,0,2,0))
dev.off()

#################### Extract P values and Corr Stats #########################

IntraTLS_R_testRes <- cor.mtest(aggr_IntraTLS_R[2:39])
IntraTLS_R_pvals <- as.data.frame(IntraTLS_R_testRes$p)
IntraTLS_NR_testRes <- cor.mtest(aggr_IntraTLS_NR[2:39])
IntraTLS_NR_pvals <- as.data.frame(IntraTLS_NR_testRes$p)

PeriTLS_R_testRes <- cor.mtest(aggr_PeriTLS_R[2:39])
PeriTLS_R_pvals <- as.data.frame(PeriTLS_R_testRes$p)
PeriTLS_NR_testRes <- cor.mtest(aggr_PeriTLS_NR[2:39])
PeriTLS_NR_pvals <- as.data.frame(PeriTLS_NR_testRes$p)

allTLS_R_testRes <- cor.mtest(aggr_allTLS_R[2:39])
allTLS_R_pvals <- as.data.frame(allTLS_R_testRes$p)
allTLS_NR_testRes <- cor.mtest(aggr_allTLS_NR[2:39])
allTLS_NR_pvals <- as.data.frame(allTLS_NR_testRes$p)

LN_R_testRes <- cor.mtest(aggr_LN_R[2:39])
LN_R_pvals <- as.data.frame(LN_R_testRes$p)
LN_NR_testRes <- cor.mtest(aggr_LN_NR[2:39])
LN_NR_pvals <- as.data.frame(LN_NR_testRes$p)

#################### Write Tables #########################
# each sheet in the excel file will be a df
#first add correlation stats
write.xlsx(IntraTLS_R_corrStats, file="../figPlots/correlations/IMC_marker_correlation_stats.xlsx", sheetName="IntraTLS_Responder_CorrStats", row.names=TRUE, col.names = TRUE)
write.xlsx(IntraTLS_NR_corrStats, file="../figPlots/correlations/IMC_marker_correlation_stats.xlsx", sheetName="IntraTLS_NonResponder_CorrStats", append = TRUE, row.names=TRUE, col.names = TRUE)
write.xlsx(PeriTLS_R_corrStats, file="../figPlots/correlations/IMC_marker_correlation_stats.xlsx", sheetName="PeriTLS_Responder_CorrStats", append = TRUE, row.names=TRUE, col.names = TRUE)
write.xlsx(PeriTLS_NR_corrStats, file="../figPlots/correlations/IMC_marker_correlation_stats.xlsx", sheetName="PeriTLS_NonResponder_CorrStats", append = TRUE, row.names=TRUE, col.names = TRUE)
write.xlsx(allTLS_R_corrStats, file="../figPlots/correlations/IMC_marker_correlation_stats.xlsx", sheetName="allTLS_Responder_CorrStats", append = TRUE, row.names=TRUE, col.names = TRUE)
write.xlsx(allTLS_NR_corrStats, file="../figPlots/correlations/IMC_marker_correlation_stats.xlsx", sheetName="allTLS_NonResponder_CorrStats", append = TRUE, row.names=TRUE, col.names = TRUE)
write.xlsx(LN_R_corrStats, file="../figPlots/correlations/IMC_marker_correlation_stats.xlsx", sheetName="LN_Responder_CorrStats", append = TRUE, row.names=TRUE, col.names = TRUE)
write.xlsx(LN_NR_corrStats, file="../figPlots/correlations/IMC_marker_correlation_stats.xlsx", sheetName="LN_NonResponder_CorrStats", append = TRUE, row.names=TRUE, col.names = TRUE)
#add in sheets for pvalues as well
write.xlsx(IntraTLS_R_pvals, file="../figPlots/correlations/IMC_marker_correlation_stats.xlsx", sheetName="IntraTLS_Responder_Pvals", append = TRUE, row.names=TRUE, col.names = TRUE)
write.xlsx(IntraTLS_NR_pvals, file="../figPlots/correlations/IMC_marker_correlation_stats.xlsx", sheetName="IntraTLS_NonResponder_Pvals", append = TRUE, row.names=TRUE, col.names = TRUE)
write.xlsx(PeriTLS_R_pvals, file="../figPlots/correlations/IMC_marker_correlation_stats.xlsx", sheetName="PeriTLS_Responder_Pvals", append = TRUE, row.names=TRUE, col.names = TRUE)
write.xlsx(PeriTLS_NR_pvals, file="../figPlots/correlations/IMC_marker_correlation_stats.xlsx", sheetName="PeriTLS_NonResponder_Pvals", append = TRUE, row.names=TRUE, col.names = TRUE)
write.xlsx(allTLS_R_pvals, file="../figPlots/correlations/IMC_marker_correlation_stats.xlsx", sheetName="allTLS_Responder_Pvals", append = TRUE, row.names=TRUE, col.names = TRUE)
write.xlsx(allTLS_NR_pvals, file="../figPlots/correlations/IMC_marker_correlation_stats.xlsx", sheetName="allTLS_NonResponder_Pvals", append = TRUE, row.names=TRUE, col.names = TRUE)
write.xlsx(LN_R_pvals, file="../figPlots/correlations/IMC_marker_correlation_stats.xlsx", sheetName="LN_Responder_Pvals", append = TRUE, row.names=TRUE, col.names = TRUE)
write.xlsx(LN_NR_pvals, file="../figPlots/correlations/IMC_marker_correlation_stats.xlsx", sheetName="LN_NonResponder_Pvals", append = TRUE, row.names=TRUE, col.names = TRUE)
