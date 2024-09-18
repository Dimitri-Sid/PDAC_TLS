## Composition Plot Statistics - whole dataset

# Meredith Wetzel - Modified from Dimitri
# Updated: August 20, 2024

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

#fully annotated/gated object from DS_IMCgating_plots script
aggr_flowframe <- readRDS("aggr_flowframe_withClusterCelltypes_Final_Updated.RDS")

## Read in individual prop data frames ##

whole_ggdf <- readRDS("../figPlots/composition/Whole_props_ggdf.RDS")
IntraTLS_ggdf <- readRDS("../figPlots/composition/IntraTLS_props_ggdf.RDS")
PeriTLS_ggdf <- readRDS("../figPlots/composition/PeriTLS_props_ggdf.RDS")
allTLS_ggdf <- readRDS("../figPlots/composition/allTLS_props_ggdf.RDS")
LN_ggdf <- readRDS("../figPlots/composition/LN_props_ggdf.RDS")

clusterlevels <- unique(aggr_flowframe$cluster_annot)

####ABUNDANCE DIFFERENTIAL ANALYSES####
#for comparing Response (% cells)

### WHOLE ###

# anova

pvalueprops<-data.frame(type=NA,value=NA) #first for AOV

for(i in 1:length(clusterlevels)){
  ggdfsubset<-whole_ggdf[whole_ggdf$cluster==clusterlevels[i],]
  
  res<-aov(proportion ~ Response, data = ggdfsubset)
  res<-summary(res)
  res<-unlist(res)["Pr(>F)1"]
  pvalueprops<-rbind(pvalueprops,data.frame(type=clusterlevels[i],value=unlist(res)["Pr(>F)1"]))
  pvalueprops<-pvalueprops[!is.na(pvalueprops$value),]
}
rownames(pvalueprops)<-1:nrow(pvalueprops)
pvalueprops$FDRadjP<-p.adjust(pvalueprops$value, method="BH") #BH = FDR
write.csv(pvalueprops,"../figPlots/composition/Whole_props_Response_stats_aov.csv")

# Ttest

pvalueprops<-data.frame(type=NA,Var1=NA,Var2=NA,value=NA)

for(i in 1:length(clusterlevels)){
  ggdfsubset<-whole_ggdf[whole_ggdf$cluster==clusterlevels[i],]
  
  res<-pairwise.t.test(ggdfsubset$proportion, ggdfsubset$Response, paired = F, p.adjust.method = "none")
  pvalueprops<-rbind(pvalueprops,cbind(type=rep(clusterlevels[i],nrow(melt(res$p.value))),melt(res$p.value)))
  pvalueprops<-pvalueprops[!is.na(pvalueprops$value),]
}
rownames(pvalueprops)<-1:nrow(pvalueprops)
pvalueprops$FDRadjP<-p.adjust(pvalueprops$value, method="BH")
write.csv(pvalueprops,"../figPlots/composition/Whole_props_Response_stats_Ttest.csv")

### IntraTLS ###

# anova

pvalueprops<-data.frame(type=NA,value=NA) #first for AOV

for(i in 1:length(clusterlevels)){
  ggdfsubset<-IntraTLS_ggdf[IntraTLS_ggdf$cluster==clusterlevels[i],]
  
  res<-aov(proportion ~ Response, data = ggdfsubset)
  res<-summary(res)
  res<-unlist(res)["Pr(>F)1"]
  pvalueprops<-rbind(pvalueprops,data.frame(type=clusterlevels[i],value=unlist(res)["Pr(>F)1"]))
  pvalueprops<-pvalueprops[!is.na(pvalueprops$value),]
}
rownames(pvalueprops)<-1:nrow(pvalueprops)
pvalueprops$FDRadjP<-p.adjust(pvalueprops$value, method="BH") #BH = FDR
write.csv(pvalueprops,"../figPlots/composition/IntraTLS_props_Response_stats_aov.csv")

# Ttest

pvalueprops<-data.frame(type=NA,Var1=NA,Var2=NA,value=NA)

for(i in 1:length(clusterlevels)){
  ggdfsubset<-IntraTLS_ggdf[IntraTLS_ggdf$cluster==clusterlevels[i],]
  
  res<-pairwise.t.test(ggdfsubset$proportion, ggdfsubset$Response, paired = F, p.adjust.method = "none")
  pvalueprops<-rbind(pvalueprops,cbind(type=rep(clusterlevels[i],nrow(melt(res$p.value))),melt(res$p.value)))
  pvalueprops<-pvalueprops[!is.na(pvalueprops$value),]
}
rownames(pvalueprops)<-1:nrow(pvalueprops)
pvalueprops$FDRadjP<-p.adjust(pvalueprops$value, method="BH")
write.csv(pvalueprops,"../figPlots/composition/IntraTLS_props_Response_stats_Ttest.csv")

### PeriTLS ###

# anova

pvalueprops<-data.frame(type=NA,value=NA) #first for AOV

for(i in 1:length(clusterlevels)){
  ggdfsubset<-PeriTLS_ggdf[PeriTLS_ggdf$cluster==clusterlevels[i],]
  
  res<-aov(proportion ~ Response, data = ggdfsubset)
  res<-summary(res)
  res<-unlist(res)["Pr(>F)1"]
  pvalueprops<-rbind(pvalueprops,data.frame(type=clusterlevels[i],value=unlist(res)["Pr(>F)1"]))
  pvalueprops<-pvalueprops[!is.na(pvalueprops$value),]
}
rownames(pvalueprops)<-1:nrow(pvalueprops)
pvalueprops$FDRadjP<-p.adjust(pvalueprops$value, method="BH") #BH = FDR
write.csv(pvalueprops,"../figPlots/composition/PeriTLS_props_Response_stats_aov.csv")

# Ttest

pvalueprops<-data.frame(type=NA,Var1=NA,Var2=NA,value=NA)

for(i in 1:length(clusterlevels)){
  ggdfsubset<-PeriTLS_ggdf[PeriTLS_ggdf$cluster==clusterlevels[i],]
  
  res<-pairwise.t.test(ggdfsubset$proportion, ggdfsubset$Response, paired = F, p.adjust.method = "none")
  pvalueprops<-rbind(pvalueprops,cbind(type=rep(clusterlevels[i],nrow(melt(res$p.value))),melt(res$p.value)))
  pvalueprops<-pvalueprops[!is.na(pvalueprops$value),]
}
rownames(pvalueprops)<-1:nrow(pvalueprops)
pvalueprops$FDRadjP<-p.adjust(pvalueprops$value, method="BH")
write.csv(pvalueprops,"../figPlots/composition/PeriTLS_props_Response_stats_Ttest.csv")

### allTLS ###

# anova

pvalueprops<-data.frame(type=NA,value=NA) #first for AOV

for(i in 1:length(clusterlevels)){
  ggdfsubset<-allTLS_ggdf[allTLS_ggdf$cluster==clusterlevels[i],]
  
  res<-aov(proportion ~ Response, data = ggdfsubset)
  res<-summary(res)
  res<-unlist(res)["Pr(>F)1"]
  pvalueprops<-rbind(pvalueprops,data.frame(type=clusterlevels[i],value=unlist(res)["Pr(>F)1"]))
  pvalueprops<-pvalueprops[!is.na(pvalueprops$value),]
}
rownames(pvalueprops)<-1:nrow(pvalueprops)
pvalueprops$FDRadjP<-p.adjust(pvalueprops$value, method="BH") #BH = FDR
write.csv(pvalueprops,"../figPlots/composition/allTLS_props_Response_stats_aov.csv")

# Ttest

pvalueprops<-data.frame(type=NA,Var1=NA,Var2=NA,value=NA)

for(i in 1:length(clusterlevels)){
  ggdfsubset<-allTLS_ggdf[allTLS_ggdf$cluster==clusterlevels[i],]
  
  res<-pairwise.t.test(ggdfsubset$proportion, ggdfsubset$Response, paired = F, p.adjust.method = "none")
  pvalueprops<-rbind(pvalueprops,cbind(type=rep(clusterlevels[i],nrow(melt(res$p.value))),melt(res$p.value)))
  pvalueprops<-pvalueprops[!is.na(pvalueprops$value),]
}
rownames(pvalueprops)<-1:nrow(pvalueprops)
pvalueprops$FDRadjP<-p.adjust(pvalueprops$value, method="BH")
write.csv(pvalueprops,"../figPlots/composition/allTLS_props_Response_stats_Ttest.csv")

### LN ###

# anova

pvalueprops<-data.frame(type=NA,value=NA) #first for AOV

for(i in 1:length(clusterlevels)){
  ggdfsubset<-LN_ggdf[LN_ggdf$cluster==clusterlevels[i],]
  
  res<-aov(proportion ~ Response, data = ggdfsubset)
  res<-summary(res)
  res<-unlist(res)["Pr(>F)1"]
  pvalueprops<-rbind(pvalueprops,data.frame(type=clusterlevels[i],value=unlist(res)["Pr(>F)1"]))
  pvalueprops<-pvalueprops[!is.na(pvalueprops$value),]
}
rownames(pvalueprops)<-1:nrow(pvalueprops)
pvalueprops$FDRadjP<-p.adjust(pvalueprops$value, method="BH") #BH = FDR
write.csv(pvalueprops,"../figPlots/composition/LN_props_Response_stats_aov.csv")

# Ttest

pvalueprops<-data.frame(type=NA,Var1=NA,Var2=NA,value=NA)

for(i in 1:length(clusterlevels)){
  ggdfsubset<-LN_ggdf[LN_ggdf$cluster==clusterlevels[i],]
  
  res<-pairwise.t.test(ggdfsubset$proportion, ggdfsubset$Response, paired = F, p.adjust.method = "none")
  pvalueprops<-rbind(pvalueprops,cbind(type=rep(clusterlevels[i],nrow(melt(res$p.value))),melt(res$p.value)))
  #pvalueprops<-pvalueprops[!is.na(pvalueprops$value),]
}
rownames(pvalueprops)<-1:nrow(pvalueprops)
#pvalueprops$FDRadjP<-p.adjust(pvalueprops$value, method="BH")
write.csv(pvalueprops,"../figPlots/composition/LN_props_Response_stats_Ttest.csv")
