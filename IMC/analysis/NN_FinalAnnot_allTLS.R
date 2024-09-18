
# Nearest Neighbor Abundance with Final Annotations
# 1.Intra TLS R vs. NR 
# 2.Peri TLS R vs. NR
# 3.All TLS R. vs NR <<<<<<<<<<<
# 4.All LN R vs. NR

# Meredith Wetzel
# August 15, 2024

# Calculate the nearest neighbors
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
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

output<-readRDS(paste0(script_dir,'output_w15.RDS')) #saved output object with 15 clusters

#fully annotated object from "Add_NN_Annotation.R" script
aggr_flowframe <- readRDS("aggr_flowframe_withClusterCelltypes_Final_Updated_NNannotated.RDS")

# colors should be in order of their matched labels
neighbor_cols <- c("forestgreen","lightgreen","cyan","lightblue","orange","purple","blue","brown","red","pink")
names(neighbor_cols) <- unique(aggr_flowframe$cluster_annot)

add_missing_rows <- function(df, vec_rows, fill_value = '0') {
  missing_names <- setdiff(vec_rows, rownames(df))
  new_rows <- data.frame(matrix(fill_value, nrow = length(missing_names), ncol = ncol(df))) #make a df with the missing rows all with 0 values
  rownames(new_rows) <- missing_names
  colnames(new_rows) <- colnames(df)
  #new_rows <- new_rows[rownames(df),] #make the df the same order of rows
  df <- rbind(df, new_rows)
  df$Var2 <- rownames(df)
  return(df)
}

########################################### BREAK UP DATA INTO SUBGROUPINGS ##########################################

aggr_allTLS <- aggr_flowframe[aggr_flowframe$Location == "Intra" | aggr_flowframe$Location == "Peri",]

aggr_allTLS_R <- aggr_allTLS[aggr_allTLS$Response == "R",]
aggr_allTLS_NR <- aggr_allTLS[aggr_allTLS$Response == "NR",]

########################################### ALL TLS BY RESPONSE FREQ ##########################################

# Create tables of counts - remember that the rows here will be the clusters and the columns will be the neighbor type

table(aggr_allTLS_R$cluster_annot, aggr_allTLS_R$N1type)

NN_1_df <- as.data.frame(table(aggr_allTLS_R$cluster_annot, aggr_allTLS_R$N1type))
names(NN_1_df) #Var1 = clusters; Var2 = neighbor type
NN_1_df <- as.data.frame(pivot_wider(NN_1_df, names_from = "Var1", values_from = "Freq"))
NN_1_df
rownames(NN_1_df) <- NN_1_df$Var2
N1_df_small <- NN_1_df[,-1] 
props_table_N1 <- (N1_df_small / rowSums(N1_df_small)) * 100 
rownames(props_table_N1) <- NN_1_df$Var2
ggdf_N1 <- melt(data.frame(cluster = rownames(props_table_N1),props_table_N1, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")

pdf("../figPlots/NN/allTLS/N1_neighbor_props_allTLS_R.pdf",height = 10, width = 10) 
ggplot(ggdf_N1, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  # coord_polar(theta = "y") +
  # facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Responder N1 Neighbors")
dev.off()


NN_2_df <- as.data.frame(table(aggr_allTLS_R$cluster_annot, aggr_allTLS_R$N2type))
names(NN_2_df)
NN_2_df <- as.data.frame(pivot_wider(NN_2_df, names_from = "Var1", values_from = "Freq"))
NN_2_df
rownames(NN_2_df) <- NN_2_df$Var2
props_table_N2 <- (NN_2_df[,-1] / rowSums(NN_2_df[,-1])) * 100 
rownames(props_table_N2) <- NN_2_df$Var2
ggdf_N2 <- melt(data.frame(cluster = rownames(props_table_N2),props_table_N2, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N2_neighbor_props_allTLS_R.pdf",height = 10, width = 10)
ggplot(ggdf_N2, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  # coord_polar(theta = "y") +
  # facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Responder N2 Neighbors")
dev.off()


NN_3_df <- as.data.frame(table(aggr_allTLS_R$cluster_annot, aggr_allTLS_R$N3type))
names(NN_3_df) #Var1 = clusters; Var2 = neighbor type 
NN_3_df <- as.data.frame(pivot_wider(NN_3_df, names_from = "Var1", values_from = "Freq"))
NN_3_df
rownames(NN_3_df) <- NN_3_df$Var2
props_table_N3 <- (NN_3_df[,-1] / rowSums(NN_3_df[,-1])) * 100 
rownames(props_table_N3) <- NN_3_df$Var2
ggdf_N3 <- melt(data.frame(cluster = rownames(props_table_N3),props_table_N3, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N3_neighbor_props_allTLS_R.pdf",height = 10, width = 10)
ggplot(ggdf_N3, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  # coord_polar(theta = "y") +
  # facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Responder N3 Neighbors")
dev.off()

cols = colnames(NN_1_df[,-1])

NN_1to3_df <- NN_1_df

NN_1to3_df[, cols] = NN_1_df[, cols] + NN_2_df[, cols] + NN_3_df[, cols]

NN_1to3_df_small <- as.data.frame(NN_1to3_df[,-1])

rownames(NN_1to3_df_small) <- NN_1to3_df$Var2

props_table_100 <- (NN_1to3_df_small / rowSums(NN_1to3_df_small)) * 100 

(NN_1to3_df_small[1,1] / sum(NN_1to3_df_small[1,])) * 100

#set up df for plotting

ggdf_1to3 <- melt(data.frame(cluster = rownames(props_table_100),props_table_100, check.names = FALSE),
                  id.vars = "cluster", value.name = "proportion", 
                  variable.name = "neighbor")

colnames(ggdf_1to3)

pdf("../figPlots/NN/allTLS/neighbor_props_NN1to3_barplot_allTLS_R.pdf",height = 10,width = 10)
ggplot(ggdf_1to3, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Proportion') +
  xlab('Cluster') +
  guides(fill=guide_legend(title="Neighbors"))+
  ggtitle("Responder N1-N3 Neighbors All TLS")
dev.off()

NN_4_df <- as.data.frame(table(aggr_allTLS_R$cluster_annot, aggr_allTLS_R$N4type))
NN_4_df <- as.data.frame(pivot_wider(NN_4_df, names_from = "Var1", values_from = "Freq"))
rownames(NN_4_df) <- NN_4_df$Var2
props_table_N4 <- (NN_4_df[,-1] / rowSums(NN_4_df[,-1])) * 100 
rownames(props_table_N4) <- NN_4_df$Var2
ggdf_N4 <- melt(data.frame(cluster = rownames(props_table_N4),props_table_N4, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N4_neighbor_props_allTLS_R.pdf",height = 10, width = 10)
ggplot(ggdf_N4, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  # coord_polar(theta = "y") +
  # facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Responder N4 Neighbors")
dev.off()


NN_5_df <- as.data.frame(table(aggr_allTLS_R$cluster_annot, aggr_allTLS_R$N5type))
NN_5_df <- as.data.frame(pivot_wider(NN_5_df, names_from = "Var1", values_from = "Freq"))
rownames(NN_5_df) <- NN_5_df$Var2
props_table_N5 <- (NN_5_df[,-1] / rowSums(NN_5_df[,-1])) * 100 
rownames(props_table_N5) <- NN_5_df$Var2
ggdf_N5 <- melt(data.frame(cluster = rownames(props_table_N5),props_table_N5, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N5_neighbor_props_allTLS_R.pdf",height = 10, width = 10)
ggplot(ggdf_N5, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  # coord_polar(theta = "y") +
  # facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Responder N5 Neighbors")
dev.off()


NN_6_df <- as.data.frame(table(aggr_allTLS_R$cluster_annot, aggr_allTLS_R$N6type))
NN_6_df <- as.data.frame(pivot_wider(NN_6_df, names_from = "Var1", values_from = "Freq"))
rownames(NN_6_df) <- NN_6_df$Var2
props_table_N6 <- (NN_6_df[,-1] / rowSums(NN_6_df[,-1])) * 100 
rownames(props_table_N6) <- NN_6_df$Var2
ggdf_N6 <- melt(data.frame(cluster = rownames(props_table_N6),props_table_N6, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N6_neighbor_props_allTLS_R.pdf",height = 10, width = 10)
ggplot(ggdf_N6, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  # coord_polar(theta = "y") +
  # facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Responder N6 Neighbors")
dev.off()


NN_7_df <- as.data.frame(table(aggr_allTLS_R$cluster_annot, aggr_allTLS_R$N7type))
NN_7_df <- as.data.frame(pivot_wider(NN_7_df, names_from = "Var1", values_from = "Freq"))
rownames(NN_7_df) <- NN_7_df$Var2
props_table_N7 <- (NN_7_df[,-1] / rowSums(NN_7_df[,-1])) * 100 
rownames(props_table_N7) <- NN_7_df$Var2
ggdf_N7 <- melt(data.frame(cluster = rownames(props_table_N7),props_table_N7, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N7_neighbor_props_allTLS_R.pdf",height = 10, width = 10)
ggplot(ggdf_N7, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  #coord_polar(theta = "y") +
  #facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Responder N7 Neighbors")
dev.off()


NN_8_df <- as.data.frame(table(aggr_allTLS_R$cluster_annot, aggr_allTLS_R$N8type))
NN_8_df <- as.data.frame(pivot_wider(NN_8_df, names_from = "Var1", values_from = "Freq"))
rownames(NN_8_df) <- NN_8_df$Var2
props_table_N8 <- (NN_8_df[,-1] / rowSums(NN_8_df[,-1])) * 100 
rownames(props_table_N8) <- NN_8_df$Var2
ggdf_N8 <- melt(data.frame(cluster = rownames(props_table_N8),props_table_N8, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N8_neighbor_props_allTLS_R.pdf",height = 10, width = 10)
ggplot(ggdf_N8, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  #coord_polar(theta = "y") +
  #facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Responder N8 Neighbors")
dev.off()


NN_9_df <- as.data.frame(table(aggr_allTLS_R$cluster_annot, aggr_allTLS_R$N9type))
NN_9_df <- as.data.frame(pivot_wider(NN_9_df, names_from = "Var1", values_from = "Freq"))
rownames(NN_9_df) <- NN_9_df$Var2
props_table_N9 <- (NN_9_df[,-1] / rowSums(NN_9_df[,-1])) * 100 
rownames(props_table_N9) <- NN_9_df$Var2
ggdf_N9 <- melt(data.frame(cluster = rownames(props_table_N9),props_table_N9, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N9_neighbor_props_allTLS_R.pdf",height = 10, width = 10)
ggplot(ggdf_N9, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  #coord_polar(theta = "y") +
  #facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Responder N9 Neighbors")
dev.off()


NN_10_df <- as.data.frame(table(aggr_allTLS_R$cluster_annot, aggr_allTLS_R$N10type))
NN_10_df <- as.data.frame(pivot_wider(NN_10_df, names_from = "Var1", values_from = "Freq"))
rownames(NN_10_df) <- NN_10_df$Var2
props_table_N10 <- (NN_10_df[,-1] / rowSums(NN_10_df[,-1])) * 100
rownames(props_table_N10) <- NN_10_df$Var2
ggdf_N10 <- melt(data.frame(cluster = rownames(props_table_N10),props_table_N10, check.names = FALSE),
                 id.vars = "cluster", value.name = "proportion", 
                 variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N10_neighbor_props_allTLS_R.pdf",height = 10, width = 10)
ggplot(ggdf_N10, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  #coord_polar(theta = "y") +
  #facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Responder N10 Neighbors")
dev.off()

cols = colnames(NN_1_df[,-1])

all_NN_df <- NN_1_df

all_NN_df[, cols] = NN_1_df[, cols] + NN_2_df[, cols] + NN_3_df[, cols] + NN_4_df[, cols] + 
  NN_5_df[, cols] + NN_6_df[, cols] + NN_7_df[, cols] + NN_8_df[, cols] +
  NN_9_df[, cols] + NN_10_df[, cols] #add all the values together

all_NN_df_small <- as.data.frame(all_NN_df[,-1])

rownames(all_NN_df_small) <- all_NN_df$Var2

props_table_100 <- (all_NN_df_small / rowSums(all_NN_df_small)) * 100
props_table <- (all_NN_df_small / rowSums(all_NN_df_small)) 

(all_NN_df_small[1,1] / sum(all_NN_df_small[1,])) * 100 

#set up df for plotting

ggdf <- melt(data.frame(cluster = rownames(props_table_100),props_table_100, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "neighbor")

saveRDS(ggdf, file = "../figPlots/NN/allTLS/NN_propsTable_ggdf.RDS")

colnames(ggdf)

pdf("../figPlots/NN/allTLS/neighbor_props_allTLS_R.pdf",height = 10,width = 10)
ggplot(ggdf, aes(x = "", y = proportion, fill = neighbor)) +
  geom_col() +
  coord_polar(theta = "y") +
  facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal()
dev.off()

pdf("../figPlots/NN/allTLS/neighbor_props_barplot_allTLS_R.pdf",height = 10,width = 10)
ggplot(ggdf, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Responder All 10 Neighbors Summed All TLS")
dev.off()

# Hole size
hsize <- 1.5

df <- ggdf %>% 
  mutate(x = hsize)

pdf("../figPlots/NN/allTLS/neighbor_props_donut_allTLS_R.pdf",height = 10,width = 10)
ggplot(df, aes(x = hsize, y = proportion, fill = neighbor)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(0.2, hsize + 0.5)) +
  theme_minimal()+
  scale_fill_manual(values = neighbor_cols) +
  facet_wrap(~cluster)
dev.off()


###### NON-RESPONDERS (ALL TLS) ###############################################

table(aggr_allTLS_NR$cluster_annot, aggr_allTLS_NR$N1type)

NN_1_df <- as.data.frame(table(aggr_allTLS_NR$cluster_annot, aggr_allTLS_NR$N1type))
names(NN_1_df) #Var1 = clusters; Var2 = neighbor type
NN_1_df <- as.data.frame(pivot_wider(NN_1_df, names_from = "Var1", values_from = "Freq"))
NN_1_df
rownames(NN_1_df) <- NN_1_df$Var2
N1_df_small <- NN_1_df[,-1] 
props_table_N1 <- (N1_df_small / rowSums(N1_df_small)) * 100 
rownames(props_table_N1) <- NN_1_df$Var2
ggdf_N1 <- melt(data.frame(cluster = rownames(props_table_N1),props_table_N1, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")

pdf("../figPlots/NN/allTLS/N1_neighbor_props_allTLS_NR.pdf",height = 10, width = 10) 
ggplot(ggdf_N1, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  # coord_polar(theta = "y") +
  # facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Non-Responder N1 Neighbors")
dev.off()


NN_2_df <- as.data.frame(table(aggr_allTLS_NR$cluster_annot, aggr_allTLS_NR$N2type))
names(NN_2_df)
NN_2_df <- as.data.frame(pivot_wider(NN_2_df, names_from = "Var1", values_from = "Freq"))
NN_2_df
rownames(NN_2_df) <- NN_2_df$Var2
props_table_N2 <- (NN_2_df[,-1] / rowSums(NN_2_df[,-1])) * 100
rownames(props_table_N2) <- NN_2_df$Var2
ggdf_N2 <- melt(data.frame(cluster = rownames(props_table_N2),props_table_N2, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N2_neighbor_props_allTLS_NR.pdf",height = 10, width = 10)
ggplot(ggdf_N2, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  # coord_polar(theta = "y") +
  # facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Non-Responder N2 Neighbors")
dev.off()


NN_3_df <- as.data.frame(table(aggr_allTLS_NR$cluster_annot, aggr_allTLS_NR$N3type))
names(NN_3_df) #Var1 = clusters; Var2 = neighbor type
NN_3_df <- as.data.frame(pivot_wider(NN_3_df, names_from = "Var1", values_from = "Freq"))
NN_3_df
rownames(NN_3_df) <- NN_3_df$Var2
props_table_N3 <- (NN_3_df[,-1] / rowSums(NN_3_df[,-1])) * 100 
rownames(props_table_N3) <- NN_3_df$Var2
ggdf_N3 <- melt(data.frame(cluster = rownames(props_table_N3),props_table_N3, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N3_neighbor_props_allTLS_NR.pdf",height = 10, width = 10)
ggplot(ggdf_N3, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  # coord_polar(theta = "y") +
  # facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Non-Responder N3 Neighbors")
dev.off()

cols = colnames(NN_1_df[,-1])

rows <- unique(aggr_allTLS_NR$cluster_annot) 

add_missing_rows <- function(df, vec_rows, fill_value = '0') {
  missing_names <- setdiff(vec_rows, rownames(df))
  new_rows <- data.frame(matrix(fill_value, nrow = length(missing_names), ncol = ncol(df))) #make a df with the missing rows all with 0 values
  rownames(new_rows) <- missing_names
  colnames(new_rows) <- colnames(df)
  #new_rows <- new_rows[rownames(df),] #make the df the same order of rows
  df <- rbind(df, new_rows)
  df$Var2 <- rownames(df)
  return(df)
}

NN_1_df <- add_missing_rows(df = NN_1_df, vec_rows = rows)
NN_2_df <- add_missing_rows(df = NN_2_df, vec_rows = rows)
NN_3_df <- add_missing_rows(df = NN_3_df, vec_rows = rows)

NN_1to3_df <- NN_1_df
NN_1to3_df[, cols] = as.numeric(unlist(NN_1_df[, cols])) + as.numeric(unlist(NN_2_df[, cols])) + as.numeric(unlist(NN_3_df[, cols]))
NN_1to3_df_small <- as.data.frame(NN_1to3_df[,-1])
rownames(NN_1to3_df_small) <- NN_1to3_df$Var2
props_table_100 <- (NN_1to3_df_small / rowSums(NN_1to3_df_small)) * 100
(NN_1to3_df_small[1,1] / sum(NN_1to3_df_small[1,])) * 100

#set up df for plotting

ggdf_1to3 <- melt(data.frame(cluster = rownames(props_table_100),props_table_100, check.names = FALSE),
                  id.vars = "cluster", value.name = "proportion", 
                  variable.name = "neighbor")

colnames(ggdf_1to3)

pdf("../figPlots/NN/allTLS/neighbor_props_NN1to3_barplot_allTLS_NR.pdf",height = 10,width = 10)
ggplot(ggdf_1to3, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Proportion') +
  xlab('Cluster') +
  guides(fill=guide_legend(title="Neighbors"))+
  ggtitle("Non-Responder N1-N3 Neighbors All TLS")
dev.off()

NN_4_df <- as.data.frame(table(aggr_allTLS_NR$cluster_annot, aggr_allTLS_NR$N4type))
NN_4_df <- as.data.frame(pivot_wider(NN_4_df, names_from = "Var1", values_from = "Freq"))
rownames(NN_4_df) <- NN_4_df$Var2
props_table_N4 <- (NN_4_df[,-1] / rowSums(NN_4_df[,-1])) * 100 
rownames(props_table_N4) <- NN_4_df$Var2
ggdf_N4 <- melt(data.frame(cluster = rownames(props_table_N4),props_table_N4, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N4_neighbor_props_allTLS_NR.pdf",height = 10, width = 10)
ggplot(ggdf_N4, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  # coord_polar(theta = "y") +
  # facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Non-Responder N4 Neighbors")
dev.off()


NN_5_df <- as.data.frame(table(aggr_allTLS_NR$cluster_annot, aggr_allTLS_NR$N5type))
NN_5_df <- as.data.frame(pivot_wider(NN_5_df, names_from = "Var1", values_from = "Freq"))
rownames(NN_5_df) <- NN_5_df$Var2
props_table_N5 <- (NN_5_df[,-1] / rowSums(NN_5_df[,-1])) * 100 
rownames(props_table_N5) <- NN_5_df$Var2
ggdf_N5 <- melt(data.frame(cluster = rownames(props_table_N5),props_table_N5, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N5_neighbor_props_allTLS_NR.pdf",height = 10, width = 10)
ggplot(ggdf_N5, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  # coord_polar(theta = "y") +
  # facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Non-Responder N5 Neighbors")
dev.off()


NN_6_df <- as.data.frame(table(aggr_allTLS_NR$cluster_annot, aggr_allTLS_NR$N6type))
NN_6_df <- as.data.frame(pivot_wider(NN_6_df, names_from = "Var1", values_from = "Freq"))
rownames(NN_6_df) <- NN_6_df$Var2
props_table_N6 <- (NN_6_df[,-1] / rowSums(NN_6_df[,-1])) * 100 
rownames(props_table_N6) <- NN_6_df$Var2
ggdf_N6 <- melt(data.frame(cluster = rownames(props_table_N6),props_table_N6, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N6_neighbor_props_allTLS_NR.pdf",height = 10, width = 10)
ggplot(ggdf_N6, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  # coord_polar(theta = "y") +
  # facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Non-Responder N6 Neighbors")
dev.off()


NN_7_df <- as.data.frame(table(aggr_allTLS_NR$cluster_annot, aggr_allTLS_NR$N7type))
NN_7_df <- as.data.frame(pivot_wider(NN_7_df, names_from = "Var1", values_from = "Freq"))
rownames(NN_7_df) <- NN_7_df$Var2
props_table_N7 <- (NN_7_df[,-1] / rowSums(NN_7_df[,-1])) * 100 
rownames(props_table_N7) <- NN_7_df$Var2
ggdf_N7 <- melt(data.frame(cluster = rownames(props_table_N7),props_table_N7, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N7_neighbor_props_allTLS_NR.pdf",height = 10, width = 10)
ggplot(ggdf_N7, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  #coord_polar(theta = "y") +
  #facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Non-Responder N7 Neighbors")
dev.off()


NN_8_df <- as.data.frame(table(aggr_allTLS_NR$cluster_annot, aggr_allTLS_NR$N8type))
NN_8_df <- as.data.frame(pivot_wider(NN_8_df, names_from = "Var1", values_from = "Freq"))
rownames(NN_8_df) <- NN_8_df$Var2
props_table_N8 <- (NN_8_df[,-1] / rowSums(NN_8_df[,-1])) * 100 
rownames(props_table_N8) <- NN_8_df$Var2
ggdf_N8 <- melt(data.frame(cluster = rownames(props_table_N8),props_table_N8, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N8_neighbor_props_allTLS_NR.pdf",height = 10, width = 10)
ggplot(ggdf_N8, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  #coord_polar(theta = "y") +
  #facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Non-Responder N8 Neighbors")
dev.off()


NN_9_df <- as.data.frame(table(aggr_allTLS_NR$cluster_annot, aggr_allTLS_NR$N9type))
NN_9_df <- as.data.frame(pivot_wider(NN_9_df, names_from = "Var1", values_from = "Freq"))
rownames(NN_9_df) <- NN_9_df$Var2
props_table_N9 <- (NN_9_df[,-1] / rowSums(NN_9_df[,-1])) * 100 
rownames(props_table_N9) <- NN_9_df$Var2
ggdf_N9 <- melt(data.frame(cluster = rownames(props_table_N9),props_table_N9, check.names = FALSE),
                id.vars = "cluster", value.name = "proportion", 
                variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N9_neighbor_props_allTLS_NR.pdf",height = 10, width = 10)
ggplot(ggdf_N9, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  #coord_polar(theta = "y") +
  #facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Non-Responder N9 Neighbors")
dev.off()


NN_10_df <- as.data.frame(table(aggr_allTLS_NR$cluster_annot, aggr_allTLS_NR$N10type))
NN_10_df <- as.data.frame(pivot_wider(NN_10_df, names_from = "Var1", values_from = "Freq"))
rownames(NN_10_df) <- NN_10_df$Var2
props_table_N10 <- (NN_10_df[,-1] / rowSums(NN_10_df[,-1])) * 100 
rownames(props_table_N10) <- NN_10_df$Var2
ggdf_N10 <- melt(data.frame(cluster = rownames(props_table_N10),props_table_N10, check.names = FALSE),
                 id.vars = "cluster", value.name = "proportion", 
                 variable.name = "neighbor")
pdf("../figPlots/NN/allTLS/N10_neighbor_props_allTLS_NR.pdf",height = 10, width = 10)
ggplot(ggdf_N10, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  #coord_polar(theta = "y") +
  #facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Non-Responder N10 Neighbors")
dev.off()

cols = colnames(NN_1_df[,-1])

NN_1_df <- add_missing_rows(df = NN_1_df, vec_rows = rows)
NN_2_df <- add_missing_rows(df = NN_2_df, vec_rows = rows)
NN_3_df <- add_missing_rows(df = NN_3_df, vec_rows = rows)
NN_4_df <- add_missing_rows(df = NN_4_df, vec_rows = rows)
NN_5_df <- add_missing_rows(df = NN_5_df, vec_rows = rows)
NN_6_df <- add_missing_rows(df = NN_6_df, vec_rows = rows)
NN_7_df <- add_missing_rows(df = NN_7_df, vec_rows = rows)
NN_8_df <- add_missing_rows(df = NN_8_df, vec_rows = rows)
NN_9_df <- add_missing_rows(df = NN_9_df, vec_rows = rows)
NN_10_df <- add_missing_rows(df = NN_10_df, vec_rows = rows)

all_NN_df <- NN_1_df

all_NN_df[, cols] = as.numeric(unlist(NN_1_df[, cols])) + as.numeric(unlist(NN_2_df[, cols])) + as.numeric(unlist(NN_3_df[, cols])) + as.numeric(unlist(NN_4_df[, cols])) + 
  as.numeric(unlist(NN_5_df[, cols])) + as.numeric(unlist(NN_6_df[, cols])) + as.numeric(unlist(NN_7_df[, cols])) + as.numeric(unlist(NN_8_df[, cols])) +
  as.numeric(unlist(NN_9_df[, cols])) + as.numeric(unlist(NN_10_df[, cols])) 

all_NN_df_small <- as.data.frame(all_NN_df[,-1])
rownames(all_NN_df_small) <- all_NN_df$Var2
props_table_100 <- (all_NN_df_small / rowSums(all_NN_df_small)) * 100 
props_table <- (all_NN_df_small / rowSums(all_NN_df_small)) 
(all_NN_df_small[1,1] / sum(all_NN_df_small[1,])) * 100 

#set up df for plotting

ggdf <- melt(data.frame(cluster = rownames(props_table_100),props_table_100, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "neighbor")

saveRDS(ggdf, file = "../figPlots/NN/allTLS/NN_propsTable_ggdf.RDS")

colnames(ggdf)

pdf("../figPlots/NN/allTLS/neighbor_props_allTLS_NR.pdf",height = 10,width = 10)
ggplot(ggdf, aes(x = "", y = proportion, fill = neighbor)) +
  geom_col() +
  coord_polar(theta = "y") +
  facet_wrap(~cluster) +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() #or theme_void
dev.off()

pdf("../figPlots/NN/allTLS/neighbor_props_barplot_allTLS_NR.pdf",height = 10,width = 10)
ggplot(ggdf, aes(x = cluster, y = proportion, fill = neighbor)) +
  geom_col() +
  scale_fill_manual(values = neighbor_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Non-Responder All 10 Neighbors Summed All TLS")
dev.off()

# Hole size
hsize <- 1.5

df <- ggdf %>% 
  mutate(x = hsize)

pdf("../figPlots/NN/allTLS/neighbor_props_donut_allTLS_NR.pdf",height = 10,width = 10)
ggplot(df, aes(x = hsize, y = proportion, fill = neighbor)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(0.2, hsize + 0.5)) +
  theme_minimal()+
  scale_fill_manual(values = neighbor_cols) +
  facet_wrap(~cluster)
dev.off()


