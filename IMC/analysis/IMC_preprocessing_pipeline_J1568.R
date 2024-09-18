# IMC preprocessing for J1568 IMC TLS panel

# .libPaths(paste0(work,"/Libraries"))
source("IMCpipeline_functions_packages.R")

####DATA LOADING####
# ## If loading from prior run
# #for from dropbox
# url <- "https://drive.google.com/uc?export=download&confirm=9iBg&id=1bplnOKUtSmrEtQCM-5srVCfTuC9WrTIQ"
# destfile <- paste0(work,"/backup_output.rds")
# download.file(url, destfile)
output<-readRDS('backup_output.rds')

## Read (skip if previously ran)
# output <- returnfcs(metaDataFile = paste0(work,"/Config/IMCCODA_metadata.xlsx"),
#                     panelDataFile = paste0(work,"/Config/IMCCODA_panel.xlsx"),
#                     dataDirectory = paste0(work,"/Data"))
allmarkers <- output$subtype_markers #subtype includes cluster and functional markers, basically all markers

## Returns:
## fcs (only marker expressions arcsin transformed), 
## fcs1 (arcsin transformed + spatial parameters), 
## fcs2 (scaled arcsin expr per flowframe)
## fcsraw (all raw data), and all marker/parameter lists

## Set up levels
samplevels= output$meta_data$sample_id
tissuelevels= output$meta_data$Location
responselevels= output$meta_data$Response

####DIAGNOSTICS####
## Spot check - number of cells per sample
cell_table <- table(output$sample_ids)
ggdf <- data.frame(sample_id = names(cell_table), 
                   cell_counts = as.numeric(cell_table))
ggdf$Structure <- output$meta_data$Location[match(ggdf$sample_id,output$meta_data$sample_id)]
#ggdf$timepoint <- factor(output$meta_data$Timepoint[match(ggdf$sample_id,output$meta_data$sample_id)], levels=timepointlevels)
#ggdf$timegroup <- factor(output$meta_data$TimeGroup[match(ggdf$sample_id,output$meta_data$sample_id)], levels=timegrouplevels)

ggp<-ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = sample_id)) + 
  geom_bar(stat = 'identity') + 
  #geom_text(aes(label = cell_counts), angle = 45, hjust = 0.5, vjust = -0.5, size = 2) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=5)) +  
  scale_x_discrete(drop = FALSE)
pdf('./figs/Diagnostics_cellcounts.pdf',width=6, height=4);ggp; dev.off()

## Multi-dimensional scaling plot to show similarities between samples
## Get the mean marker expression per sample
expr_mean_sample_tbl <- data.frame(sample_id = output$sample_ids, fsApply(output$fcs,exprs)) %>%
  group_by(sample_id) %>%  summarize_all(funs(mean))
expr_mean_sample <- t(expr_mean_sample_tbl[, -1])
colnames(expr_mean_sample) <- expr_mean_sample_tbl$sample_id
mds <- plotMDS(expr_mean_sample, plot = FALSE)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                   sample_id = colnames(expr_mean_sample))
ggdf$Location <- output$meta_data$Location[match(ggdf$sample_id,output$meta_data$sample_id)]
ggp<-ggplot(ggdf, aes(x = MDS1, y = MDS2, color = Location)) +
  #ggp<-ggplot(ggdf, aes(x = MDS1, y = MDS2, color = tumor, shape = timegroup)) +
  geom_point(size = 2.5) +
  #geom_text(aes(label = patient_id)) +
  theme_bw()+
  theme(plot.background = element_rect(fill="black"),
        panel.background = element_rect(fill="black"),
        panel.grid = element_blank(),
        axis.line = element_line(color="white"),
        axis.text = element_text(color="white"),
        axis.title = element_text(color="white"),
        legend.background = element_rect(fill="black"),
        legend.key = element_rect(fill="white"),
        legend.text = element_text(color="white"))+stat_ellipse()
pdf('./figs/Diagnostics_MDS_sample.pdf',width=6, height=6);ggp; dev.off()

## Colors for the heatmap
color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
pdf('./figs/Diagnostics_Heatmap.pdf',width=8, height=8)
pheatmap(expr_mean_sample_tbl[,allmarkers], color = color, display_numbers = TRUE,
         number_color = "black", fontsize_number = 3, 
         clustering_method = "average")
dev.off()
