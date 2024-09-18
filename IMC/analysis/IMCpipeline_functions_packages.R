##REQUIRED LIBRARIES##

library(premessa)
library(reshape2)
library(pals)
library(ggpattern)
library(viridis)
library(CoGAPS)
library(ggplot2)
library(Hmisc)
library(ComplexHeatmap)
library(ggiraphExtra)
library(diffcyt)
library(projectR)
library(corrplot)
library(psych)
library(ggvoronoi)
library(ggpubr)
library(multcomp)
library(sf)
library(clusterSim)
library(circlize)
library(RColorBrewer)
library(stringr)
library(igraph)
library(readxl)
library(dplyr)
library(packcircles)
library(gridExtra)
library(limma)
library(qgraph)
library(flowCore)
library(plot3D)
library(akima)
library(basetheme)
library(pheatmap)
library(premessa)
library(FlowSOM)
rm(list = ls())

# paneleditor_GUI(launch.browser = TRUE) #rename files to remove sampleID if required 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
work<-getwd()


####READ and CLUSTER FUNCTIONS####
returnfcs <- function(FDR_cutoff=.05,
                      metaDataFile=NULL,
                      panelDataFile=NULL,
                      dataDirectory=NULL,
                      shape_conditions=NULL,
                      color_conditions=NULL){
  ## This function generates an fcs file, subtype_markers, colors and shapes for clustering 
  require(scales);require(readxl);require(dplyr);require(flowCore)
  ## Directory and metadatafile checking
  if(!dir.exists(dataDirectory)) {stop('ERR: cannot find data directory')}
  if(!file.exists(metaDataFile)) {stop('ERR: cannot find metadata.xlsx or .csv file')}
  ## Read-in metadata and clean
  ifelse(grepl(metaDataFile,pattern='.xls'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
  md$file_name <- factor(md$file_name)
  md$File_order <- factor(md$File_order)
  md$Tumor <- factor(md$Tumor)
  
  
  rownames(md) = md$sample_id
  md$sample_id <- md$sample_id
  
  ## Make sure all files in metadata present in datadirectory
  if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])){
    print(paste('ERR: not all filenames in metadata present in data folder - missing',md$file_name[!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])],'Subsetting...'))
    md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])),]
  }
  
  ## Read fcs into fcs_raw
  fcs_raw <- read.flowSet(paste0(dataDirectory,"/",md$file_name), transformation = FALSE, truncate_max_range = FALSE)
  panel <- read_xlsx(panelDataFile)
  head(data.frame(panel))
  panel$Parameter <- gsub('-', '_', panel$Parameter)
  
  
  ## Export out the parameter/panel data from the flowFrame to edit
  ## use panel$Antigen to fix description in panel_fcs
  ## use metal+isotope as mapping between panel from xlsx and panel from the fcs files
  
  panel_fcs <- pData(parameters(fcs_raw[[1]]))
  panel_fcs$desc <- gsub('-', '_', panel_fcs$desc)
  panel_fcs$desc[is.na(panel_fcs$desc)] <- paste0('NA_',which(is.na(panel_fcs$desc)))  
  
  rownames(panel_fcs) = panel_fcs$name
  
  ## Replace desc with revised Name
  panel_fcs[panel$Parameter,]$desc<-panel$Name
  
  ## Replace parameter data in flowSet with edits
  pData(parameters(fcs_raw[[1]])) <- panel_fcs
  
  ## Assign objects to marker lists
  subtype_markers <- panel$Name[panel$Subtype == 1]
  functional_markers <- panel$Name[panel$Functional == 1]
  otherparameters <- panel$Name[panel$Other ==1]
  cluster_by <- panel$Name[panel$Cluster == 1]
  
  ## Check marker lists
  if(!all(subtype_markers %in% panel_fcs$desc)){stop('ERR: Not all subtype_markers in panel_fcs$desc')}
  if(!all(functional_markers %in% panel_fcs$desc)){stop('ERR: Not all functional_markers in panel_fcs$desc')}
  if(!all(otherparameters %in% panel_fcs$desc)){stop('ERR: Not all otherparameters in panel_fcs$desc')}
  if(!all(cluster_by %in% panel_fcs$desc)){stop('ERR: Not all cluster markers in panel_fcs$desc')}
  
  fcs <- fsApply(fcs_raw, function(x, cofactor = 0.8){
    colnames(x) <- panel_fcs$desc
    expr <- exprs(x)
    expr <- asinh(expr[, union(subtype_markers,functional_markers)] / cofactor)
    exprs(x) <- expr
    x
  })
  
  ## Save out the original rownames from the parameter list from the fcs flowFrame
  panellist <- rownames(pData(parameters(fcs[[1]])))
  
  ## Create dummy list to save all expression data
  exprTr_list<-c()
  
  ## Save arc sin transformed expression + spatial data for each flowFrame
  for(i in 1:length(md$file_name)){
    
    ## Expression data is combined with spatial parameter data
    
    exprRaw<-exprs(fcs_raw[[i]])
    
    colnames(exprRaw)<-panel_fcs$desc
    
    expr<-cbind(exprs(fcs[[i]])[, union(subtype_markers,functional_markers)],exprRaw[,otherparameters])
    
    ## Combine other (spatial) data with the protein data
    colnames(expr)<-c(colnames(exprs(fcs[[i]])),otherparameters)
    
    ## Filter out any event that is 95th percentile for BOTH CD29 and CD45 (antibody aggregates)
    ##expr<-expr[expr[,"CD29"] < quantile(expr[,"CD29"], probs=0.95) & expr[,"CD45"] < quantile(expr[,"CD45"], probs=0.95),]
    
    ## Add to list
    exprTr_list[[i]]<-expr
    
  }
  
  ## Create a flowSet based on the list of expression data
  fcs1<-flowSet(sapply(exprTr_list,flowFrame))
  
  ## Change parameter rownames
  panel_fcs1 <- pData(parameters(fcs1[[1]]))
  rownames(pData(parameters(fcs1[[1]]))) <- rownames(panel_fcs[panel_fcs$desc %in% pData(parameters(fcs1[[1]]))$name,])
  
  ###to scale every flowframe
  ## Save out the original rownames from the parameter list from the fcs flowFrame
  panellist <- rownames(pData(parameters(fcs[[1]])))
  
  ## Create dummy list to save all expression data
  exprTr_list<-c()
  
  ## Save arc sin transformed expression + spatial data for each flowFrame
  for(i in 1:length(md$file_name)){
    
    ## Expression data is combined with spatial parameter data
    
    expr<-exprs(fcs[[i]])
    
    expr<-t(scale(t(expr)))
    
    ## Add to list
    exprTr_list[[i]]<-expr
    
  }
  
  ## Create a flowSet based on the list of expression data
  fcs2<-flowSet(sapply(exprTr_list,flowFrame))
  
  
  ## Get sample ids
  sample_ids <- rep(md$sample_id, fsApply(fcs1, nrow))
  
  ## Return: 
  ## fcs (only marker expressions arcsin transformed), 
  ## fcs1 (arcsin transformed + spatial parameters), 
  ## fcs2 (scaled arcsin expr per flowframe)
  ## fcsraw (all raw data), and all marker/parameter lists
  return(list('fcs'=fcs,'fcs1'=fcs1,'fcs2'=fcs2,'fcsraw'=fcs_raw,'subtype_markers'=subtype_markers,'functional_markers'=functional_markers,'otherparameters'=otherparameters,'cluster_by'=cluster_by,'sample_ids'=sample_ids,'meta_data'=md))
}



clusterfcs <- function(fcs=output$fcs,
                       cluster_by = output$cluster_by,
                       seed=1234,plottitle='consensus_plots',
                       scaleoption,
                       numclusters=30){
  ## Cell population identification with FlowSOM and ConsensusClusterPlus
  require(dplyr);require(FlowSOM);require(ConsensusClusterPlus)
  set.seed(seed)
  som <- ReadInput(fcs, 
                   transform = FALSE, 
                   scale = scaleoption) %>% BuildSOM(colsToUse = cluster_by)
  
  ## Get the cell clustering into 100 SOM codes
  cell_clustering_som <- som$map$mapping[,1]
  
  ## Metaclustering into numclusters with ConsensusClusterPlus
  codes <- som$map$codes
  mc <- ConsensusClusterPlus(t(codes), maxK = numclusters, reps = 100,
                             pItem = 0.9, pFeature = 1, title = plottitle, 
                             plot = "png", clusterAlg = "hc", 
                             innerLinkage = "average", finalLinkage = "average",
                             distance = "euclidean", seed = 1234)
  
  ## Get cluster ids for each cell
  code_clustering <- mc[[numclusters]]$consensusClass#metaclusters consensus
  cell_clustering <- code_clustering[cell_clustering_som]#cell clustering from som
  return(list('code_clustering'=code_clustering,'cell_clustering'=cell_clustering,'metaclusters'=mc))
}


####DIAGNOSTIC HEATMAP FUNCTIONS ####
plot_clustering_heatmap_wrapper <- function(fcs, cell_clustering, nclusters,
                                            number_clusters, cluster_merging = NULL, 
                                            cluster_by=output$cluster_by,
                                            clusterMergeFile=NULL,
                                            fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales);require(pals);require(ComplexHeatmap)
  ## Will output the heatmap object and print it 
  color_clusters = kovesi.rainbow_bgyrm_35_85_c69(number_clusters)
  
  
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,cluster_by]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,cluster_by]
  
  
  ## Calculate the mean expression##################################################
  
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  ## Annotation for the merged clusters
  
  if(!is.null(clusterMergeFile)){
    ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Merged <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Merged <- color_clusters2
  }
  
  ## Colors for the heatmap
  
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  
  clusternames<-sort(unique(cluster_merging$new_cluster))
  
  colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
  
  names(colorassigned)<-clusternames
  
  color_list = list(clusters=colorassigned)
  
  color_list_byoriginal = colorassigned[match(cluster_merging$new_cluster,names(colorassigned))]
  
  cp<-rowAnnotation(clusters=cluster_merging$new_cluster,
                    col=color_list,
                    gp = gpar(col = "white", lwd = .5),
                    prop=anno_barplot(
                      clustering_prop, 
                      gp = gpar(fill=color_list_byoriginal, col=F),
                      border = F,
                      bar_width = 0.75, 
                      width = unit(2,"cm")))
  
  q <- Heatmap(expr_heat, name="scaled",
               col=rev(brewer.rdbu(100)),
               row_order = cluster_merging[order(cluster_merging$new_cluster),]$original_cluster,
               cluster_columns = T,
               cluster_rows = T,
               border = NA,
               rect_gp = gpar(col = "white", lwd = .5),
               right_annotation = cp,
               show_row_names = T,
               row_names_gp = gpar(fontsize=7),
               column_names_gp = gpar(fontsize=10),
               heatmap_legend_param = list(at=seq(from = 0, to = 1, by = 0.2)),
               width = unit(10, "cm"))
  
  print('Colors:')
  print(color_clusters)
  
  pdf(fileName, width=8, height=6) 
  
  return(q)
  
  dev.off() 
  
}


plot_clustering_heatmap_wrapper2 <- function(fcs, cell_clustering, nclusters=40,
                                             color_clusters=colorassigned,
                                             colorbar=rev(brewer.rdbu(100)),
                                             subtype_markers,
                                             fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  #labels_row <- expr01_mean$cell_clustering
  
  labels_row <- paste0(expr01_mean$cell_clustering, " ")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  p <- pheatmap(expr_heat, 
                color = colorbar, 
                cluster_cols = F,
                cluster_rows = F, 
                labels_row = labels_row,
                #scale="row",
                display_numbers = F, 
                number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8,
                border_color = "white",
                annotation_legend = F
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}


plot_clustering_heatmap_wrapper3 <- function(fcs, cell_clustering, nclusters,
                                             cluster_by=output$cluster_by,
                                             colorassigned,
                                             sampleno,
                                             clusterMergeFile=NULL,
                                             fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales);require(pals);require(ComplexHeatmap)
  ## Will output the heatmap object and print it 
  
  number_clusters <- length(unique(cell_clustering[which(output$sample_ids==output$meta_data$sample_id[sampleno])]))
  
  color_clusters = kovesi.rainbow_bgyrm_35_85_c69(number_clusters)
  
  #get expression
  expr <- exprs(fcs);expr <-expr[,cluster_by]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,cluster_by]
  
  ## Calculate the mean expression##################################################
  
  cell_clustering <- cell_clustering[which(output$sample_ids==output$meta_data$sample_id[sampleno])]
  
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  ## Annotation for the merged clusters
  
  if(!is.null(clusterMergeFile)){
    ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging1 <- read_excel(clusterMergeFile),cluster_merging1 <- read.csv(clusterMergeFile,header = TRUE))
    cluster_merging1 <- cluster_merging1[cluster_merging1$original_cluster %in% unique(cell_clustering),]
    cluster_merging1$new_cluster <- factor(cluster_merging1$new_cluster)
    annotation_row$Merged <- cluster_merging1$new_cluster
    color_clusters2 <- colorassigned[c(annotation_row$Merged)]
    names(color_clusters2) <- levels(cluster_merging1$new_cluster)
    annotation_colors$Merged <- color_clusters2
  }
  
  ## Colors for the heatmap
  
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  
  clusternames<-sort(unique(cluster_merging1$new_cluster))
  
  colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging1$new_cluster)))
  
  names(colorassigned)<-clusternames
  
  color_list = list(clusters=colorassigned)
  
  color_list_byoriginal = colorassigned[match(cluster_merging1$new_cluster,names(colorassigned))]
  
  cp<-rowAnnotation(clusters=cluster_merging1$new_cluster,
                    col=color_list,
                    gp = gpar(col = "white", lwd = .5),
                    prop=anno_barplot(
                      clustering_prop, 
                      gp = gpar(fill=color_list_byoriginal, col=F),
                      border = F,
                      bar_width = 0.75, 
                      width = unit(2,"cm")))
  
  q <- Heatmap(expr_heat,
               right_annotation = cp)
  
  r <- Heatmap(expr_heat, name="scaled",
               col=rev(brewer.rdbu(100)),
               #row_order = cluster_merging[order(cluster_merging$new_cluster),]$original_cluster,
               cluster_columns = T,
               cluster_rows = T,
               border = NA,
               rect_gp = gpar(col = "white", lwd = .5),
               right_annotation = cp,
               show_row_names = T,
               row_names_gp = gpar(fontsize=7),
               column_names_gp = gpar(fontsize=10),
               heatmap_legend_param = list(at=seq(from = 0, to = 1, by = 0.2)),
               width = unit(10, "cm"))
  
  print('Colors:')
  print(color_clusters)
  
  pdf(fileName, width=8, height=6) 
  
  return(q)
  
  dev.off() 
  
}



plot_clustering_heatmap_wrapper4 <- function(fcs, cell_clustering, nclusters=40,
                                             color_clusters=colorassigned,
                                             colorbar=rev(brewer.rdbu(100)),
                                             rowcut=rowcut,
                                             colcut=colcut,
                                             subtype_markers,
                                             fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.1, 0.9))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  #labels_row <- expr01_mean$cell_clustering
  
  labels_row <- paste0(expr01_mean$cell_clustering, " ")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  p <- pheatmap(expr_heat, 
                color = colorbar, 
                cluster_cols = F,
                cluster_rows = F, 
                labels_row = labels_row,
                #scale="row",
                display_numbers = F, 
                number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 10,
                cellheight = 10,
                border_color = NA,
                annotation_legend = F,
                gaps_col = colcut,
                gaps_row = rowcut
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}



####UMAP####
#separate UMAP also created in Giotto
do_umap <- function(fcs,subtype_markers,sample_ids,cell_clustering,metadata,
                    clusterMergeFile='~/Desktop/ViralHCC/ViralHCC_merging.xlsx',
                    seed = 1234, ncells=2000,sample_subset=NULL){
  require(umap);require(flowCore);require(readxl)
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Create vector to later find and skip duplicates
  dups <- duplicated(expr[, subtype_markers])
  dups <- which(!(dups))## Find and skip duplicates
  ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
  ## New clustering1m
  mm <- match(cell_clustering, cluster_merging$original_cluster)
  cell_clustering1m <- cluster_merging$new_cluster[mm]
  ## Create a data frame of sample_ids and cell_clustering1m
  dtf<-data.frame(ids=sample_ids,type=cell_clustering1m)
  #dtf$B<- dtf$type!="B"#add a column that indicates whether the cell is a B cell or not; TRUE is non-B
  ##Why exclude B CELLS?
  ## WE HAVE NO B CELLS bc dtf$B depends on dtf$type depends on cellclustering1m which is just numbers 1:40 so..?
  #should it be the named parts cluster in merge file corresponding to it like if 30 -> grepl(cluster_merging[30,3],pattern='^B')??
  #Im blocking this out till we know why we have to do this
  ## Create subset columns without B cells (ONLY to generate the correct # of columns in inds2 object)
  #sample_ids2 <- dtf$ids[dtf$type!="B"] #sampleids without B cells
  #cell_clustering1m2 <- dtf$type[dtf$type!="B"] #clusters without B cells
  ## Data subsampling: create indices by sample
  inds <- split(1:length(sample_ids), sample_ids) #to get original indexes belonging to each cluster
  #inds2 <- split(1:length(sample_ids2), sample_ids2) #to fill in the original indexes that do not have B cells
  samplenames <- names(inds) #create a name vector of the files
  #FOR THIS TO WORK MUST BE IN FORMAT PRE/POST Tx
  # for (i in 1:(length(samplenames)/2)){#1:15 was here because ids in dtf was 30 factors and could be divided into 0 and 6wk for each so changed it
  #   templength <- length(inds2[[i]])
  #   inds2[[i]] <- inds[[i]][dtf$B[dtf$ids==samplenames[i]]] #subsets the "B cell or not a B cell" column for each sample by the name
  #   inds2[[i]] <- inds2[[i]][1:templength]
  # }
  
  custom.settings = umap.defaults
  custom.settings$seed = seed
  
  #custom.settings$n.neighbors = neighbors
  ####umapindex generation####
  #umap ncells = table of sample ids with how many to downsample to by default col = id, row = ncells
  #sample ids = chr [1:2851129] "4927_0wk" "4927_0wk" "4927_0wk" "4927_0wk" ...
  #^ from ## Generate sample IDs corresponding to each cell in the 'expr' matrix sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
  #can subset sample_ids and rerun umap 
  #can do timepoint or patient number or pre/post if you know corresp sample ids
  #sample_subset = '02_0wk' or c('02_0wk','2973_0wk') for example picking all 0 wk ones or use regex sample_ids[(grepl(sample_ids,pattern = '0wk'))]
  ifelse(is.null(sample_subset),
         umap_ncells <- pmin(table(sample_ids), ncells),
         umap_ncells <- pmin(table(sample_ids), ncells)[sample_subset]
  )
  if(!is.null(sample_subset)){inds <- inds[sample_subset]}
  umap_inds <- lapply(names(inds), function(i){
    s <- sample(inds[[i]], umap_ncells[i], replace = FALSE)
    intersect(s, dups)
  })
  set.seed(seed)
  umap_inds <- unlist(umap_inds)
  umap_out <- umap(expr[umap_inds, subtype_markers], config = custom.settings, method = 'naive')
  umapRes2D = data.frame(umap1 = umap_out$layout[, 1], umap2 = umap_out$layout[, 2], 
                         expr[umap_inds, subtype_markers],
                         sample_id = sample_ids[umap_inds], cell_clustering = factor(cell_clustering1m[umap_inds]), check.names = FALSE)
  
  #exclude any unassigned cluster post umap if needed--this has to be done by looking at the two columns 
  #metadata$sampleid is just a number in this metadatafile so to make unique ones combine samp_id col with timepoint (0wk)
  #to get in format of umapres2d$sample_id which looks like "02_0wk" do:
  return(umapRes2D)
}


plotUmap <- function(umapRes,seed=1234,neighbors=10,midpoint,color_clusters=colorassigned,code_clustering,subtype_markers=NULL)
{require(umap);require(ggplot2);require(viridis);require(ggrepel)
  custom.settings = umap.defaults
  custom.settings$seed = seed
  custom.settings$n.neighbors = neighbors
  ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = cell_clustering)) +
    #geom_point(size = 1) +
    geom_text(aes(label=cell_clustering), size=0.5)+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    
    scale_color_manual(values = color_clusters, name="CLUSTERS") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp)
  
  # print(ggp + facet_wrap(~ case, ncol=4)+ggtitle('CASE'))
  print(ggp + facet_wrap(~ location, ncol=4)+ggtitle('Location'))
  # print(ggp + facet_wrap(~ timepoint, ncol=3)+ggtitle('TIMEPOINT'))
  # print(ggp + facet_wrap(~ timegroup, ncol=2)+ggtitle('TIMEGROUP'))
  # 
  
  ggp2 <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = sample_id)) +
    geom_point(size = 1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank()
    ) +
    
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp2)
  #can specify which markers to display
  if(!is.null(subtype_markers)){
    for(i in subtype_markers)
    {
      ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = umapRes[,i])) +
        geom_point(size = 1) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        scale_color_gradient2(i, low="dark blue",mid="white",high="dark red", midpoint = mean(unlist(umapRes[,i])))
      print(ggp)
    }
  }
}


plot_cluster_composition <- function(ggdf, group, grouplabel){
  gp <- ggplot(ggdf, aes(x = group, y = proportion, fill=cluster, order=cluster)) +
    geom_bar(stat = "identity", position="fill", width=0.8) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
          axis.text.y = element_text(color="black"),
          axis.ticks = element_line(size=0.25),
          strip.text=element_text(size=8),
          strip.background = element_rect(fill=NA, color=NA),
          legend.key.size = unit(.75,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white"),
    ) +
    ylab("% of Total Cells")+
    xlab(grouplabel)+
    scale_fill_manual(values = colorassigned,
                      breaks = clusterlevels,
                      labels = clusterlevels)+
    scale_y_continuous(expand = c(0,0))+
    guides(fill=guide_legend(ncol=2))
  return(gp)
}

plotCells <- function(object,marker,title,save){ # plots UMAPs or saves UMAP file
  if(save=="yes"){
    tiff(paste0(title," ",marker,".tiff"))
    print(ggplot(data.frame(reducedDims(object)$UMAP),  aes(x = X1, y = X2, color = t(exprs(object))[,marker] ))+
            geom_point(size = 1) +
            theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())+
            scale_color_viridis(discrete = FALSE, option = "D")+
            scale_fill_viridis(discrete = FALSE)+ggtitle(marker)+ylab("UMAP2") + xlab("UMAP1")+labs(color = "Expression"))
    dev.off()
  } else {
    ggplot(data.frame(reducedDims(object)$UMAP),  aes(x = X1, y = X2, color = t(exprs(object))[,marker] ))+
      geom_point(size = 1) +
      theme_minimal()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      scale_color_viridis(discrete = FALSE, option = "D")+
      scale_fill_viridis(discrete = FALSE)+ggtitle(title)+ylab("UMAP2") + xlab("UMAP1")+labs(color = "Expression")
  }
}

trajInfer <- function(object,title){ #learns pseudotime, user picks start node, saves figure
  object <-cluster_cells(object,resolution=1e-5) 
  object<-learn_graph(object)
  object<-order_cells(object)
  pData(object)$pseudotime <- pseudotime(object)
  tiff(paste0(title,".tiff"))
  print(plot_cells(object, color_cells_by = "pseudotime",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE, graph_label_size=1.5))
  dev.off()
  pData(object)$subset <- title
  return(object)
}

densityAUCratio <- function(data,sample){ #returns the ratios of AUCs of each pole 
  r <- data.frame(pseudotime=density(data[data$sample_id == sample,]$pseudotime)$x,
                  density=density(data[data$sample_id == sample,]$pseudotime)$y)
  timebins <- sort(unique(cut(r$pseudotime,breaks=2))) #two bins for two poles 
  timebins <- gsub("]","",timebins)
  timebins<- gsub("[()]", "", timebins)
  timebins <- str_split(timebins, ",", simplify = TRUE)
  return(auc(r$pseudotime, r$density,from = as.numeric(timebins[1,2]), to = max(r$pseudotime))/auc(r$pseudotime, r$density,from = min(r$pseudotime), to = as.numeric(timebins[1,2]))) #ratio indicating cell density bias 
}

plot_pseudotime_umap <- function(object,pseudotime,title){
  ggplot(data.frame(reducedDims(object)$UMAP),  aes(x = X1, y = X2, color = as.numeric(pseudotime)))+
    geom_point(size = 1) +
    theme_minimal()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_color_viridis(discrete = FALSE, option = "D")+
    scale_fill_viridis(discrete = FALSE)+ggtitle(title)+ylab("UMAP2") + xlab("UMAP1")+labs(color = "Pseudotime")
}

save_pseudotime_umap <- function(object,pseudotime,title){
  tiff(paste0(title," ",".tiff"))
  print(ggplot(data.frame(reducedDims(object)$UMAP),  aes(x = X1, y = X2, color = as.numeric(pseudotime)))+
          geom_point(size = 1) +
          theme_minimal()+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())+
          scale_color_viridis(discrete = FALSE, option = "D")+
          scale_fill_viridis(discrete = FALSE)+ggtitle(title)+ylab("UMAP2") + xlab("UMAP1")+labs(color = "Pseudotime"))
  dev.off()
}

DE_model <- function(object,formula,expression_fam){
  gene_fits <- fit_models(object, model_formula_str = formula,expression_family = expression_fam) #anything you want
  fit_coefs <- coefficient_table(gene_fits)
  terms <- fit_coefs %>% filter(!term == "(Intercept)") ## remove intercept rows
  terms <- terms %>% mutate(q_value = p.adjust(p_value))
  terms <- data.frame(terms)
  #terms <- terms %>% filter (q_value < 0.05) 
  remove<-which(colnames(terms) == "model_summary" | colnames(terms) == "model")
  terms<-terms[,-c(remove)] 
  return(terms)
}

rsq <- function (x, y) cor(x, y) ^ 2

plot_cluster_composition_pie <- function(ggdf, group, grouplabel){
  gp <- ggplot(ggdf, aes(x = group, fill=cluster)) +
    theme_bw() +
    geom_bar(width = 1, position = "fill") +coord_polar(theta = "y")
  return(gp)
}

plotNMF <- function(object,pattern,title,save){ #  UMAP with NMF pattern weights overlayed 
  if(save=="yes"){
    tiff(paste0(title,".tiff"))
    print(ggplot(data.frame(reducedDims(object)$UMAP),  aes(x = X1, y = X2, color = pattern))+
            geom_point(size = 1) +
            theme_classic()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())+
            scale_color_viridis(discrete = FALSE, option = "D")+
            scale_fill_viridis(discrete = FALSE)+ggtitle(title)+ylab("UMAP2") + xlab("UMAP1")+labs(color = "Pattern Weight"))
    dev.off()
  } else {
    ggplot(data.frame(reducedDims(object)$UMAP),  aes(x = X1, y = X2, color = pattern))+
      geom_point(size = 1) +
      theme_classic()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      scale_color_viridis(discrete = FALSE, option = "D")+
      scale_fill_viridis(discrete = FALSE)+ggtitle(title)+ylab("UMAP2") + xlab("UMAP1")+labs(color = "Pattern Weight")
  }
}


arrangeSinglePatientROIs <- function(aggr_flowframe, patient){
  
  # ARRANGE THE ROIS FROM A SINGLE SECTION TO PLOT
  # ORIGINAL CODE BEFORE LOOP:
  # ## These are all the ROIs from A11 
  # sdfA11.1 <- aggr_flowframe[aggr_flowframe$ROI_ID == "1.1-1",]
  # sdfA11.2 <- aggr_flowframe[aggr_flowframe$ROI_ID == "1.2-1",]
  # sdfA11.2$X_position <- sdfA11.2$X_position + max(sdfA11.1$X_position) +50
  # sdfA11.3 <- aggr_flowframe[aggr_flowframe$ROI_ID == "1.2-2",]
  # sdfA11.3$X_position <- sdfA11.3$X_position  + max(sdfA11.2$X_position) +50
  # sdfA11 <- rbind(sdfA11.1,sdfA11.2,sdfA11.3)
  
  #NEW:
  
  ptff <- aggr_flowframe[aggr_flowframe$PatientID == patient,]
  numrois <- length(unique(ptff$ROI_ID))
  rois <- unique(ptff$ROI_ID)
  roivec <- c(1:numrois)[-1]# need this so you don't change the first ROI coords
  if (numrois <= 3){
    for (r in roivec){
      ptff[ptff$ROI_ID == rois[r],]$X_position <- ptff[ptff$ROI_ID == rois[r],]$X_position + max(ptff[ptff$ROI_ID == rois[r-1],]$X_position) + 50
    }
  }
  
  if (numrois  == 5){
    roivec <- c(2,3)# need this so you don't change the first ROI coords
    for (r in roivec){
      ptff[ptff$ROI_ID == rois[r],]$X_position <- ptff[ptff$ROI_ID == rois[r],]$X_position + max(ptff[ptff$ROI_ID == rois[r-1],]$X_position) + 50
    }
    ptff[ptff$ROI_ID == rois[5],]$X_position <- ptff[ptff$ROI_ID == rois[5],]$X_position + max(ptff[ptff$ROI_ID == rois[4],]$X_position) + 50
    ptff[ptff$ROI_ID == rois[1]| ptff$ROI_ID == rois[2]|ptff$ROI_ID ==  rois[3],]$Y_position <- ptff[ptff$ROI_ID == rois[1]|ptff$ROI_ID ==  rois[2]|ptff$ROI_ID ==  rois[3],]$Y_position + max(ptff[ptff$ROI_ID == rois[4],]$Y_position) + 50
    
  }
  
  # if (numrois  = 5){
  # 
  #   
  # }
  return(ptff)
}


arrangeMultiplePatientROIs <- function(aggr_flowframe, rois){
  
  # ARRANGE THE ROIS FROM A SINGLE SECTION TO PLOT
  # ORIGINAL CODE BEFORE LOOP:
  # ## These are all the ROIs from A11 
  # sdfA11.1 <- aggr_flowframe[aggr_flowframe$ROI_ID == "1.1-1",]
  # sdfA11.2 <- aggr_flowframe[aggr_flowframe$ROI_ID == "1.2-1",]
  # sdfA11.2$X_position <- sdfA11.2$X_position + max(sdfA11.1$X_position) +50
  # sdfA11.3 <- aggr_flowframe[aggr_flowframe$ROI_ID == "1.2-2",]
  # sdfA11.3$X_position <- sdfA11.3$X_position  + max(sdfA11.2$X_position) +50
  # sdfA11 <- rbind(sdfA11.1,sdfA11.2,sdfA11.3)
  
  #NEW:
  
  ptff <- aggr_flowframe[aggr_flowframe$ROI_ID %in% rois,]
  numrois <- length(unique(ptff$ROI_ID))
  rois <- unique(ptff$ROI_ID)
  roivec <- c(1:numrois)[-1]# need this so you don't change the first ROI coords
  if (numrois <= 3){
    for (r in roivec){
      ptff[ptff$ROI_ID == rois[r],]$X_position <- ptff[ptff$ROI_ID == rois[r],]$X_position + max(ptff[ptff$ROI_ID == rois[r-1],]$X_position) + 50
    }
  }
  
  if (numrois  == 4){
    ptff[ptff$ROI_ID == rois[2],]$X_position <- ptff[ptff$ROI_ID == rois[2],]$X_position + max(ptff[ptff$ROI_ID == rois[1],]$X_position) + 50
    ptff[ptff$ROI_ID == rois[4],]$X_position <- ptff[ptff$ROI_ID == rois[4],]$X_position + max(ptff[ptff$ROI_ID == rois[3],]$X_position) + 50
    ptff[ptff$ROI_ID == rois[1]| ptff$ROI_ID == rois[2],]$Y_position <- ptff[ptff$ROI_ID == rois[1]|ptff$ROI_ID ==  rois[2],]$Y_position + max(ptff[ptff$ROI_ID == rois[4],]$Y_position) + 50
    
  }
  
  if (numrois  == 5){
    roivec <- c(2,3)# need this so you don't change the first ROI coords
    for (r in roivec){
      ptff[ptff$ROI_ID == rois[r],]$X_position <- ptff[ptff$ROI_ID == rois[r],]$X_position + max(ptff[ptff$ROI_ID == rois[r-1],]$X_position) + 50
    }
    ptff[ptff$ROI_ID == rois[5],]$X_position <- ptff[ptff$ROI_ID == rois[5],]$X_position + max(ptff[ptff$ROI_ID == rois[4],]$X_position) + 50
    ptff[ptff$ROI_ID == rois[1]| ptff$ROI_ID == rois[2]|ptff$ROI_ID ==  rois[3],]$Y_position <- ptff[ptff$ROI_ID == rois[1]|ptff$ROI_ID ==  rois[2]|ptff$ROI_ID ==  rois[3],]$Y_position + max(ptff[ptff$ROI_ID == rois[4],]$Y_position) + 50
    
  }
  
  if (numrois  == 6){
    roivec <- c(2,3)# need this so you don't change the first ROI coords
    for (r in roivec){
      ptff[ptff$ROI_ID == rois[r],]$X_position <- ptff[ptff$ROI_ID == rois[r],]$X_position + max(ptff[ptff$ROI_ID == rois[r-1],]$X_position) + 50
    }
    ptff[ptff$ROI_ID == rois[5],]$X_position <- ptff[ptff$ROI_ID == rois[5],]$X_position + max(ptff[ptff$ROI_ID == rois[4],]$X_position) + 50
    ptff[ptff$ROI_ID == rois[6],]$X_position <- ptff[ptff$ROI_ID == rois[6],]$X_position + max(ptff[ptff$ROI_ID == rois[5],]$X_position) + 50
    ptff[ptff$ROI_ID == rois[1]| ptff$ROI_ID == rois[2]|ptff$ROI_ID ==  rois[3],]$Y_position <- ptff[ptff$ROI_ID == rois[1]|ptff$ROI_ID ==  rois[2]|ptff$ROI_ID ==  rois[3],]$Y_position + max(ptff[ptff$ROI_ID == rois[4],]$Y_position) + 50
    
  }
  
  # if (numrois  = 5){
  # 
  #   
  # }
  return(ptff)
}

clusterToPattern <- function(
  clusters, # a vector of character cluster assignments
  data # data used to make clusters object
){
  nD<-length(unique(clusters))
  nG<-dim(data)[1]
  tempP<-matrix(data=rep(0,nD*nG),nrow = nG,ncol =nD)
  rownames(tempP)<-rownames(data)
  colnames(tempP)<-unique(clusters)
  #for(x in 1:nD) {tempP[Patterns$cluster==x,x]<-rowMeans(data[Patterns$cluster==x,])}
  for(x in unique(clusters)) 
  {
    mean_cluster <- colMeans(data[,clusters==x])
    tempP[, x]<-apply(data[,clusters==x], 1, cor ,mean_cluster)
  }
  Patterns<-tempP
  return(Patterns)
}

