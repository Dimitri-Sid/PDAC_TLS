# - - - - - - - - - - - - 
# functions and packages
# Author: Dimitri Sidiropoulos
# - - - - - - - - - - - - 

Sys.unsetenv("GITHUB_PAT")

##REQUIRED LIBRARIES##
library("rtracklayer")
library("pryr")
library(readxl)
library(stringr)
library(msigdbr)
library(raster)
library(plotrix)
library(dendextend)
library(circlize)
library(stringdist)
library(dplyr)
library(ggplot2)
library(viridis)
library(Rfast)
library(DESeq2)
library(msigdbr)
library(ggpubr)
library(ggrepel)
library(SpaceMarkers)
library(fgsea)
library(corrplot)
#library(Matrix.utils)
library(dplyr)
library("rsconnect")
library("GEOquery")
library("Matrix")
library("scater")
library("plyr")
library("ggplot2")
library(escape)
library(dittoSeq)
library("Seurat")
library(harmony)
library(boot)
library(tibble)
library(dplyr)
library(reshape2)
library(ggplot2)
library(argparse)
library(patchwork)
library(cowplot)
library(rhdf5)
library(sctransform)
library(tidyverse)
library(ggrepel)
library(projectR)
library(Rfast)
library("STdeconvolve")
library("splitstackshape")
library("limma")
library("ComplexHeatmap")
#library(spacexr)
library("ggpubr")
library("pheatmap")
library("grid")
library("gridExtra")
library("UCell")
library("splitstackshape")
library("corrplot")
library("BiocFileCache")
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)

enrichGOintersectingGenes <- function(smarkers){
  smgenes <- unique(smarkers[smarkers[,2]=="vsBoth",]$Gene)
  gene.df <- bitr(smgenes, fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  ego <- enrichGO(gene  = gene.df$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENTREZID',
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  maxGSSize = 300)
  dotplot(ego, showCategory=10,  font.size = 10)
}

enrichKEGGintersectingGenes <- function(smarkers){
  smgenes <- unique(smarkers[smarkers[,2]=="vsBoth",]$Gene)
  gene.df <- bitr(smgenes, fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  kk <- enrichKEGG(gene = gene.df$ENTREZID)
  dotplot(kk, showCategory=10)
}


gseGOintersectingGenes <- function(smarkers){
  smarkers <- smarkers[smarkers[,2]=="vsBoth",]
  smarkers_subset <- smarkers[, c('Gene', 'SpaceMarkersMetric' )]
  smarkers_collapsed <- aggregate(.~Gene,data=smarkers_subset, FUN=mean )

  geneList <- smarkers_collapsed$SpaceMarkersMetric
  names(geneList) <- smarkers_collapsed$Gene

  geneList <- sort(geneList, decreasing = TRUE)
  
  gene.df <- bitr(names(geneList), fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  
  rownames(gene.df) <- gene.df$SYMBOL
  
  names(geneList) <- gene.df[names(geneList),2]

  gseaResult <- gseGO(geneList     = geneList,
                      ont          = "BP", # Biological Process. You can also use "CC" for Cellular Component, or "MF" for Molecular Function
                      OrgDb        = org.Hs.eg.db,
                      nPerm        = 1000, # Number of permutations. Increase for more accurate results, but it will take longer.
                      minGSSize    = 10, # Minimum size of gene sets considered
                      maxGSSize    = 500, # Maximum size of gene sets considered
                      pAdjustMethod = "BH", # Method for adjusting p values
                      verbose      = FALSE) # Set to TRUE to see progress
  dotplot(gseaResult) + ggplot2::theme_minimal()
}



unitVector <- function(n, length)
{
  vec <- rep(0, length)
  vec[n] <- 1
  return(vec)
}

patternMarkers_corrected <- function(object, threshold, lp, axis)
{
  ## check inputs to the function
  if (!(threshold %in% c("cut", "all")))
    stop("threshold must be either 'cut' or 'all'")
  if (!is.na(lp) & length(lp) != ncol(object@featureLoadings))
    stop("lp length must equal the number of patterns")
  if (!(axis %in% 1:2))
    stop("axis must be either 1 or 2")
  ## need to scale each row of the matrix of interest so that the maximum is 1
  resultMatrix <- if (axis == 1) object@featureLoadings else object@sampleFactors
  normedMatrix <- t(apply(resultMatrix, 1, function(row) row / max(row)))
  ## handle the case where a custom linear combination of patterns was passed in
  if (!is.na(lp))
  {
    markerScores <- apply(normedMatrix, 1, function(row) sqrt(sum((row-lp)^2)))
    markersByPattern <- names(sort(markerScores, decreasing=FALSE, na.last=TRUE))
    return(list(
      "PatternMarkers"=markersByPattern,
      "PatternMarkerRanks"=rank(markerScores),
      "PatternMarkerScores"=markerScores
    ))
  }
  ## default pattern marker calculation, each pattern has unit weight
  markerScores <- sapply(1:ncol(normedMatrix), function(patternIndex)
    apply(normedMatrix, 1, function(row)
    {
      lp <- unitVector(patternIndex, ncol(normedMatrix))
      return(sqrt(sum((row-lp)^2)))
    })
  )
  markerRanks <- apply(markerScores, 2, rank)
  colnames(markerScores) <- colnames(markerRanks) <- colnames(normedMatrix)
  ## keep only a subset of markers for each pattern depending on the type of threshold
  if (threshold == "cut") # all markers which achieve minimum rank
  {
    simplicityGENES <- function(As, Ps) {
      # rescale p's to have max 1
      pscale <- apply(Ps,1,max)
      
      # rescale A in accordance with p's having max 1
      As <- sweep(As, 2, pscale, FUN="*")
      
      # find the A with the highest magnitude
      Arowmax <- t(apply(As, 1, function(x) x/max(x)))
      pmax <- apply(As, 1, max)
      
      # determine which genes are most associated with each pattern
      ssl <- matrix(NA, nrow=nrow(As), ncol=ncol(As),
                    dimnames=dimnames(As))
      for (i in 1:ncol(As)) {
        lp <- rep(0, ncol(As))
        lp[i] <- 1
        ssl.stat <- apply(Arowmax, 1, function(x) sqrt(t(x-lp)%*%(x-lp)))
        ssl[order(ssl.stat),i] <- 1:length(ssl.stat)
      }
      
      return(ssl)
      
    }
    simGenes <- simplicityGENES(As=object@featureLoadings,
                                Ps=object@sampleFactors)
    
    patternMarkers <- list()
    
    nP <- ncol(simGenes)
    
    for (i in 1:nP) {
      
      sortSim <- names(sort(simGenes[,i],decreasing=F))
      
      geneThresh <- min(which(simGenes[sortSim,i] > 
                                apply(simGenes[sortSim,],1,min)))
      
      markerGenes <- sortSim[1:geneThresh]
      
      markerGenes <- unique(markerGenes)
      
      patternMarkers[[i]] <- markerGenes
      
      markersByPattern <- patternMarkers
      
    }
  }
  else if (threshold == "all") # only the markers with the lowest scores
  {
    patternsByMarker <- colnames(markerScores)[apply(markerScores, 1, which.min)]
    markersByPattern <- sapply(colnames(markerScores), USE.NAMES=TRUE, simplify=FALSE,
                               function(pattern) rownames(markerScores)[which(patternsByMarker==pattern)])
    
  }
  return(list(
    "PatternMarkers"=markersByPattern,
    "PatternMarkerRanks"=markerRanks,
    "PatternMarkerScores"=markerScores
  ))
}


plot_dendrogram <- function(file_path, clone_id, num_colors) {
  #a <- load('C:/Users/sunny/Downloads/22_06_01_pdac-visium-v1.rda')
  #a <- load('C:/Users/sunny/Downloads/changeo_annotated.rda')
  a <- load(file_path)
  ##igh <- res[res$locus == 'IGH',]
  #igk_df <- changeo[changeo$locus == "IGK",]
  #igk_6013_df <- igk_df[igk_df$clone_id == "IGK_7199",]
  cid_df <- changeo[changeo$clone_id == clone_id,]
  # Take counts and cdr3 to collapse rows
  #igk_6013_subset <- igk_6013_df[, c('cdr3', 'count')]
  cid_subset <- cid_df[, c('cdr3', 'count')]
  #igk_6013_collapsed <- aggregate(.~cdr3,data=igk_6013_subset, FUN=sum )
  cid_collapsed <- aggregate(.~cdr3,data=cid_subset, FUN=sum )
  
  #d <- data.frame(matrix(ncol=(nrow(igk_6013_collapsed) + 1), nrow=nrow(igk_6013_collapsed)))
  d <- data.frame(matrix(ncol=(nrow(cid_collapsed) + 1), nrow=nrow(cid_collapsed)))
  ##cols <- igk_6013_df$cdr3_aa
  #cols <- igk_6013_collapsed$cdr3
  cols <- cid_collapsed$cdr3
  ##colnames(d) <- c('cdr3_aa', cols)
  ##d$cdr3_aa <- igk_6013_df$cdr3_aa
  
  colnames(d) <- c('cdr3', cols)
  #d$cdr3 <- igk_6013_collapsed$cdr3
  d$cdr3 <- cid_collapsed$cdr3
  for (i in 2:ncol(d)) {
    d[, i] <- stringdist(cols[i - 1], d[, 1], method='lv')
  }
  # Load distance matrix
  #d <- dist('dist matrix')
  
  dist <- as.dist(d[, 2:ncol(d)])
  
  #igk$cdr3[igk$barcode == 'AATTGCAGCAATCGAC']
  # Create distance matrix using hierarchical clustering
  hc <- as.dendrogram(hclust(dist))
  
  # Perform labeling of counts, branches, etc.
  hc <- color_branches(hc)
  hc <- hc %>% set('branches_lwd', 3)
  
  #counts = igk_6013_collapsed$count
  counts = cid_collapsed$count
  #  counts = c(1, 30, 50, 75, 100, 120, 150, 170, 200, 210)
  
  #breaks = seq(min(counts), max(counts), length.out=11)
  breaks = seq(min(counts), max(counts), length.out=num_colors+1)
  #colors <- viridis(10)
  colors <- viridis(num_colors)
  labels <- hc %>% labels
  label_color <- c()
  label_ct <- 1
  for (label in labels) {
    color <- NULL
    break_ct <- 2
    while (is.null(color)) {
      if (counts[label_ct] <= breaks[break_ct]) {
        color <- colors[break_ct - 1]
      }
      break_ct <- break_ct + 1
    }
    label_color <- append(label_color, color)
    label_ct <- label_ct + 1
  }
  cid_collapsed$color <- label_color
  
  hc <- hc %>% set('leaves_pch', 19) %>% set('leaves_cex', 2) %>% set('leaves_col', label_color)
  # Circular dendrogram
  #hcircle <- circlize_dendrogram(hc)
  
  #pdf('changeo_dendrogram.pdf')
  circlize_dendrogram(hc)
  return(list('dend' = hc, 'cid_collapsed' = cid_collapsed))
}

# Function to find top cell types
find_top_cell_types <- function(row) {
  sorted_indices <- order(row, decreasing = TRUE)
  top_value <- row[sorted_indices[1]]
  
  if (top_value > 90) {
    return(names(row)[sorted_indices[1]])
  } else {
    return(paste(names(row)[sorted_indices[1:2]], collapse = "-"))
  }
}

distanceDE <- function(tls_subset, targetD) {
  
  if (targetD == "distance2PDAC"){
    ddsSE <- DESeqDataSetFromMatrix(as.matrix(as.matrix(tls_subset@assays$Spatial@counts)+1), DataFrame(tls_subset@meta.data), design = ~distance2PDAC + slide)
  }
  
  if (targetD == "distance2Islets"){
    ddsSE <- DESeqDataSetFromMatrix(as.matrix(as.matrix(tls_subset@assays$Spatial@counts)+1), DataFrame(tls_subset@meta.data), design = ~distance2Islets + slide)
  }
  
  ddsSE <- DESeq(ddsSE)
  res2 <- results(ddsSE, tidy = TRUE,name= targetD )
  coefs <- coef(ddsSE)
  res2$gene_short_name <- res2$row
  res2 <- res2[order(res2$stat, decreasing = TRUE),]
  res2$index <- 1:length(res2$gene_short_name)
  print("Top Genes:")
  print(topGenes <- res2$gene_short_name[c(1:10,(length(res2$gene_short_name)-10):length(res2$gene_short_name))]) 
  return(res2)
}

distanceGSEA <- function(tls_subset, res, human.genes, genesets.interest ) {
  pathways.interest <- genesets.interest %>% split(x = .$gene_symbol, f = .$gs_name)
  stats <- res$stat 
  names(stats) <- res$gene_short_name 
  gsres <- fgsea(pathways = pathways.interest,
                 stats = stats) %>%
    arrange(desc(NES)) %>% 
    mutate(log10pval = -log10(pval)) %>%
    dplyr::select(pathway, pval, log10pval, everything())
  gsres$fdr <- p.adjust(gsres$pval)
  return(gsres)
}


# Function to find the top cell type
find_top_cell_type <- function(row) {
  sorted_indices <- order(row, decreasing = TRUE)
  return(names(row)[sorted_indices[1]])
}

cell_type_color2 <- c("coda_tissue_props_acinar" = "purple",
                     "coda_tissue_props_ecm" = "#53C0AD",
                     "coda_tissue_props_fat" = "gray",
                     "coda_tissue_props_islet" = "#F0C61F",
                     "coda_tissue_props_nerve" = "black",
                     "coda_tissue_props_pdac" = "red",
                     "coda_tissue_props_TLS" = "#3333FF",
                     "coda_tissue_props_vasculature" = "white",
                     "coda_tissue_props_whitespace.fiducials" = "#00CCCC")
                     # "coda_tissue_props_acinar" = "gray",
                     # "coda_tissue_props_acinar" = "#3333FF")


cell_type_color <- c("Stroma" = "purple",
                     "GC" = "#53C0AD",
                     "PNI" = "#EFAD1E",
                     "PDAC" = "#F0C61F",
                     "PDAC Mets" = "#FF8000",
                     "LN" = "red",
                     "PANIN" = "#FF6F6F",
                     "TLS" = "#0E2E76",
                     "Acinar" = "#00CCCC",
                     "TBD" = "gray",
                     "Normal Pancreas" = "#3333FF")

arm_color <- c("A" = "black",
               "B" = "red",
               "C" = "blue")

min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

bootCorr <- function(data,
                     replicates = 100,
                     alpha = 0.05
) {
  
  prepArr <- function() {
    z <- array(0,dim = c(ncol(data),
                         ncol(data), replicates ))
    return(z)
  }
  
  setNames <- function(x) {
    colnames(x) <- colnames(data)
    rownames(x) <- colnames(data)
    return(x)
  }
  
  
  nsamples <- nrow(data)
  corrVals <- prepArr()
  
  for (ii in 1:replicates) {
    corrVals[,,ii] <- cor(data %>% sample_n(size = nsamples,
                                            replace = T),
                          method = "pearson")
    
  }
  
  res <- list()
  
  res$mean <- apply(corrVals,
                    c(1,2),
                    mean)
  
  res$ci.lower <- apply(corrVals,
                        c(1,2),
                        quantile,
                        alpha / 2.0)
  res$ci.upper <- apply(corrVals,
                        c(1,2),
                        quantile,
                        1-alpha / 2.0)
  
  res$outside <- ifelse(res$ci.lower * res$ci.upper > 0,
                        1,
                        0)
  
  res <- lapply(res,setNames)
  
  return(res)
  
}


plotVisiumRGB <- function(spe, vars, ...) {
  plt_df <- data.frame(colData(spe), spatialCoords(spe))
  
  if (any(!vars %in% names(plt_df))) {
    stop("One or more variables not found in the data.")
  }
  
  if (length(vars) > 4) {
    stop("A maximum of 4 variables is allowed.")
  }
  
  for (var in vars) {
    plt_df[[var]] <- scales::rescale(plt_df[[var]], to = c(0, 1))
  }
  
  num_vars <- length(vars)
  
  # Initialize channels based on the number of variables:
  # Magenta (R and B channels)
  # Yellow (R and G channels)
  # Green (G channel)
  # Blue (B channel)
  
  if (num_vars >= 1) {
    plt_df$R <- plt_df[[vars[1]]] # Part of Magenta and Yellow
    plt_df$B <- plt_df[[vars[1]]] # Part of Magenta
  }
  
  if (num_vars >= 2) {
    plt_df$R <- plt_df$R + (1 - plt_df$R) * plt_df[[vars[2]]] # Yellow component
    plt_df$G <- plt_df[[vars[2]]] # Green
  }
  
  if (num_vars >= 3) {
    plt_df$G <- plt_df$G + (1 - plt_df$G) * plt_df[[vars[3]]] # Green component
  }
  
  if (num_vars == 4) {
    plt_df$B <- plt_df$B + (1 - plt_df$B) * plt_df[[vars[4]]] # Blue component
  }
  
  spe$RGB <- rgb(plt_df$R, plt_df$G, plt_df$B, maxColorValue = 1)
  plotVisium(spe, fill = "RGB", ...)+scale_fill_identity()
}


prepVisiumRGBNMF <- function(sobj, vars, ...) {
  plt_df <- data.frame(sobj@meta.data)
  
  if (any(!vars %in% names(plt_df))) {
    stop("One or more variables not found in the data.")
  }
  
  if (length(vars) > 4) {
    stop("A maximum of 4 variables is allowed.")
  }
  
  for (var in vars) {
    plt_df[[var]] <- scales::rescale(plt_df[[var]], to = c(0, 1))
  }
  
  num_vars <- length(vars)
  
  # Initialize channels based on the number of variables:
  # Magenta (R and B channels)
  # Yellow (R and G channels)
  # Green (G channel)
  # Blue (B channel)
  
  if (num_vars >= 1) {
    plt_df$R <- plt_df[[vars[1]]] # Part of Magenta and Yellow
    plt_df$B <- plt_df[[vars[1]]] # Part of Magenta
  }
  
  if (num_vars >= 2) {
    plt_df$R <- plt_df$R + (1 - plt_df$R) * plt_df[[vars[2]]] # Yellow component
    plt_df$G <- plt_df[[vars[2]]] # Green
  }
  
  if (num_vars >= 3) {
    plt_df$G <- plt_df$G + (1 - plt_df$G) * plt_df[[vars[3]]] # Green component
  }
  
  if (num_vars == 4) {
    plt_df$B <- plt_df$B + (1 - plt_df$B) * plt_df[[vars[4]]] # Blue component
  }
  
  RGB <- rgb(plt_df$R, plt_df$G, plt_df$B, maxColorValue = 1)
  return(RGB)
}

