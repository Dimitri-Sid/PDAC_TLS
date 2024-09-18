# - - - - - - - - - - - - 
# Visium spatial Analysis
# Author: Dimitri Sidiropoulos
# - - - - - - - - - - - - 

source("./0_Visium_config_functions.R")

# load j1568sobj
load("../DataObjects/j1568sobj_hiresImages_transformed_projected.rda")
# load j0810sobj
load("../DataObjects/j0180sobj_hiresImages_transformed_projected.rda")

# - - -  - - - - - - - - -
# Various Spatial  Plots to viz the data ####
# - - - - - - - - - - - - 

p0 <- SpatialDimPlot(j1568sobj,  label = FALSE, group.by = "Anns", label.size = 2,
                     label.box = FALSE,stroke = NA,image.alpha = 1,
                     pt.size.factor = 7, combine = FALSE)
for(i in 1:length(p0)) {
  p0[[i]] <- p0[[i]]  + labs(title = NULL)+ NoLegend() 
}

legend <- cowplot::get_legend( SpatialDimPlot(j1568sobj, stroke = NA, label = FALSE, group.by = "top_cell_type"))

pdf("CODA top cell type.pdf", height = 12, width = 12)
print(cowplot::plot_grid(plotlist = p0))
print(as_ggplot(legend))
dev.off()

p1 <- SpatialFeaturePlot(j1568sobj, features = c("ANXA2"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p2 <- SpatialFeaturePlot(j1568sobj, features = c("coda_tissue_props_TLS"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE)
p3 <- SpatialFeaturePlot(j1568sobj, features = c("coda_cell_count_TLS"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE)
pdf("CODA tls props counts tissue.pdf", height = 12, width = 12)
print(p1)
print(cowplot::plot_grid(plotlist = p2))
print(cowplot::plot_grid(plotlist = p3))
dev.off()

p <- SpatialFeaturePlot(j1568sobj, features = c("Pattern_3"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE)
p <- SpatialFeaturePlot(j1568sobj, features = c("sc_CAF"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE)

DefaultAssay(j1568sobj) <- "module_scores"
Idents(j1568sobj) <- paste(j1568sobj$Arm, j1568sobj$TLSvsGC)
DotPlot(object = j1568sobj, features = c('fDC',"NK","Tcyto","plasma","HEV","KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY","HALLMARK_COMPLEMENT"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# - - -  - - - - - - - - -
# Compute module scores ####
# - - - - - - - - - - - - 

ChemokineScore <- c("CCL21","CCL19","CXCL13","CXCL11","CCL8","CXCL10","CXCL9","CCL2","CCL3","CCL18","CCL5","CCL4") #https://onlinelibrary.wiley.com/doi/10.1002/ijc.32982
iggscore <- c("IGHG1",              "IGHG2",              "IGHG3",              "IGHG4", "IGHM",
              "IGHGP",              'IGKC',              'IGLC1',              'IGLC2',
              'IGLC3',              'IGLC6',              'IGLC7',              "IGHA1") #alpha chain os IGHA#
# IGKV1-5 (LOTS MORE VARIABLE GENES PRESENT)
# FCER2 = CD23, MS4A1 = CD20, CR2 = CD21, important for activation of B cells
Tcyto <- c("CD3E","CD3D","CD8A","GZMB")
NaiveTh <- c("CD3E","CD3D","CCR7","SELL","CD4")
NaiveTc <- c("CD3E","CD3D","CCR7","SELL","CD8A")
ActiveTh <- c("CD3E","CD3D","IL2RA","TFRC", "CD4")
ActiveTc <- c("CD3E","CD3D","IL2RA","TFRC", "CD8A")
Treg <- c("CD3E","CD3D","CD4","FOXP3","CTLA4","PDCD1")
DC <- c("CLEC9A","BTLA","ITGAX","CD1C","LAMP3")
fDC <- c("CR2","CXCL13")
profAPC <- c("HLA-DRA","HLA-DPB1","HLA-DQA1","ITGAX","ITGAM","ITGAE")
tlsMature <- c('RAG1',"AICDA","SELL","CHST4","CLU") #class switching: AICDA
follicularB <- c("CD19","CXCR5","FCER2","TNFRSF13B","PAX5") # TACI = TNFRSF13B
marginalZoneB <- c("CD19","CR2","TNFRSF13B","CD27", "TCF3","SLC22A2","PAX5","EBF1")  
memoryB <- c("CD40","CD80","CD86","CD19","FAS","PTPRJ","TNFRSF13B","POU2AF1","PAX5","CD79A", "MS4A1", "CD27")
regB <- c("CD1D","CD5","CD19","CD24","CD38","CD40","PAX5","TCF3","EBF1", "IL10", "CD79A", "MS4A1")
panB <- c("CD19","MS4A1", "CD79A", "CD22")
naiveB <- c("CD19","MS4A1", "IGHD", "IGHM", "CD79A")
germinalB <- c("CD19","MS4A1","CD79A","BCL6")
plasmablast <- c("CD19","MS4A1","CD27","CD38","MKI67","PWWP3A","CD79A")
plasma <- c("CD19", "CD38", "TNFRSF17", "PWWP3A", "SDC1","MZB1","IRF4","CD27","PRDM1","XBP1","CD38","CXCR4")
NK <- c("NCAM1","KLRD1","NCR1","PRF1","KLRK1")
HEV <- c("CD34", "GLYCAM1", "CCL21", "MADCAM1", "VCAM1", "ICAM1")

human.genes <- msigdbr(species = "Homo sapiens")
genesets.interest <- filter(human.genes, gs_cat == "H" | gs_subcat == "CP:KEGG" | gs_subcat == "CP:HALLMARK"|  gs_subcat == "IMMUNESIGDB" | gs_subcat == "CP:REACTOME")
pathways.interest <- genesets.interest %>% split(x = .$gene_symbol, f = .$gs_name)
# names(pathways.interest[grepl("NK",names(pathways.interest))]) # check which gene sets exist

                      
geneSetsList <- list( Tcyto =Tcyto,
                      NaiveTh= NaiveTh,
                      NaiveTc =NaiveTc,
                      ActiveTh=ActiveTh,
                      ActiveTc =ActiveTc,
                      Treg=Treg,
                      DC =DC,
                      profAPC= profAPC,
                      tlsMature= tlsMature,
                      follicularB=follicularB,
                      marginalZoneB =marginalZoneB,
                      memoryB =memoryB,
                      regB =regB,
                      panB =panB,
                      naiveB=naiveB,
                      germinalB= germinalB,
                      ChemokineScore = ChemokineScore, 
                      iggscore = iggscore,
                      plasmablast= plasmablast,
                      plasma =plasma,
                      fDC = fDC,
                      NK= NK,
                      HEV = HEV,
                      HALLMARK_COMPLEMENT = pathways.interest[grepl("HALLMARK_COMPLEMENT",names(pathways.interest))][[1]],
                      REACTOME_COMPLEMENT_CASCADE= pathways.interest[grepl("REACTOME_COMPLEMENT_CASCADE",names(pathways.interest))][[1]],
                      REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT= pathways.interest[grepl("REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT",names(pathways.interest))][[1]],
                      KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY= pathways.interest[grepl("KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",names(pathways.interest))][[1]],
                      REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION= pathways.interest[grepl("REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",names(pathways.interest))][[1]],
                      KEGG_PANCREATIC_CANCER= pathways.interest[grepl("KEGG_PANCREATIC_CANCER",names(pathways.interest))][[1]],
                      HALLMARK_PANCREAS_BETA_CELLS = pathways.interest[grepl("HALLMARK_PANCREAS_BETA_CELLS",names(pathways.interest))][[1]])

DefaultAssay(j1568sobj) <- "SCT"
j1568sobj <- AddModuleScore(j1568sobj, geneSetsList, name = "module")
colnames(j1568sobj@meta.data)[grepl("module",colnames(j1568sobj@meta.data))] <- names(geneSetsList)
module_scores <- CreateAssayObject(data = t(j1568sobj@meta.data[,names(geneSetsList)]))
j1568sobj[["module_scores"]] <- module_scores #also save as an assay to use in FindSpatiallyVariableFeatures
NMFpatterns <- CreateAssayObject(data = t(j1568sobj@meta.data[,grepl("Pattern",colnames(j1568sobj@meta.data))])) # do the same for cogaps patterns
j1568sobj[["NMFpatterns"]] <- NMFpatterns 

DefaultAssay(j0810sobj) <- "SCT"
j0810sobj <- AddModuleScore(j0810sobj, geneSetsList, name = "module")
colnames(j0810sobj@meta.data)[grepl("module",colnames(j0810sobj@meta.data))] <- names(geneSetsList)
module_scores <- CreateAssayObject(data = t(j0810sobj@meta.data[,names(geneSetsList)]))
j0810sobj[["module_scores"]] <- module_scores #also save as an assay to use in FindSpatiallyVariableFeatures
NMFpatterns <- CreateAssayObject(data = t(j0810sobj@meta.data[,grepl("Pattern",colnames(j0810sobj@meta.data))])) # do the same for cogaps patterns
j0810sobj[["NMFpatterns"]] <- NMFpatterns 

# - - -  - - - - - - - - -
# Spatially Variable Features ####
# - - - - - - - - - - - - 

# J1568
## Identify genes, modules, or patterns that have a spatial enrichment 
j1568sobj <- FindSpatiallyVariableFeatures(j1568sobj, assay = "SCT", features = VariableFeatures(j1568sobj)[1:1000],selection.method = "moransi")
top.features = rownames(dplyr::slice_min(j1568sobj[["SCT"]]@meta.features,moransi.spatially.variable.rank,n = 6))
# repeat but with module scores to see which are spatially enriched
moduleNames <- gsub("_","-",names(geneSetsList))
j1568sobj <- FindSpatiallyVariableFeatures(j1568sobj, assay = "module_scores", features = moduleNames,selection.method = "moransi", slot = "data")
j1568sobj[["module_scores"]]@meta.features
top.modules = rownames(dplyr::slice_min(j1568sobj[["module_scores"]]@meta.features, moransi.spatially.variable.rank,n = 2))
# repeat but with patterns to see which are spatially enriched
patsNames <- paste0("Pattern","-",c(1:15))
j1568sobj <- FindSpatiallyVariableFeatures(j1568sobj, assay = "NMFpatterns", features = patsNames,selection.method = "moransi", slot = "data")
j1568sobj[["NMFpatterns"]]@meta.features
top.patterns = rownames(dplyr::slice_min(j1568sobj[["NMFpatterns"]]@meta.features,moransi.spatially.variable.rank,n = 5))

# J0810
## Identify genes, modules, or patterns that have a spatial enrichment 
j0810sobj <- FindSpatiallyVariableFeatures(j0810sobj, assay = "SCT", features = VariableFeatures(j0810sobj)[1:1000],selection.method = "moransi")
top.features = rownames(dplyr::slice_min(j0810sobj[["SCT"]]@meta.features,moransi.spatially.variable.rank,n = 6))
# repeat but with module scores to see which are spatially enriched
moduleNames <- gsub("_","-",names(geneSetsList))
j0810sobj <- FindSpatiallyVariableFeatures(j0810sobj, assay = "module_scores", features = moduleNames,selection.method = "moransi", slot = "data")
j0810sobj[["module_scores"]]@meta.features
top.modules = rownames(dplyr::slice_min(j0810sobj[["module_scores"]]@meta.features, moransi.spatially.variable.rank,n = 2))
# repeat but with patterns to see which are spatially enriched
patsNames <- paste0("Pattern","-",c(1:10))
j0810sobj <- FindSpatiallyVariableFeatures(j0810sobj, assay = "NMFpatterns", features = patsNames,selection.method = "moransi", slot = "data")
j0810sobj[["NMFpatterns"]]@meta.features
top.patterns = rownames(dplyr::slice_min(j0810sobj[["NMFpatterns"]]@meta.features,moransi.spatially.variable.rank,n = 5))

# - - - - - - - - - - - - 
# Plots patterns on a TLS ####
# - - - - - - - - - - - - 

p0 <- SpatialDimPlot(j1568sobj,  label = FALSE, group.by = "RGB", label.size = 2,
                                          label.box = FALSE,stroke = NA,image.alpha = 1,
                                          pt.size.factor = 3.5, combine = FALSE)
p1 <- SpatialFeaturePlot(j1568sobj, features = c("h_oxphos"),stroke=NA,image.alpha=1,pt.size.factor=2.5,slot ="counts",combine = FALSE)


singletls <- j1568sobj[,j1568sobj$panel == "C9"]
singletls@images <- singletls@images[9]

singletls$RBG <- prepVisiumRGBNMF(singletls, c("Pattern_3","Pattern_9")) # tls, ilset
SpatialDimPlot(singletls, group.by = "RBG",pt.size.factor = 10)+scale_fill_identity()

singletls$RBG <- prepVisiumRGBNMF(singletls, c("Pattern_3","Pattern_4")) # tls pdac1
SpatialDimPlot(singletls, group.by = "RBG",pt.size.factor = 10)+scale_fill_identity()

singletls$RBG <- prepVisiumRGBNMF(singletls, c("Pattern_3","Pattern_7")) #tls pdac2
SpatialDimPlot(singletls, group.by = "RBG",pt.size.factor = 10)+scale_fill_identity()

singletls$RBG <- prepVisiumRGBNMF(singletls, c("Pattern_3","Pattern_4","Pattern_7")) 
SpatialDimPlot(singletls, group.by = "RBG",pt.size.factor = 10)+scale_fill_identity()

singletls$RBG <- prepVisiumRGBNMF(singletls, c("coda_tissue_props_TLS","coda_tissue_props_islet", "coda_tissue_props_pdac"))
SpatialDimPlot(singletls, group.by = "RBG",pt.size.factor = 10)+scale_fill_identity()

                         
p0 <- SpatialDimPlot(j1568sobj,  label = FALSE, group.by = "RGB", label.size = 2,
                     label.box = FALSE,stroke = NA,image.alpha = 1,
                     pt.size.factor = 3.5, combine = FALSE)

print(cowplot::plot_grid(plotlist = trem1))

# - - - - - - - - - - - - 
# Run PatternMarkers ####
# - - - - - - - - - - - - 

# get pattern markers for all the cogaps results
patternMarkerResultsj1568<-patternMarkers(p15,threshold = "cut", lp = NA, axis = 1)
scores <- data.frame(patternMarkerResultsj1568$PatternMarkerScores)
#scores$genes <- rownames(scores)
ranks <- data.frame(patternMarkerResultsj1568$PatternMarkerRanks)
#ranks$genes <- rownames(ranks)
order(ranks$Pattern_9,decreasing=FALSE)[1:30]
j1568 <- patternMarkerResultsj1568$PatternMarkers[[3]]


patternMarkerResultsj0810<-patternMarkers(p10,threshold = "cut", lp = NA, axis = 1)
scores <- data.frame(patternMarkerResultsj0810$PatternMarkerScores)
#scores$genes <- rownames(scores)
ranks <- data.frame(patternMarkerResultsj0810$PatternMarkerRanks)
#ranks$genes <- rownames(ranks)
order(ranks$Pattern_9,decreasing=FALSE)[1:30]
j0810 <- patternMarkerResultsj0810$PatternMarkers[[9]]

# - - - - - - - - - - - - 
# Run SpaceMarkers between TLS and pdac, acinar, and islet pattern for j0810 and j1568  ####
# - - - - - - - - - - - - 

#### j1568 ####
library(SpaceMarkers)
## Set these parameters
# SpInMarkersMode: defaut mode is "residual". You can also set "DE" mode for Differential Expression mode.
SpInMarkersMode = "residual"  
# SpinMarkersRefPattern is the pattern whose "interaction" with every other pattern we want to study. If refPattern is not explicitly assigned, the code assumes Pattern_1 to be refPattern.
SpinMarkersRefPattern = "Pattern_3"
#### Run neighbor identification across all slides:
counts <- j1568sobj@assays$Spatial@counts
images <- j1568sobj@images
panels <- names(images)
patterns <- c(paste0("Pattern_",rep(1:15))) #if you want all to be used
gapsResult <- readRDS("../DataObjects/j1568visium_15p_10kiterations_genome-wide_CogapsParams_backup.rds")
pdat <- j1568sobj@meta.data
spPatterns_allSlides <- data.frame()

for (panel in panels){
  print(panel)
  skip_to_next <- FALSE
  spCoords <- images[[panel]]@coordinates[,c(4,5)]
  colnames(spCoords) <- c("x","y")
  spPatterns <- cbind(spCoords,pdat[rownames(spCoords),patterns]) # Cell type patterns
  cgMat <- gapsResult@featureLoadings %*% t(gapsResult@sampleFactors[rownames(spCoords),])
  #colnames(spPatterns)[-c(1,2)] <- paste0("Pattern_",colnames(spPatterns)[-c(1,2)])
  fullMat <- counts[rownames(cgMat),rownames(spCoords)]
  #optParams <- getSpatialParameters(spPatterns)
  tryCatch({
    SpInMarkers <- getInteractingGenes(data = fullMat, reconstruction = cgMat, spPatterns = spPatterns, refPattern = SpinMarkersRefPattern, mode = SpInMarkersMode, analysis="overlap", minOverlap = 10)
    #SpInMarkers <- getInteractingGenes(data = fullMat, reconstruction = cgMat, spatialPatterns = spPatterns, refPattern = SpinMarkersRefPattern, mode = SpInMarkersMode, minOverlap =10)
    #SpInMarkers$optParams <- optParams
    
    spPatterns <- cbind(spPatterns, SpInMarkers$hotspots)
    #spPatterns <- cbind(spPatterns, SpInMarkers$hotspotRegions)
    spPatterns$panel <- panel
    spPatterns_allSlides <- rbind(spPatterns_allSlides, spPatterns)
    
    panelmarkers <- SpInMarkers$interacting_gene # save all interacting genes for each panel/slide in case we need to retrieve other interactions with pat 3
    names(panelmarkers) <-  unlist(lapply(SpInMarkers$interacting_genes, function(x) colnames(x)[2]))
    saveRDS(panelmarkers, file = paste0("../DataObjects/SpaceMarkers_j1568/",panel, "_interacting_Genes",".RDS"))
  }
  , error = function(e) {skip_to_next <<- TRUE})
}
 
#save(spPatterns_allSlides, file = "../DataObjects/SpaceMarkers_j1568/spPatterns_allSlides.rda")

#### j0810 ####
## Set these parameters
# SpInMarkersMode: defaut mode is "residual". You can also set "DE" mode for Differential Expression mode.
SpInMarkersMode = "residual"  
# SpinMarkersRefPattern is the pattern whose "interaction" with every other pattern we want to study. If refPattern is not explicitly assigned, the code assumes Pattern_1 to be refPattern.
SpinMarkersRefPattern = "Pattern_9" # tls pattern
#### Run neighbor identification across all slides:
counts <- j0810sobj@assays$Spatial@counts
images <- j0810sobj@images
panels <- names(images)
patterns <- c(paste0("Pattern_",rep(1:10))) #if you want all to be used
gapsResult <- readRDS("C:/Users/dsidiro1/OneDrive - Johns Hopkins/FertigLab/Spatial/J0810_trial/dataObjects/cogapsObjects/aggregate/10patterns_aggregated_pdac_visium.rds")
cellnamessplit <- str_split_fixed(rownames(gapsResult@sampleFactors),":",2)
rownames(gapsResult@sampleFactors) <- paste0(cellnamessplit[,2], "_", cellnamessplit[,1])
pdat <- j0810sobj@meta.data
spPatterns_allSlides <- data.frame()

for (panel in panels){
  print(panel)
  skip_to_next <- FALSE
  spCoords <- images[[panel]]@coordinates[,c(4,5)]
  colnames(spCoords) <- c("x","y")
  spPatterns <- cbind(spCoords,pdat[rownames(spCoords),patterns]) # Cell type patterns
  cgMat <- gapsResult@featureLoadings %*% t(gapsResult@sampleFactors[rownames(spCoords),])
  fullMat <- counts[rownames(cgMat),rownames(spCoords)]
  tryCatch({
    SpInMarkers <- getInteractingGenes(data = fullMat, reconstruction = cgMat, spPatterns = spPatterns, refPattern = SpinMarkersRefPattern, mode = SpInMarkersMode, analysis="overlap", minOverlap = 10)
    spPatterns <- cbind(spPatterns, SpInMarkers$hotspots)
    spPatterns$panel <- panel
    spPatterns_allSlides <- rbind(spPatterns_allSlides, spPatterns)
    panelmarkers <- SpInMarkers$interacting_gene # save all interacting genes for each panel/slide in case we need to retrieve other interactions with pat 3
    names(panelmarkers) <-  unlist(lapply(SpInMarkers$interacting_genes, function(x) colnames(x)[2]))
    saveRDS(panelmarkers, file = paste0("../DataObjects/SpaceMarkers_j0810/",panel, "_interacting_Genes",".RDS"))
  }
  , error = function(e) {skip_to_next <<- TRUE})
}

#save(spPatterns_allSlides, file = "../DataObjects/SpaceMarkers_j0810/spPatterns_allSlides.rda")

#### Go over the results:

load("./DataObjects/SpaceMarkers_j1568/spPatterns_allSlides.rda")
load("./DataObjects/SpaceMarkers_j0810/spPatterns_allSlides.rda")

## separated LN from tumor samples and combined both trials into two directories:
# cross-trial-spaceMarkers-no-lymph-nodes
# cross-trial-spaceMarkers-lymph-node-only

# colleect spacemarkers from all slides for a given interaction

SpInMarkersGenes_3_7_allTumors_j1568 <- data.frame() #pdac
SpInMarkersGenes_3_4_allTumors_j1568 <- data.frame() #pdac
SpInMarkersGenes_3_9_allTumors_j1568 <- data.frame() #islet
SpInMarkersGenes_3_6_allTumors_j1568<- data.frame() #acinar
SpInMarkersGenes_3_2_allTumors_j1568 <- data.frame() #ecm/stroma

j1568SMfiles <- list.files("../DataObjects/cross-trial-spaceMarkers-no-lymph-nodes/j1568", full.names= TRUE) # these are lymph node spaceMarkers
for (file in j1568SMfiles){
  
  print(file)
  SpInMarkers <- readRDS(file)
  panel <- str_split_fixed(file, "_", 3)[2]
  
  if ('Pattern_2 x Pattern_3' %in% names(SpInMarkers)){
    SpInMarkersGenes_3_2 <- SpInMarkers[['Pattern_2 x Pattern_3']]
    SpInMarkersGenes_3_2$panel <- panel
    SpInMarkersGenes_3_2_allTumors_j1568 <- rbind(SpInMarkersGenes_3_2_allTumors_j1568, SpInMarkersGenes_3_2)
  }

  if ('Pattern_3 x Pattern_4' %in% names(SpInMarkers)){
    SpInMarkersGenes_3_4 <- SpInMarkers[['Pattern_3 x Pattern_4']]
    SpInMarkersGenes_3_4$panel <- panel
    SpInMarkersGenes_3_4_allTumors_j1568 <- rbind(SpInMarkersGenes_3_4_allTumors_j1568, SpInMarkersGenes_3_4)
  }
  
  if ('Pattern_3 x Pattern_7' %in% names(SpInMarkers)){
    SpInMarkersGenes_3_7 <- SpInMarkers[['Pattern_3 x Pattern_7']]
    SpInMarkersGenes_3_7$panel <- panel
    SpInMarkersGenes_3_7_allTumors_j1568 <- rbind(SpInMarkersGenes_3_7_allTumors_j1568, SpInMarkersGenes_3_7)
  }
  
  if ('Pattern_3 x Pattern_9' %in% names(SpInMarkers)){
    SpInMarkersGenes_3_9 <- SpInMarkers[['Pattern_3 x Pattern_9']]
    SpInMarkersGenes_3_9$panel <- panel
    SpInMarkersGenes_3_9_allTumors_j1568 <- rbind(SpInMarkersGenes_3_9_allTumors_j1568, SpInMarkersGenes_3_9)
  }
  
  if ('Pattern_3 x Pattern_6' %in% names(SpInMarkers)){
    SpInMarkersGenes_3_6 <- SpInMarkers[['Pattern_3 x Pattern_6']]
    SpInMarkersGenes_3_6$panel <- panel
    SpInMarkersGenes_3_6_allTumors_j1568 <- rbind(SpInMarkersGenes_3_6_allTumors_j1568, SpInMarkersGenes_3_6)
  }
}

# Run clusterprofiler:
R.utils::setOption("clusterProfiler.download.method","auto")

enrichGOintersectingGenes(SpInMarkersGenes_3_9_allTumors_j1568)
enrichKEGGintersectingGenes(SpInMarkersGenes_3_7_allTumors_j1568)

## extract Interaction SpaceMarkers to illustrate the interaction component vs each Pattern Alone

SM_p3_p2 <- unique(SpInMarkersGenes_3_2_allTumors_j1568[SpInMarkersGenes_3_2_allTumors_j1568[,2]=="vsBoth",]$Gene)
SM_p3_p7 <- unique(SpInMarkersGenes_3_7_allTumors_j1568[SpInMarkersGenes_3_7_allTumors_j1568[,2]=="vsBoth",]$Gene)
SM_p3_p9 <- unique(SpInMarkersGenes_3_9_allTumors_j1568[SpInMarkersGenes_3_9_allTumors_j1568[,2]=="vsBoth",]$Gene)
SM_p3_p6 <- unique(SpInMarkersGenes_3_6_allTumors_j1568[SpInMarkersGenes_3_6_allTumors_j1568[,2]=="vsBoth",]$Gene)
SM_p3_p4 <- unique(SpInMarkersGenes_3_4_allTumors_j1568[SpInMarkersGenes_3_4_allTumors_j1568[,2]=="vsBoth",]$Gene)

geneSetsList <- list( SM_p3_p2 =SM_p3_p2,
                      SM_p3_p7=SM_p3_p7,
                      SM_p3_p9=SM_p3_p9,
                      SM_p3_p4=SM_p3_p4,
                      SM_p3_p6=SM_p3_p6)
DefaultAssay(j1568sobj) <- "SCT"
j1568sobj <- AddModuleScore(j1568sobj, geneSetsList, name = "module")
colnames(j1568sobj@meta.data)[grepl("module",colnames(j1568sobj@meta.data))] <- names(geneSetsList)
p1 <- SpatialFeaturePlot(j1568sobj, features = c("SM_p3_p2"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p2 <- SpatialFeaturePlot(j1568sobj, features = c("SM_p3_p7"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p3 <- SpatialFeaturePlot(j1568sobj, features = c("SM_p3_p9"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p4 <- SpatialFeaturePlot(j1568sobj, features = c("SM_p3_p6"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p5 <- SpatialFeaturePlot(j1568sobj, features = c("SM_p3_p4"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )

p1[[9]] + p2[[9]] + p3[[9]] + p4[[9]] + p5[[9]]

### Plot hotspots of patterns 2, 5 and 9 on all panels

load("./dataObjects/spPatterns_allSlides_p259.rda")
load("./J0810_trial/dataObjects/SpInMarkers_allSlides_p259.rda")
hotspots <- spPatterns_allSlides[,-c(1,2,3,4,5,9)]

pat2 <- rownames(hotspots[!is.na(hotspots$Pattern_2),])
pat5 <- rownames(hotspots[!is.na(hotspots$Pattern_5),])
pat9 <- rownames(hotspots[!is.na(hotspots$Pattern_9),])
pat2.9 <- intersect(pat2, pat9)
pat5.9 <- intersect(pat9, pat5)
pat2.5 <- intersect(pat2, pat5)

allCoords$hotspots <- "Coldspot"
#allCoords[pat2,]$hotspots <- "p2_hotspot"
allCoords[pat5,]$hotspots <- "p5_hotspot"
allCoords[pat9,]$hotspots <- "p9_hotspot"
# allCoords[pat2.9,]$hotspots <- "p2_p9_intersection"
# allCoords[pat2.5,]$hotspots <- "p2_p5_intersection"
allCoords[pat5.9,]$hotspots <- "p5_p9_intersection"

# need to set diff to make sure assignments are correct
ggplot(allCoords, aes(x = imagecolscaled, y = -imagerowscaled, color = hotspots))+ geom_point(size = 1) + theme_void() + theme(legend.position="right",legend.title = element_blank()) # NEED TO SET DIFF TO INTERSECTED TO MAKE SURE THE SPOTS ARE COLORED CORRECTLY 

### Plot differential expression results on the pattern 9 vs pattern 5 interactions
# extract vs. both for interaction 

vsBoth <- SpInMarkers_allSlides[SpInMarkers_allSlides$`Pattern_5 x Pattern_9` == "vsBoth",]
p5 <- SpInMarkers_allSlides[SpInMarkers_allSlides$`Pattern_5 x Pattern_9` == "vsPattern_5",]
p9 <- SpInMarkers_allSlides[SpInMarkers_allSlides$`Pattern_5 x Pattern_9` == "vsPattern_9",]

# Violin of top genes comparing interacting vs not interacting spots
sobj$hotspots <- allCoords[colnames(sobj),]$hotspots
DotPlot(sobj, features = c("HLA-DPB1","HLA-DRA", "HLA-A", "HLA-B","CD74","MS4A1","INS","IGFBP5", "B2M"), group.by = "hotspots")+RotatedAxis()

vsBothego <- enrichGO(gene         = vsBoth$Gene,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'SYMBOL',
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      maxGSSize = 300)
dotplot(vsBothego, showCategory=20,  font.size = 10)

# need to switch to entrezid for kegg
gene.df <- bitr(vsBoth$Gene, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
R.utils::setOption("clusterProfiler.download.method","auto")
kk <- enrichKEGG(gene = gene.df$ENTREZID)
barplot(kk, showCategory=20)


########################
# TLS composition using projections ####
##########################

j1568sc <- readRDS("C:/Users/dsidiro1/OneDrive - Johns Hopkins/FertigLab/Spatial/J1568_trial/J1568_Visium/dataObjects/20221206_03_J1568singlecell_seurat_providedByJoe.rds")
DimPlot(j1568sc, group.by = "cell_type")

data_matrix <- t(j1568sobj@assays$SeuratPredictions_PDAC_Atlas@data)
top_cell_type <- apply(data_matrix, 1, find_top_cell_type) # regardless of mixture
j1568sobj$top_cell_type_PDAC_Atlas <- top_cell_type

data_matrix <- t(j1568sobj@assays$SeuratPredictions_j1568sc@data)
top_cell_type <- apply(data_matrix, 1, find_top_cell_type) # regardless of mixture
j1568sobj$top_cell_type_j1568sc <- top_cell_type

p <- SpatialDimPlot(j1568sobj,  label = FALSE, group.by = "top_cell_type_PDAC_Atlas", label.size = 2,
                    label.box = FALSE,stroke = NA,image.alpha = 1,
                    pt.size.factor = 5, combine = FALSE)

p <- SpatialDimPlot(j1568sobj,  label = FALSE, group.by = "top_cell_type_j1568sc", label.size = 2,
                    label.box = FALSE,stroke = NA,image.alpha = 1,
                    pt.size.factor = 5, combine = FALSE)


### Compare in an alluvial plot:

dat <- j1568sobj@meta.data[,c("Anns", "top_cell_type_j1568sc")]
dat <- dat[dat$top_cell_type_j1568sc != "CAF",]
# Creating the alluvial plot
pdf("alluvial_tls_gc.pdf")
ggplot(dat[dat$Anns =="TLS" | dat$Anns =="GC",], aes(axis1 = Anns, axis2 = top_cell_type_j1568sc, y = 1)) +
  geom_alluvium(aes(fill = top_cell_type_j1568sc), width = 0) +  # set width to 0 for individual flows
  geom_stratum(width = 0.2) +  # adjust width to separate the strata
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Histology", "Cell Deconvolution"), expand = c(0.05, 0.05)) +
  theme_minimal()
dev.off()


# Calculate proportions
dat <- j1568sobj@meta.data[,c("Anns", "top_cell_type_j1568sc")]
dat <- dat[dat$Anns =="TLS" | dat$Anns =="GC",]
ggplot(dat, aes(x = Anns, fill = top_cell_type_j1568sc)) +
  geom_bar(position="fill") +
  theme_minimal() +
  labs(title = "Stacked Bar Plot of Label2 within Label1",
       x = "Label1 Categories",
       y = "Count",
       fill = "Label2 Categories")    

# Calculate proportions
dat <- j1568sobj@meta.data[,c("Anns","LymphoidStructure", "top_cell_type_j1568sc")]
dat <- dat[dat$Anns =="TLS" | dat$Anns =="GC",]
dat <- dat[dat$top_cell_type_j1568sc !="CAF",]
dat$Categories <- paste(dat$LymphoidStructure, dat$Anns)
ggplot(dat, aes(x = Categories, fill = top_cell_type_j1568sc)) +
  geom_bar(position="fill") +
  theme_minimal() 

ggplot(pdat,aes(y = gene, x=cloupe, fill = cloupe)) +
  geom_violin()+
  theme_classic2() + ylab("Weights") + xlab("Distance")+theme(legend.position="top") + labs(color='Markers')+stat_summary()+
  stat_compare_means(aes(label = ..p.signif..),method="wilcox",paired = FALSE,hide.ns=FALSE, comparisons = list(c("PDAC","Normal Pancreas")))



prop.table(table(dat$Categories, dat$top_cell_type_j1568sc),1)



### repeat in j0810

data_matrix <- t(j0810sobj@assays$SeuratPredictions_PDAC_Atlas@data)
top_cell_type <- apply(data_matrix, 1, find_top_cell_type) # regardless of mixture
j0810sobj$top_cell_type_PDAC_Atlas <- top_cell_type


p <- SpatialDimPlot(j0810sobj,  label = FALSE, group.by = "top_cell_type_PDAC_Atlas", label.size = 2,
                    label.box = FALSE,stroke = NA,image.alpha = 1,
                    pt.size.factor = 7, combine = FALSE)

p[[13]]

# Calculate proportions
dat <- j0810sobj@meta.data[,c("top_cell_type","LN", "top_cell_type_PDAC_Atlas")]
#dat <- dat[dat$Anns =="TLS" | dat$Anns =="GC",]
dat$Categories <- paste(dat$LymphoidStructure, dat$Anns)
ggplot(dat, aes(x = top_cell_type, fill = top_cell_type_PDAC_Atlas)) +
  geom_bar(position="fill") +
  theme_minimal() +  scale_x_discrete(guide = guide_axis(angle = 45))

prop.table(table(dat$Categories, dat$top_cell_type_PDAC_Atlas),1)


########################
# IGG GENE SPATIAL CORRELATIONS ####
##########################

mat1 <- t(as.matrix(j1568sobj@assays$Spatial@data[igscore,]))
mat2 <- t(as.matrix(j0810sobj@assays$Spatial@data[igscore,]))
ggcorrplot::ggcorrplot(cor(rbind(mat1, mat2)),   hc.order = TRUE)

moraniigg <- cbind(
  j0810sobj[["SCT"]]@meta.features[igscore,c(1,2)],
  j1568sobj[["SCT"]]@meta.features[igscore,c(1,2)]
)
moraniigg <- moraniigg[,c(1,3)]
colnames(moraniigg) <- c("J0810", "J1568")
pheatmap(moraniigg, cluster_cols = FALSE, cluster_rows = FALSE)

Idents(j1568sobj) <- j1568sobj$Anns
DefaultAssay(j1568sobj) <- "SCT"
Seurat::DotPlot(j1568sobj, features = igscore)+  scale_x_discrete(guide = guide_axis(angle = 90))


########################
# PLOT PLASMA AND PLASMABLAST CELLS ####
##########################

pdf("plasma_cells_all_samples.pdf")

p <- SpatialFeaturePlot(j1568sobj, features = c("SDC1"),stroke=NA,image.alpha=1,pt.size.factor=5,combine = FALSE, ncol = 4 )
legend <- cowplot::get_legend( p[[1]])
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + labs(title = NULL)+ NoLegend() #
}
print(cowplot::plot_grid(plotlist = p))
print(as_ggplot(legend))

p <- SpatialFeaturePlot(j0810sobj, features = c("SDC1"),stroke=NA,image.alpha=1,pt.size.factor=5,combine = FALSE, ncol = 4 )
legend <- cowplot::get_legend( p[[1]])
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + labs(title = NULL)+ NoLegend() #
}
print(cowplot::plot_grid(plotlist = p))
print(as_ggplot(legend))


p <- SpatialFeaturePlot(j1568sobj, features = c("plasma"),stroke=NA,image.alpha=1,pt.size.factor=5,combine = FALSE, ncol = 4 )
legend <- cowplot::get_legend( p[[1]])
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + labs(title = NULL)+ NoLegend() #
}
print(cowplot::plot_grid(plotlist = p))
print(as_ggplot(legend))

p <- SpatialFeaturePlot(j0810sobj, features = c("plasma"),stroke=NA,image.alpha=1,pt.size.factor=5,combine = FALSE, ncol = 4 )
legend <- cowplot::get_legend( p[[1]])
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + labs(title = NULL)+ NoLegend() #
}
print(cowplot::plot_grid(plotlist = p))
print(as_ggplot(legend))

p <- SpatialFeaturePlot(j1568sobj, features = c("plasmablast"),stroke=NA,image.alpha=1,pt.size.factor=5,combine = FALSE, ncol = 4 )
legend <- cowplot::get_legend( p[[1]])
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + labs(title = NULL)+ NoLegend() #
}
print(cowplot::plot_grid(plotlist = p))
print(as_ggplot(legend))

p <- SpatialFeaturePlot(j0810sobj, features = c("plasmablast"),stroke=NA,image.alpha=1,pt.size.factor=5,combine = FALSE, ncol = 4 )
legend <- cowplot::get_legend( p[[1]])
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + labs(title = NULL)+ NoLegend() #
}
print(cowplot::plot_grid(plotlist = p))
print(as_ggplot(legend))

dev.off()


VlnPlot(j1568sobj, features = "plasma", group.by = "LymphoidStructure",pt.size = 0, ncol = 5)

VlnPlot(j1568sobj, features = "plasma", group.by = "Response",pt.size = 0, ncol = 5)

VlnPlot(j1568sobj, features = "plasma", group.by = "Arm",pt.size = 0, ncol = 5)

VlnPlot(j1568sobj, features = "plasmablast", group.by = "LymphoidStructure",pt.size = 0, ncol = 5)

VlnPlot(j1568sobj, features = "plasmablast", group.by = "Response",pt.size = 0, ncol = 5)

VlnPlot(j1568sobj, features = "plasmablast", group.by = "Arm",pt.size = 0, ncol = 5)


########################
# PLOT CD8A CELLS ####
##########################

pdf("cd8_cells_all_samples.pdf")

p <- SpatialFeaturePlot(j1568sobj, features = c("CD8A"),stroke=NA,image.alpha=1,pt.size.factor=5,combine = FALSE, ncol = 4 )
legend <- cowplot::get_legend( p[[1]])
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + labs(title = NULL)+ NoLegend() #
}
print(cowplot::plot_grid(plotlist = p))
print(as_ggplot(legend))

p <- SpatialFeaturePlot(j0810sobj, features = c("CD8A"),stroke=NA,image.alpha=1,pt.size.factor=5,combine = FALSE, ncol = 4 )
legend <- cowplot::get_legend( p[[1]])
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + labs(title = NULL)+ NoLegend() #
}
print(cowplot::plot_grid(plotlist = p))
print(as_ggplot(legend))


p <- SpatialFeaturePlot(j1568sobj, features = c("Tcyto"),stroke=NA,image.alpha=1,pt.size.factor=5,combine = FALSE, ncol = 4 )
legend <- cowplot::get_legend( p[[1]])
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + labs(title = NULL)+ NoLegend() #
}
print(cowplot::plot_grid(plotlist = p))
print(as_ggplot(legend))

p <- SpatialFeaturePlot(j0810sobj, features = c("Tcyto"),stroke=NA,image.alpha=1,pt.size.factor=5,combine = FALSE, ncol = 4 )
legend <- cowplot::get_legend( p[[1]])
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + labs(title = NULL)+ NoLegend() #
}
print(cowplot::plot_grid(plotlist = p))
print(as_ggplot(legend))

dev.off()


VlnPlot(j1568sobj, features = "Tcyto", group.by = "LymphoidStructure",pt.size = 0, ncol = 5)

VlnPlot(j1568sobj, features = "Tcyto", group.by = "Response",pt.size = 0, ncol = 5)

VlnPlot(j1568sobj, features = "Tcyto", group.by = "Arm",pt.size = 0, ncol = 5)


########################
# IMMUNOGLOBULIN DISTRIBUTION IN J1568 ####
##########################

igscore <- c("IGHG1",              "IGHG2",              "IGHG3",              "IGHG4",
             'IGKC',              'IGLC1',              "IGHM",               'IGLC7',              "IGHA1", "JCHAIN") #alpha chain os IGHA#

geneSetsList <- list( igscore =igscore)
DefaultAssay(j1568sobj) <- "SCT"
j1568sobj <- AddModuleScore(j1568sobj, geneSetsList, name = "module")
colnames(j1568sobj@meta.data)[grepl("module",colnames(j1568sobj@meta.data))] <- names(geneSetsList)

p1 <- SpatialFeaturePlot(j1568sobj, features = c("leukocyte_migration"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p1[[9]]

pdf("igscore_genes_scale.data.pdf")
for (gene in c(igscore, "igscore")){
  print(gene)
  p <- SpatialFeaturePlot(j1568sobj, features = gene,stroke=NA,image.alpha=1,pt.size.factor=7,slot ="scale.data",combine = FALSE, ncol = 4 )
  
  legend <- cowplot::get_legend( p[[1]])
  
  for(i in 1:length(p)) {
    p[[i]] <- p[[i]]  + labs(title = NULL)+ NoLegend() #
  }
  
  print(cowplot::plot_grid(plotlist = p))
  print(as_ggplot(legend))
}
p <- SpatialFeaturePlot(j1568sobj, features = gene,stroke=NA,image.alpha=1,pt.size.factor=0,slot ="counts",combine = FALSE, ncol = 4 )
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + labs(title = NULL)+ NoLegend() #
}
print(cowplot::plot_grid(plotlist = p))
dev.off()

DefaultAssay(j1568sobj) <- "SCT"
VlnPlot(j1568sobj, features = igscore, group.by = "LymphoidStructure",pt.size = 0, ncol = 5)

VlnPlot(j1568sobj, features = igscore, group.by = "Response",pt.size = 0, ncol = 5)

VlnPlot(j1568sobj, features = igscore, group.by = "Arm",pt.size = 0, ncol = 5)


########################
# IMMUNOGLOBULIN DISTRIBUTION IN J0810 ####
##########################

igscore <- c("IGHG1",              "IGHG2",              "IGHG3",              "IGHG4",
             'IGKC',              'IGLC1',              "IGHM",               'IGLC7',              "IGHA1", "JCHAIN") #alpha chain os IGHA#

geneSetsList <- list( igscore =igscore)
DefaultAssay(j0810sobj) <- "SCT"
j0810sobj <- AddModuleScore(j0810sobj, geneSetsList, name = "module")
colnames(j0810sobj@meta.data)[grepl("module",colnames(j0810sobj@meta.data))] <- names(geneSetsList)

p1 <- SpatialFeaturePlot(j0810sobj, features = c("IGHG1"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p1[[9]] + p1[[9]]

pdf("igscore_genes_j0810_counts.pdf")
for (gene in c(igscore, "igscore")){
  print(gene)
  p <- SpatialFeaturePlot(j0810sobj, features = gene,stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
  
  legend <- cowplot::get_legend( p[[1]])
  
  for(i in 1:length(p)) {
    p[[i]] <- p[[i]]  + labs(title = NULL)+ NoLegend() #
  }
  
  print(cowplot::plot_grid(plotlist = p))
  print(as_ggplot(legend))
}
p <- SpatialFeaturePlot(j0810sobj, features = gene,stroke=NA,image.alpha=1,pt.size.factor=0,slot ="counts",combine = FALSE, ncol = 4 )
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + labs(title = NULL)+ NoLegend() #
}
print(cowplot::plot_grid(plotlist = p))
dev.off()

DefaultAssay(j0810sobj) <- "SCT"
VlnPlot(j0810sobj, features = igscore, group.by = "LN",pt.size = 0, ncol = 5)+stat_summary()


# composition 
isotypes <- c("IGHG1",              "IGHG2",              "IGHG3",              "IGHG4",
              "IGHM",           "IGHA1") #alpha chain os IGHA#
mat <- t(as.matrix(j0810sobj@assays$SCT@data[isotypes,]))
df <- data.frame(mat)
df$Anns <- j0810sobj$
  df <- melt(df)

df <- df[df$value >=0,]
df.agg <- aggregate(df, FUN = mean, by = list(df$Anns,df$variable ))

df.agg <- df.agg[df.agg$Group.1 != "TBD",]
df.agg <- df.agg[df.agg$Group.1 != "LN",]
ggplot(df.agg,aes( fill = Group.2, y = value,
                   x = Group.1)) + 
  geom_bar(position="fill", stat="identity",  width=.7) +
  theme_classic()+  scale_x_discrete(guide = guide_axis(angle = 45))


########################
# ECM Remodeling &  IL2/STAT5 in TLS ####
##########################

human.genes <- msigdbr(species = "Homo sapiens")
genesets.interest <- filter(human.genes, gs_cat == "H" | gs_subcat == "CP:KEGG" | gs_subcat == "CP:HALLMARK"|  gs_subcat == "IMMUNESIGDB" | gs_subcat == "CP:REACTOME")
pathways.interest <- genesets.interest %>% split(x = .$gene_symbol, f = .$gs_name)
DefaultAssay(j1568sobj) <- "SCT"
geneSetsList <- list(ecm_degradation = pathways.interest[grepl("REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",names(pathways.interest))][[1]],
                     neutrophil_degran = pathways.interest[grepl("REACTOME_NEUTROPHIL_DEGRANULATION",names(pathways.interest))][[1]],
                     il2_stat5 = pathways.interest[grepl("HALLMARK_IL2_STAT5_SIGNALING",names(pathways.interest))][[1]],
                     leukocyte_migration= pathways.interest[grepl("KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION",names(pathways.interest))][[1]],
                     ecm_organization= pathways.interest[grepl("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",names(pathways.interest))][[1]])
j1568sobj <- AddModuleScore(j1568sobj, geneSetsList, name = "module")
colnames(j1568sobj@meta.data)[grepl("module",colnames(j1568sobj@meta.data))] <- names(geneSetsList)

j1568sobj@meta.data$TREM2 <- as.matrix(j1568sobj@assays$SCT@data["TREM2",])
j1568sobj@meta.data$TREM1 <- as.matrix(j1568sobj@assays$SCT@data["TREM1",])
j1568sobj@meta.data$TREM2_raw <- as.matrix(j1568sobj@assays$SCT@counts["TREM2",])
j1568sobj@meta.data$TREM1_raw <- as.matrix(j1568sobj@assays$SCT@counts["TREM1",])
scpred <-  t(j1568sobj@assays$SeuratPredictions_j1568sc$data)
colnames(scpred) <- paste0("sc_",colnames(scpred))
j1568sobj@meta.data<- cbind(j1568sobj@meta.data,scpred)

pdat <- data.frame(j1568sobj@meta.data)
write.table(pdat,file="j1568pdat.txt")
read.table("j1568pdat.txt")

Idents(j1568sobj) <- j1568sobj$Anns
v1 <- VlnPlot(object = j1568sobj, features = 'ecm_organization', split.by = 'Arm', pt.size = 0.001)
v2 <- VlnPlot(object = j1568sobj, features = 'ecm_degradation', split.by = 'Arm', pt.size = 0.001)
v3 <- VlnPlot(object = j1568sobj, features = 'neutrophil_degran', split.by = 'Arm', pt.size = 0.001)
v4 <- VlnPlot(object = j1568sobj, features = 'TREM2', split.by = 'Arm', pt.size = 0.001)
v5 <- VlnPlot(object = j1568sobj, features = 'sc_TAM', split.by = 'Arm', pt.size = 0.001)

pdf("ecm_violins.pdf", height = 8, width = 12)
print(v1)
print(v2)
print(v3)
print(v4)
print(v5)
dev.off()

p0 <- SpatialFeaturePlot(j1568sobj, features = c("TREM2_raw"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p1 <- SpatialFeaturePlot(j1568sobj, features = c("TREM2"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p2 <- SpatialFeaturePlot(j1568sobj, features = c("ecm_degradation"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p3 <- SpatialFeaturePlot(j1568sobj, features = c("ecm_organization"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p4 <- SpatialFeaturePlot(j1568sobj, features = c("sc_TAM"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p5 <- SpatialFeaturePlot(j1568sobj, features = c("sc_NK"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p6 <- SpatialFeaturePlot(j1568sobj, features = c("sc_cancer"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p7 <- SpatialFeaturePlot(j1568sobj, features = c("sc_CAF"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p8 <- SpatialFeaturePlot(j1568sobj, features = c("sc_Bprolif"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p9 <- SpatialFeaturePlot(j1568sobj, features = c("iggscore"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p10 <- SpatialFeaturePlot(j1568sobj, features = c("sc_Tprolif"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p11 <- SpatialFeaturePlot(j1568sobj, features = c("neutrophil_degran"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )
p12 <- SpatialFeaturePlot(j1568sobj, features = c("coda_tissue_props_TLS"),stroke=NA,image.alpha=1,pt.size.factor=7,slot ="counts",combine = FALSE, ncol = 4 )


pdf("12panelFigs_allSamples.pdf", height = 10, width = 10)
for(num in 1:12){
  print(num)
  print(p1[[num]]+p2[[num]]+p3[[num]]+p4[[num]]+p5[[num]]+p6[[num]]+p7[[num]]+p8[[num]]+p9[[num]]+p10[[num]]+p11[[num]]+p12[[num]])
}
dev.off()


pb <- aggregate.Matrix(j1568sobj@meta.data[,c("TREM2","ecm_degradation","ecm_organization")], 
                       groupings = j1568sobj@meta.data[,c("top_cell_type","Arm")], fun = "mean") 
pheatmap(as.matrix(pb), scale = "column")

features <- c("ecm_degradation","leukocyte_migration" , "iggscore" ,"sc_Bprolif","sc_CAF", "sc_TAM","sc_NK","sc_CD8")
mat <- as.matrix(j1568sobj@meta.data[,features])
ggcorrplot::ggcorrplot(cor(mat),   hc.order = TRUE, p.mat = p.mat)

ct <- corr.test(mat, adjust = "none")
corr <- ct$r
p.mat <- ct$p

# Label significant coefficients with asterisks (*, default) denoting the significance level
ggcorrplot.mixed(cor(mat), upper = "ellipse", lower = "number", p.mat = p.mat,
                 insig = "label_sig", sig.lvl = c(0.05, 0.01, 0.001))





