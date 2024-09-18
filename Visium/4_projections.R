# - - - - - - - - - - - - 
# Visium projections & integrations with single cell data
# Author: Dimitri Sidiropoulos
# - - - - - - - - - - - - 

# - - - - - - - - - - - - 
# Integrate PDAC ATLAS single cell with J1568 Visium ####
# - - - - - - - - - - - - 

load("./dataObjects/PDACAtlas_downsampled_seuratObject.rda")
DimPlot(pdacatlas_sobj, group.by = "anns")

anchors <- FindTransferAnchors(reference = pdacatlas_sobj, query = j1568sobj, normalization.method = "SCT")

Spredictions.assay <- TransferData(anchorset = anchors, refdata = pdacatlas_sobj$anns, prediction.assay = TRUE, weight.reduction = j1568sobj[["pca"]],   dims =1:30 )
j1568sobj[["SeuratPredictions_PDAC_Atlas"]] <- Spredictions.assay
SpatialFeaturePlot(j1568sobj, features = c("T cell", "B cell", "Epithelial-cancer","Treg cell"), ncol = 2, alpha = c(0.1, 1), pt.size.factor = 10)

# - - - - - - - - - - - - 
# Integrate PDAC ATLAS single cell with J0810 Visium ####
# - - - - - - - - - - - - 

load("./dataObjects/PDACAtlas_downsampled_seuratObject.rda")
DimPlot(pdacatlas_sobj, group.by = "anns")

anchors <- FindTransferAnchors(reference = pdacatlas_sobj, query = j0810sobj, normalization.method = "SCT")

Spredictions.assay <- TransferData(anchorset = anchors, refdata = pdacatlas_sobj$anns, prediction.assay = TRUE, weight.reduction = j0810sobj[["pca"]],   dims =1:30 )
j0810sobj[["SeuratPredictions_PDAC_Atlas"]] <- Spredictions.assay
DefaultAssay(j0810sobj) <- "SeuratPredictions_PDAC_Atlas"
SpatialFeaturePlot(j0810sobj, features = c("T cell", "B cell", "Epithelial-cancer","Treg cell"), ncol = 2, alpha = c(0.1, 1), pt.size.factor = 10)

# - - - - - - - - - - - - 
# Create annotations based on top CODA classifications and top single cell projections for J1568 ####
# - - - - - - - - - - - - 

data_matrix <- j1568sobj@meta.data[,grepl("coda_tissue_props", colnames(j1568sobj@meta.data))]
top_cell_types <- apply(data_matrix, 1, find_top_cell_types) # this is at 90% threshold
top_cell_type <- apply(data_matrix, 1, find_top_cell_type) # regardless of mixture

j1568sobj$top_cell_types <- top_cell_types
plots <- SpatialDimPlot(j1568sobj,  group.by="top_cell_types",stroke=NA,image.alpha=1,pt.size.factor=7,combine = FALSE)
plots[[3]]
j1568sobj$top_cell_type <- top_cell_type

# - - - - - - - - - - - - 
# Incorporate CoGAPS ####
# - - - - - - - - - - - - 

p15 <- readRDS("./dataObjects/j1568visium_15p.rds")
j1568sobj@meta.data <- cbind(j1568sobj@meta.data, p15@sampleFactors)

NMFpatterns <- CreateAssayObject(data = t(j1568sobj@meta.data[,grepl("Pattern",colnames(j1568sobj@meta.data))])) # also add as a slot
j1568sobj[["NMFpatterns"]] <- NMFpatterns 

# pseudobulk and see which patterns overlap with which cell types or layers
pb <- aggregate.Matrix(j1568sobj@meta.data[,paste0("Pattern_",rep(1:numpat))], 
                       groupings = j1568sobj$Layer, fun = "mean") %>% data.frame()
pheatmap(scale(pb))

# plot patterns of interest on spatial plots for demo:
plots <- SpatialFeaturePlot(j1568sobj, features = c("Pattern_1"),stroke=NA,image.alpha=1,pt.size.factor=2.5,slot ="counts",combine = FALSE)
plots[[1]]


p10 <- readRDS("./dataObjects/10patterns_aggregated_pdac_visium.rds")
j0810sobj@meta.data <- cbind(j0810sobj@meta.data, p15@sampleFactors)

NMFpatterns <- CreateAssayObject(data = t(j0810sobj@meta.data[,grepl("Pattern",colnames(j0810sobj@meta.data))])) # also add as a slot
j0810sobj[["NMFpatterns"]] <- NMFpatterns 

# - - - - - - - - - - - - 
# ProjectR J0810 CoGAPS ONTO J1568 ####
# - - - - - - - - - - - - 

p10 <- readRDS("./dataObjects/10patterns_aggregated_pdac_visium.rds")
spatial_matrix <- as.matrix(j1568sobj@assays[["Spatial"]]@counts)
projection <- projectR(data = spatial_matrix, loadings = p10, full=FALSE)
new.assay <- CreateAssayObject(data = projection)
j1568sobj[["j0810_CoGAPS_projections"]] <-  new.assay

# - - - - - - - - - - - - 
# SAVE OBJECTs ####
# - - - - - - - - - - - - 

# save j1568sobj
save(j1568sobj,file="./dataObjects/j1568sobj_hiresImages_transformed_projected.rda")
# save j0810sobj
save(j0810sobj,file="./dataObjects/j0180sobj_hiresImages_transformed_projected.rda")
