# - - - - - - - - - - - - 
# Process raw human Visium data
# Author: Dimitri Sidiropoulos
# - - - - - - - - - - - - 

source("./0_Visium_config_functions.R")

# - - - - - - - - - - - - 
# J1568 SEURAT INIT ####
# - - - - - - - - - - - - 

datadir <- "./rawData/J1568_Visium/count" #SpaceRanger output
samples <- c("J1568_C1","J1568_C2","J1568_C3","J1568_C4_manual_alignment","J1568_C5","J1568_C6","J1568_C7_manual_alignment","J1568_C8","J1568_C9","J1568_C10","J1568_C11","J1568_C12")
samplespaths <- paste0(datadir,"/",samples)

c1 <- Load10X_Spatial(data.dir = samplespaths[1],slice = samples[1],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[1],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE))

# merge all samples
j1568sobj <- merge(c1, y = c(Load10X_Spatial(data.dir = samplespaths[2],slice = samples[2],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[2],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                   Load10X_Spatial(data.dir = samplespaths[3],slice = samples[3],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[3],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                   Load10X_Spatial(data.dir = samplespaths[4],slice = samples[4],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[4],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                   Load10X_Spatial(data.dir = samplespaths[5],slice = samples[5],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[5],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                   Load10X_Spatial(data.dir = samplespaths[6],slice = samples[6],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[6],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                   Load10X_Spatial(data.dir = samplespaths[7],slice = samples[7],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[7],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                   Load10X_Spatial(data.dir = samplespaths[8],slice = samples[8],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[8],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                   Load10X_Spatial(data.dir = samplespaths[9],slice = samples[9],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[9],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                   Load10X_Spatial(data.dir = samplespaths[10],slice = samples[10],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[10],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                   Load10X_Spatial(data.dir = samplespaths[11],slice = samples[11],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[11],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                   Load10X_Spatial(data.dir = samplespaths[12],slice = samples[12],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[12],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE))),
                      add.cell.ids = samples, project = "J1568Visium")

# Fix all images scalefactors to match hires
for (i in 1:12){
  j1568sobj@images[[i]]@scale.factors$lowres<- j1568sobj@images[[i]]@scale.factors$hires #unless you do this the plotting wont work
}

# add panel info to j1568sobj:
j1568sobj@meta.data$panel <- ""
for (panel in samples){
  j1568sobj@meta.data[grepl(panel,colnames(j1568sobj)),]$panel <- panel
}

# add spatial coordinates to metadata
spatCord <- NULL 
for (i in 1:12){
  spatCord <- rbind(spatCord, j1568sobj@images[[i]]@boundaries$centroids@coords)
}
j1568sobj@meta.data <- cbind(j1568sobj@meta.data, spatCord)

# add treatment, tissue metadata
j1568sobj$panel <- str_split_fixed(colnames(j1568sobj), "_", 3)[,2]
j1568sobj$Arm <- ""
j1568sobj$Tissues <- ""
j1568sobj$Response <- ""
j1568sobj$patientIDs <- ""
j1568sobj$LymphoidStructure <- ""

# Create clinical and experimental metadata dataframe denovo
panels <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12")
patientIDs <- c("A11","A11","C4","C4","A3","A15","C10","C1","B14","B10","C5","C7")
Arms <- c("A","A","C","C","A","A","C","C","B","B","C","C")
Tissues <- c("LN","Tumor","LN","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","LN","Tumor","LN")
LymphoidStructure <- c("LN","IntraTLS","LN","PeriTLS","IntraTLS","IntraTLS","IntraTLS","PeriTLS","IntraTLS","LN","PeriTLS","LN")
Response <- c("R","R","R","R","R","NR","R","NR","R","NR","NR","R")
seqslide <- c(1,1,1,1,2,2,2,2,3,3,3,3)
j1568meta <- data.frame(patientIDs = patientIDs, seqslide= seqslide, Panels = panels,Arms = Arms, Tissues = Tissues, Response = Response, LymphoidStructure = LymphoidStructure)
write.csv(j1568meta, "j1568meta.csv")

j1568sobj$Arm <- j1568meta$Arm[match(j1568sobj$panel,j1568meta$Panels)]
j1568sobj$patientIDs <- j1568meta$patientIDs[match(j1568sobj$panel,j1568meta$Panels)]
j1568sobj$Tissues <- j1568meta$Tissues[match(j1568sobj$panel,j1568meta$Panels)]
j1568sobj$Response <- j1568meta$Response[match(j1568sobj$panel,j1568meta$Panels)]
j1568sobj$seqslide <- j1568meta$seqslide[match(j1568sobj$panel,j1568meta$Panels)]
j1568sobj$LymphoidStructure <- j1568meta$LymphoidStructure[match(j1568sobj$panel,j1568meta$Panels)]

## Include CODA, remove non-tissue spots from c4 and c7 
# Remove non-tissue spots from manually aligned samples c4 and c7 
coda_cell_props <- read.table("coda_cell_props_all_samples.txt")
coda_cell_count <- read.table("coda_cell_count_all_samples.txt")
coda_tissue_props <- read.table("coda_tissue_props_all_samples.txt")

colnames(coda_cell_props) <- paste0("coda_cell_props_",colnames(coda_cell_props))
colnames(coda_cell_count) <- paste0("coda_cell_count_",colnames(coda_cell_count))
colnames(coda_tissue_props) <- paste0("coda_tissue_props_",colnames(coda_tissue_props))

data2add <- cbind(coda_cell_props[,7:16],coda_cell_count[,7:16],coda_tissue_props[,7:16])
j1568sobj@meta.data <- cbind(j1568sobj@meta.data, data2add[colnames(j1568sobj),])

# Remove non-tissue spots from manually aligned samples c4 and c7 
plots <- SpatialFeaturePlot(j1568sobj, features = c("coda_tissue_props_whitespace.fiducials"),stroke=NA,image.alpha=1,pt.size.factor=7,combine = FALSE)

c4nontissue <- j1568sobj@meta.data[j1568sobj@meta.data$panel == "C4",] %>% filter(coda_tissue_props_whitespace.fiducials > 90)
c7nontissue <- j1568sobj@meta.data[j1568sobj@meta.data$panel == "C7",] %>% filter(coda_tissue_props_whitespace.fiducials > 90)
nontissue <- c(rownames(c4nontissue),rownames(c7nontissue))
j1568sobj <- subset(j1568sobj, cells = setdiff(colnames(j1568sobj),nontissue))

plots <- SpatialFeaturePlot(j1568sobj, features = c("coda_tissue_props_whitespace.fiducials"),stroke=NA,image.alpha=1,pt.size.factor=7,combine = FALSE)
plots[[4]]
plots[[7]]

## QC and normalize 
j1568sobj$percent.mt <- PercentageFeatureSet(object = j1568sobj, 
                                             pattern = "^MT-") #calculate MT gene expression
VlnPlot(j1568sobj, features = c("nCount_Spatial","nFeature_Spatial","percent.mt"), pt.size = 0)+NoLegend()
j1568sobj <- j1568sobj[,Matrix::colSums(j1568sobj@assays$Spatial@counts) != 0]
j1568sobj <- SCTransform(j1568sobj, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)#normalize, scale.

## Correlate nuclei with UMI counts and genes
colSums(j1568sobj@meta.data[,grep("coda_cell_count",colnames(j1568sobj@meta.data))])

j1568sobj$total_cell_count <- rowSums(j1568sobj@meta.data[,grep("coda_cell_count",colnames(j1568sobj@meta.data))])
total_cell_count <- colSums(j1568sobj@meta.data[,grep("coda_cell_count",colnames(j1568sobj@meta.data))])

par(mar=c(7, 4, 4, 2) + 0.5)
barplot(total_cell_count, main="Total Cell Count", ylab="Count", col="blue", xaxt='n')
# Create a barplot and capture the bar midpoints
midpoints <- barplot(total_cell_count, main="Total Cell Count", ylab="Count", col="blue", xaxt='n')
# Use the midpoints to add rotated x-axis labels
text(x=midpoints, y=-5000, labels=names(total_cell_count), srt=45, adj=1, xpd=TRUE, cex=0.8)


# - - - - - - - - - - - - 
# J0810 SEURAT REINIT ####
# - - - - - - - - - - - - 

datadir <- "./rawData/J0810_Visium/count" #SpaceRanger output
samples <- list.files(path = datadir) #indivdual panels
samplespaths <- paste0(datadir,"/",samples)

s1 <- Load10X_Spatial(data.dir = samplespaths[1],slice = samples[1],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[1],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE))

j0810sobj <- merge(s1, 
              y = c(Load10X_Spatial(data.dir = samplespaths[2],slice = samples[2],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[2],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                    Load10X_Spatial(data.dir = samplespaths[3],slice = samples[3],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[3],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                    Load10X_Spatial(data.dir = samplespaths[4],slice = samples[4],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[4],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                    Load10X_Spatial(data.dir = samplespaths[5],slice = samples[5],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[5],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                    Load10X_Spatial(data.dir = samplespaths[6],slice = samples[6],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[6],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                    Load10X_Spatial(data.dir = samplespaths[7],slice = samples[7],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[7],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                    Load10X_Spatial(data.dir = samplespaths[8],slice = samples[8],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[8],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                    Load10X_Spatial(data.dir = samplespaths[9],slice = samples[9],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[9],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                    Load10X_Spatial(data.dir = samplespaths[10],slice = samples[10],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[10],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                    Load10X_Spatial(data.dir = samplespaths[11],slice = samples[11],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[11],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                    Load10X_Spatial(data.dir = samplespaths[12],slice = samples[12],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[12],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                    Load10X_Spatial(data.dir = samplespaths[13],slice = samples[13],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[13],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                    Load10X_Spatial(data.dir = samplespaths[14],slice = samples[14],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[14],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                    Load10X_Spatial(data.dir = samplespaths[15],slice = samples[15],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[15],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE)),
                    Load10X_Spatial(data.dir = samplespaths[16],slice = samples[16],image = Read10X_Image(image.dir = file.path(paste0(samplespaths[16],"/spatial")), image.name = "tissue_hires_image.png", filter.matrix = TRUE))),
              add.cell.ids = samples, project = "J0810Visium")

# Fix all images scalefactors to match hires
for (i in 1:16){
  j0810sobj@images[[i]]@scale.factors$lowres<- j0810sobj@images[[i]]@scale.factors$hires #unless you do this the plotting wont work
}


## QC, looks good
j0810sobj$percent.mt <- PercentageFeatureSet(object = j0810sobj, pattern = "^MT-") #calculate MT gene expression
VlnPlot(j0810sobj, features = c("nCount_Spatial","nFeature_Spatial","percent.mt"), pt.size = 0)+NoLegend()

## normalize, scale sobj 
j0810sobj <- j0810sobj[,Matrix::colSums(j0810sobj@assays$Spatial@counts) != 0] # filter out spots with zero counts, skip since J0810 has been preprocessed already.
j0810sobj <- SCTransform(j0810sobj, assay = "Spatial", verbose = FALSE) #normalize, scale.
j0810sobj <- j0810sobj %>%  RunPCA(verbose = FALSE)

# - - - - - - - - - - - - 
# SAVE OBJECTs ####
# - - - - - - - - - - - - 

# save j1568sobj
save(j1568sobj,file="./dataObjects/j1568sobj_hiresImages_transformed.rda")
# save j0810sobj
save(j0810sobj,file="./dataObjects/j0180sobj_hiresImages_transformed_projected.rda")

