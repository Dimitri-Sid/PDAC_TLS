####
#### Load MIXCR results, create prelim immunarch report and embded clonotypes into Seurat object #### 
####

library(psych)
library(corrplot)
library(Seurat)
load("./runMiXCR/22_06_01_pdac-visium-v1.rda")
load("../j0180sobj_hiresImages_transformed_projected.rda") #seurat object
sobj$log10GenesPerUMI <- log10(sobj$nFeature_Spatial) / log10(sobj$nCount_Spatial)

View(res)
dim(res)
length(unique(res$barcode))
#2234

table(res$sample)
table(res$locus)

## match batcodes
res$visium_barcode <- paste0("VI_",res$sample,"_",res$barcode, "-1")

## add res data as separate assays for each locus: "IGL" "IGK" "IGH" "TRB" "TRA"
## add metadata to res
res$Anns <- ""
barcodes <- unique(res$visium_barcode)
for (barcode in barcodes){
  res[res$visium_barcode == barcode,]$Anns <- sobj@meta.data[barcode,]$Anns
}

#### Create a data structure that separates all chains by patient into separate tibbles and metadata
res$CDR3.aa <- res$cdr3_aa
res$CDR3.nt <- res$cdr3
res$Clones <- res$cloneCount
repTCR <- res[res$locus == "TRB" | res$locus == "TRA",] #all tcr no bcr
repBCR <- res[res$locus != "TRB" & res$locus != "TRA",] #all bcr no tcr
repIGH <- res[res$locus == "IGH",]
repIGK <- res[res$locus == "IGK",]
repIGL <- res[res$locus == "IGL",]
visiumTCR <- list()
visiumBCR <- list()

## Get unique clones per spot
repBCR$uniqueClones <- 0
barcodes <- unique(repBCR$visium_barcode)
for (barcode in barcodes){
  repBCR[repBCR$visium_barcode == barcode,]$uniqueClones <- dim(repBCR[repBCR$visium_barcode == barcode,])[1]
}

repIGL$uniqueClones <- 0
barcodes <- unique(repIGL$visium_barcode)
for (barcode in barcodes){
  repIGL[repIGL$visium_barcode == barcode,]$uniqueClones <- dim(repIGL[repIGL$visium_barcode == barcode,])[1]
}

repIGH$uniqueClones <- 0
barcodes <- unique(repIGH$visium_barcode)
for (barcode in barcodes){
  repIGH[repIGH$visium_barcode == barcode,]$uniqueClones <- dim(repIGH[repIGH$visium_barcode == barcode,])[1]
}

repIGK$uniqueClones <- 0
barcodes <- unique(repIGK$visium_barcode)
for (barcode in barcodes){
  repIGK[repIGK$visium_barcode == barcode,]$uniqueClones <- dim(repIGK[repIGK$visium_barcode == barcode,])[1]
}

#### BCR clonotypes
library(immunarch); data("immdata") #load package and example data
visiumBCR$allBCR <- list()
visiumBCR$IGH <- list()
visiumBCR$IGL <- list()
visiumBCR$IGK <- list()

for (sample in unique(repBCR$sample)){
  visiumBCR$allBCR[[sample]] <- repBCR[repBCR$sample == sample,] 
}

for (sample in unique(repIGH$sample)){
  visiumBCR$IGH[[sample]] <- repIGH[repIGH$sample == sample,] 
}

for (sample in unique(repIGL$sample)){
  visiumBCR$IGL[[sample]] <- repIGL[repIGL$sample == sample,] 
}

for (sample in unique(repIGK$sample)){
  visiumBCR$IGK[[sample]] <- repIGK[repIGK$sample == sample,] 
}

visiumBCR$meta <- tibble(Sample = unique(repBCR$sample))
visiumBCR$meta$treatment <- sobj@meta.data$Arm[match(paste0("VI_",visiumBCR$meta$Sample),sobj@meta.data$panel)]

## run stats from immunoarch for each chain:
pdf("immunarch_output.pdf", height = 20, width = 20)
p1 <- repOverlap(visiumBCR$IGL)%>% vis() + ggtitle("IGL") 
#geneUsage(visiumBCR$IGL[[1]]) %>% vis() # does not work unless we load directly from mixcr using repLoad
p2 <- repDiversity(visiumBCR$IGL) %>% vis(.by = "treatment", .meta = visiumBCR$meta)   # by annotation
p3 <- repDiversity(visiumBCR$IGL) %>% vis()
p4 <- repDiversity(visiumBCR$IGL,.method = "hill") %>% vis(.by = "treatment", .meta = visiumBCR$meta)    # by annotation
p5 <- repDiversity(visiumBCR$IGL, .method = "hill") %>% vis() 
p6 <- repExplore(visiumBCR$IGL, "lens") %>% vis()   # distribution of CDR3, can be indicative of a clone family
p7 <- repClonality(visiumBCR$IGL, "homeo") %>% vis()  # Visualise the relative abundance of clonotypes
p8 <- repClonality(visiumBCR$IGL, "clonal.prop") %>% vis()  # Visualise the relative abundance of clonotypes
print(p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8)

p1 <- repOverlap(visiumBCR$IGH)%>% vis() + ggtitle("IGH") 
#geneUsage(visiumBCR$IGH[[1]]) %>% vis() # does not work unless we load directly from mixcr using repLoad
p2 <- repDiversity(visiumBCR$IGH) %>% vis(.by = "treatment", .meta = visiumBCR$meta)   # by annotation
p3 <- repDiversity(visiumBCR$IGH) %>% vis()
p4 <- repDiversity(visiumBCR$IGH,.method = "hill") %>% vis(.by = "treatment", .meta = visiumBCR$meta)    # by annotation
p5 <- repDiversity(visiumBCR$IGH, .method = "hill") %>% vis() 
p6 <- repExplore(visiumBCR$IGH, "lens") %>% vis()   # distribution of CDR3, can be indicative of a clone family
p7 <- repClonality(visiumBCR$IGH, "homeo") %>% vis()  # Visualise the relative abundance of clonotypes
p8 <- repClonality(visiumBCR$IGH, "clonal.prop") %>% vis()  # Visualise the relative abundance of clonotypes
print(p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8)

p1 <- repOverlap(visiumBCR$IGK)%>% vis() + ggtitle("IGK") 
#geneUsage(visiumBCR$IGK[[1]]) %>% vis() # does not work unless we load directly from mixcr using repLoad
p2 <- repDiversity(visiumBCR$IGK) %>% vis(.by = "treatment", .meta = visiumBCR$meta)   # by annotation
p3 <- repDiversity(visiumBCR$IGK) %>% vis()
p4 <- repDiversity(visiumBCR$IGK,.method = "hill") %>% vis(.by = "treatment", .meta = visiumBCR$meta)    # by annotation
p5 <- repDiversity(visiumBCR$IGK, .method = "hill") %>% vis() 
p6 <- repExplore(visiumBCR$IGK, "lens") %>% vis()   # distribution of CDR3, can be indicative of a clone family
p7 <- repClonality(visiumBCR$IGK, "homeo") %>% vis()  # Visualise the relative abundance of clonotypes
p8 <- repClonality(visiumBCR$IGK, "clonal.prop") %>% vis()  # Visualise the relative abundance of clonotypes
print(p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8)

dev.off()

pdf("other_results.pdf", height = 8, width = 15)
## Clone count by annotations
p1 <- ggplot(repIGH, aes(x = Anns, y = cloneCount)) + geom_boxplot() + geom_jitter() + theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ggtitle("IGH UMI Clone Count")
p2 <- ggplot(repIGL, aes(x = Anns, y = cloneCount)) + geom_boxplot() + geom_jitter() + theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ggtitle("IGL UMI Clone Count")
p3 <- ggplot(repIGK, aes(x = Anns, y = cloneCount)) + geom_boxplot() + geom_jitter() + theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ggtitle("IGK UMI Clone Count")
p4 <- ggplot(repIGH, aes(x = Anns, y = uniqueClones)) + geom_boxplot() + geom_jitter() + theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ggtitle("IGH Unique Clone Count")
p5 <- ggplot(repIGL, aes(x = Anns, y = uniqueClones)) + geom_boxplot() + geom_jitter() + theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ggtitle("IGL Unique Clone Count")
p6 <- ggplot(repIGK, aes(x = Anns, y = uniqueClones)) + geom_boxplot() + geom_jitter() + theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ggtitle("IGK Unique Clone Count")
print(p1+p2+p3+p4+p5+p6)

# clonotype characterization:
top(visiumBCR$allBCR[["113_1"]]) #top clonotypes
coding(visiumBCR$allBCR[[1]])
noncoding(visiumBCR$allBCR[[1]]) # alex filtered out noncoding hence there are none
nrow(inframes(visiumBCR$allBCR[[1]]))
nrow(outofframes(visiumBCR$allBCR[[1]]))
## add BCR to seurat misc slot
sobj@misc$BCR <- visiumBCR

#### Aggregate clone UMI counts, and append unique clones; paste to sobj metadata
aggrIGHumi <- aggregate(cloneCount ~ visium_barcode, data = repIGH, FUN = sum)
aggrIGLumi <- aggregate(cloneCount ~ visium_barcode, data = repIGL, FUN = sum)
aggrIGKumi <- aggregate(cloneCount ~ visium_barcode, data = repIGK, FUN = sum)

aggrIGHuniq <- aggregate(uniqueClones ~ visium_barcode, data = repIGH, FUN = mean)
aggrIGLuniq <- aggregate(uniqueClones ~ visium_barcode, data = repIGL, FUN = mean)
aggrIGKuniq <- aggregate(uniqueClones ~ visium_barcode, data = repIGK, FUN = mean)

sobj@meta.data$aggrIGHumi <- 0
sobj@meta.data$aggrIGLumi <- 0
sobj@meta.data$aggrIGKumi <- 0
sobj@meta.data$aggrIGHuniq <- 0
sobj@meta.data$aggrIGLuniq <- 0
sobj@meta.data$aggrIGKuniq <- 0

sobj@meta.data[aggrIGHumi$visium_barcode,]$aggrIGHumi <- as.numeric(aggrIGHumi$cloneCount)
sobj@meta.data[aggrIGLumi$visium_barcode,]$aggrIGLumi <- as.numeric(aggrIGLumi$cloneCount)
sobj@meta.data[aggrIGKumi$visium_barcode,]$aggrIGKumi <- as.numeric(aggrIGKumi$cloneCount)

sobj@meta.data[aggrIGHuniq$visium_barcode,]$aggrIGHuniq <- as.numeric(aggrIGHuniq$uniqueClones)
sobj@meta.data[aggrIGLuniq$visium_barcode,]$aggrIGLuniq <- as.numeric(aggrIGLuniq$uniqueClones)
sobj@meta.data[aggrIGKuniq$visium_barcode,]$aggrIGKuniq <- as.numeric(aggrIGKuniq$uniqueClones)

## Plot clone counts by chain on slides
p1 <- SpatialFeaturePlot(sobj, features = c("aggrIGHumi"),stroke=NA,image.alpha=2,pt.size.factor=1.5,combine = FALSE)
# for(i in 1:length(p)) {
#   p[[i]] <- p[[i]]  + labs(title = NULL)#+ NoLegend()
# }
p2 <- SpatialFeaturePlot(sobj, features = c("aggrIGLumi"),stroke=NA,image.alpha=2,pt.size.factor=1.5,combine = FALSE)
p3 <- SpatialFeaturePlot(sobj, features = c("aggrIGKumi"),stroke=NA,image.alpha=2,pt.size.factor=1.5,combine = FALSE)
p4 <- SpatialFeaturePlot(sobj, features = c("aggrIGHuniq"),stroke=NA,image.alpha=2,pt.size.factor=1.5,combine = FALSE)
p5 <- SpatialFeaturePlot(sobj, features = c("aggrIGLuniq"),stroke=NA,image.alpha=2,pt.size.factor=1.5,combine = FALSE)
p6 <- SpatialFeaturePlot(sobj, features = c("aggrIGKuniq"),stroke=NA,image.alpha=2,pt.size.factor=1.5,combine = FALSE)

## plot top CDR3 sequences like the renal spatial BCR paper
top116_1 <- top(visiumBCR$allBCR[["116_1"]]) #top clonotypes from all chains
sobj@meta.data$top116_1 <- ""
sobj@meta.data[top116_1$visium_barcode,]$top116_1 <- top116_1$CDR3.aa
p <- SpatialDimPlot(sobj,group.by  = "top116_1",stroke=NA,image.alpha=2,pt.size.factor=1.5,combine = FALSE)

print(p1[[13]] + p2[[13]] + p3[[13]] + p4[[13]] + p5[[13]] + p6[[13]])
print(p[[13]]+scale_fill_manual(values = c("NA", "blue", "green","red","orange","black","darkblue","gray","yellow","purple")))


## plot correlations of clone counts
cor.test(sobj@meta.data$aggrIGH,
         sobj@meta.data$aggrIGK,
         sobj@meta.data$aggrIGL)

mat <- as.matrix(sobj@meta.data[,c("aggrIGLumi","aggrIGHumi","aggrIGKumi","aggrIGLuniq","aggrIGHuniq","aggrIGKuniq","nCount_Spatial","nFeature_Spatial","Size_Factor","log10GenesPerUMI")])
corrs <- cor(mat)
cor_test_mat <- corr.test(mat,method = "spearman")
corr.p(cor_test_mat$r,cor_test_mat$n,adjust="holm",alpha=.05,minlength=5,ci=TRUE)
corrplot(corrs, type = "upper", order = "hclust",diag = FALSE, tl.col = "black", tl.srt = 45)
corrplot(corrs, type = "full", order = "hclust",diag = TRUE, tl.col = "black", tl.srt = 45)
dev.off()

#### At this point I have a seurat object with all chains in sobj@misc$BCR and UMI/Clone counts in metadata (cloneCount,uniqueClones --> aggrIG*umi, aggrIG*uniq)
# now add TCR and save
sobj@misc$allTCR <- list()
sobj@misc$TCR$compiledTCR <- repTCR #not sep by sample
sobj@misc$TCR$allTCR <- list()
sobj@misc$TCR$TRA <- list()
sobj@misc$TCR$TRB <- list()

repTRA <- repTCR[repTCR$locus == "TRA",] # separate by chain
repTRB <- repTCR[repTCR$locus == "TRB",]

for (sample in unique(repTCR$sample)){
  sobj@misc$TCR$allTCR[[sample]] <- repTCR[repTCR$sample == sample,] 
}

for (sample in unique(repTRA$sample)){
  sobj@misc$TCR$TRA[[sample]] <- repTRA[repTRA$sample == sample,] 
}

for (sample in unique(repTRB$sample)){
  sobj@misc$TCR$TRB[[sample]] <- repTRB[repTRB$sample == sample,] 
}

save(sobj, file = "sobj_transformed_projected_TBCR.rda")
