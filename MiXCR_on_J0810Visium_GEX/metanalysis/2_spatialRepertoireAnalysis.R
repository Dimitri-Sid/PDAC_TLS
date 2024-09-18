
####
#### Run alluvial plots and other available thru scRepertoire #### 
####

load("sobj_transformed_projected_TBCR.rda")
library(Seurat)
library(stringr)
library(dplyr)
library(ape)
library(ggalluvial)
library(scRepertoire) #  scRepertoire: ncborcherding.github.io/vignettes/vignette.html
source("clonotypesAlluvials.R") # functions I modified

pdf("prelimAlluvials.pdf", height = 8, width = 15)

compareClonotypes(sobj@misc$BCR$IGL, 
                  numbers = 50, 
                  samples = names(sobj@misc$BCR$IGL), 
                  cloneCall="CDR3.aa", 
                  graph = "alluvial")+ theme(axis.text.x = element_text(angle = 45, hjust =1)) +ggtitle("top 50 IGL")

compareClonotypes(sobj@misc$BCR$IGH, 
                  numbers = 50, 
                  samples = names(sobj@misc$BCR$IGL), 
                  cloneCall="CDR3.aa", 
                  graph = "alluvial")+ theme(axis.text.x = element_text(angle = 45, hjust =1)) +ggtitle("top 50 IGH")

compareClonotypes(sobj@misc$BCR$IGK, 
                  numbers = 50, 
                  samples = names(sobj@misc$BCR$IGK), 
                  cloneCall="CDR3.aa", 
                  graph = "alluvial")+ theme(axis.text.x = element_text(angle = 45, hjust =1)) +ggtitle("top 50 IGK")

compareClonotypes(sobj@misc$TCR$TRA, 
                  numbers = 30, 
                  samples = names(sobj@misc$TCR$TRA), 
                  cloneCall="CDR3.aa", 
                  graph = "alluvial")+ theme(axis.text.x = element_text(angle = 45, hjust =1)) +ggtitle("TRA")

compareClonotypes(sobj@misc$TCR$TRB, 
                  numbers = 30, 
                  samples = names(sobj@misc$TCR$TRB), 
                  cloneCall="CDR3.aa", 
                  graph = "alluvial")+ theme(axis.text.x = element_text(angle = 45, hjust =1)) +ggtitle("TRB")

### Focus on 116_1

compareClonotypes_modified(sobj@misc$BCR$IGL[9], 
                           numbers = 50, 
                           metadata = "Anns",
                           samples = names(sobj@misc$BCR$IGL)[9], 
                           cloneCall="CDR3.aa",
                           graph = "alluvial", exportTable = FALSE)+ggtitle("IGL top 50 in 116_1")

compareClonotypes_modified(sobj@misc$BCR$IGH[9], 
                           numbers = 50, 
                           metadata = "Anns",
                           samples = names(sobj@misc$BCR$IGH)[9], 
                           cloneCall="CDR3.aa",
                           graph = "alluvial", exportTable = FALSE)+ggtitle("IGH top 50 in 116_1")


compareClonotypes_modified(sobj@misc$BCR$IGK[9], 
                           numbers = 50, 
                           metadata = "Anns",
                           samples = names(sobj@misc$BCR$IGK)[9], 
                           cloneCall="CDR3.aa",
                           graph = "alluvial", exportTable = FALSE)+ggtitle("IGK top 50 in 116_1")

compareClonotypes_modified(sobj@misc$TCR$TRA[2], 
                           numbers = 50, 
                           metadata = "Anns",
                           samples = names(sobj@misc$TCR$TRA)[2], 
                           cloneCall="CDR3.aa",
                           graph = "alluvial", exportTable = FALSE)+ggtitle("TRA in 116_1")

dev.off()


# - - - - - - - - - - - - 
# Spatial Clonotype Analysis ####
# - - - - - - - - - - - - 

source("distMatrix_plot_dendrogram.R") # functions I modified

op <- plot_dendrogram("../runMiXCR/changeo_annotated.rda", 'IGK_7199', 10)

num_colors <- 10

cid_df <- changeo[changeo$locus == "IGH",]
cid_df <- cid_df[cid_df$sample == "116_1",]

# Take counts and cdr3_aa to collapse rows
cid_df$cdr3_aa_ann <- paste0(cid_df$cdr1_aa, "_", cid_df$spatial_annotations)
cid_subset <- cid_df[, c('cdr3_aa', 'count' )]
cid_collapsed <- aggregate(.~cdr3_aa,data=cid_subset, FUN=sum )

d <- data.frame(matrix(ncol=(nrow(cid_collapsed) + 1), nrow=nrow(cid_collapsed)))
cols <- cid_collapsed$cdr3_aa
colnames(d) <- c('cdr3_aa', cols)
d$cdr3_aa <- cid_collapsed$cdr3_aa
for (i in 2:ncol(d)) {
  d[, i] <- stringdist(cols[i - 1], d[, 1], method='lv')
}

dist <- stats::as.dist(d[, 2:ncol(d)])

# Create distance matrix using hierarchical clustering
hc <- as.dendrogram(hclust(dist))

# Perform labeling of counts, branches, etc.
hc <- color_branches(hc)
hc <- hc %>% set('branches_lwd', 4)

labels <- hc %>% labels
rownames(cid_collapsed) <- cid_collapsed$cdr3_aa
counts = cid_collapsed[labels,]$count

label_color<- color.scale(counts, extremes = c("#0000ff", "#ff0000"))

# which are in TLS?
cid_subset2 <- cid_df[, c('cdr3_aa', 'count',"spatial_annotations" )]
cid_collapsed2 <- aggregate(.~cdr3_aa+spatial_annotations,data=cid_subset2, FUN=sum )
tlslabels <- cid_collapsed2[cid_collapsed2$spatial_annotations == "TLS",]$cdr3_aa
pdaclabels <- cid_collapsed2[cid_collapsed2$spatial_annotations == "PDAC",]$cdr3_aa

pdf("dendrogram_116_1_IGH_TLSred.pdf")
hc <- hc %>% set('leaves_pch', 19) %>% set("labels_colors") %>% set('leaves_cex', 1) %>% set('labels_cex', 0.5) %>% set("branches_col", 1) %>% set('leaves_col', label_color) %>% set("by_labels_branches_col", value = tlslabels)
circlize_dendrogram(hc, dend_track_height = 0.4, labels_track_height = 0.5)
dev.off()

pdf("dendrogram_116_1_IGH_PDACred.pdf")
hc <- hc %>% set('leaves_pch', 19) %>% set("labels_colors") %>% set('leaves_cex', 1) %>% set('labels_cex', 0.5) %>% set("branches_col", 1) %>% set('leaves_col', label_color) %>% set("by_labels_branches_col", value = pdaclabels)
circlize_dendrogram(hc, dend_track_height = 0.4, labels_track_height = 0.5)
dev.off()

#Plot a sample clone; for example clone21, a TCR with CDR3 AA sequence CASSFRPGQPDNEQFF, TRBV7-8, and TRBJ2-1
#In scRNA-seq data, this clone is 82% in clusters B/C (less-exhausted)
SpatialPlot(brain16, features = "clone21") &
  scale_fill_gradient2(low = "white", mid = "blue", high = "red", midpoint = 0.5) & 
  scale_alpha_continuous(range=c(0,1))

## GET MORANI:
d <- cid_collapsed[6,1] # extract barcode 1 by 1
cid_df[cid_df$cdr3_aa == d,]

# Assuming this is your matrix of coordinates

# Calculate the distance matrix
dist_matrix <- as.matrix(dist) ## Assuming this is your matrix of coordinates

# Calculate Moran's I using the ape package
moran_i <- Moran.I(counts, dist_matrix) ## i think i need the dist_matrix to be the entire capture area and the counts to be the total for that clone across all spots

# Print the result
print(moran_i)

## plot top CDR3 sequences 
cid_df <- changeo[changeo$locus == "IGH",]
cid_df <- cid_df[cid_df$sample == "116_1",]

j0810sobj_116_1 <- j0810sobj[,j0810sobj$panel == "VI_116_1"]
clones <- cid_df$cdr3_aa

for (clone in clones){
  barcodes <- cid_df[cid_df$cdr3_aa == clone,]$visium_barcode
  j0810sobj_116_1@meta.data[,clone] <- 0
  for (barcode in barcodes){
    j0810sobj_116_1@meta.data[barcode,clone] <- j0810sobj_116_1@meta.data[barcode,clone] + sum(cid_df[ cid_df$visium_barcode == barcode & cid_df$cdr3_aa == clone,]$count)
  }
}

j0810sobj_116_1$cloneCountSums <- rowSums(j0810sobj_116_1@meta.data[,c(clones)])

mat <-j0810sobj_116_1@meta.data[,unique(clones)]
most_expressed_clone <- apply(mat, 1, function(x) colnames(mat)[which.max(x)])
j0810sobj_116_1$most_expressed_clone <- most_expressed_clone

nokeeps <- rowSums(j0810sobj_116_1@meta.data[,c(unique(clones))])==0
keeps <- rowSums(j0810sobj_116_1@meta.data[,c(unique(clones))])>0
j0810sobj_116_1@meta.data[nokeeps,]$most_expressed_clone <- "N/A" # filter out those with no clone counts

SpatialFeaturePlot(j0810sobj_116_1,features =c("cloneCountSums"),stroke=NA,image.alpha=1,pt.size.factor=6,combine = FALSE)[[13]]

SpatialDimPlot(j0810sobj_116_1[,j0810sobj_116_1$most_expressed_clone != "N/A"],group.by=c("most_expressed_clone"),stroke=NA,image.alpha=1,pt.size.factor=15,combine = FALSE)[[13]]

ggplot(j0810sobj_116_1@meta.data,aes( fill = top_cell_type, y = cloneCountSums,
                                      x = top_cell_type)) + 
  geom_bar(stat="identity",  width=.7) +
  theme_classic()+  scale_x_discrete(guide = guide_axis(angle = 45))

kbl(table(j0810sobj_116_1$most_expressed_clone))

mat <- as.matrix(j0810sobj_116_1@meta.data[keeps,unique(clones)])
anns <- j0810sobj_116_1@meta.data[keeps,]$top_cell_type
names(anns) <- keeps

annotation = data.frame(Annotations = anns)
rownames(annotation) = rownames(j0810sobj_116_1@meta.data[keeps,])

pheatmap(scale(mat),show_rownames= FALSE, annotation_row = annotation)

col_clust <- hclust(dist) #calculated in dendrogram section above

# Use this order in pheatmap without performing additional clustering
pheatmap(scale(mat), 
         cluster_cols = col_clust, show_rownames= FALSE, annotation_row = annotation,fontsize_col = 8)

# Label significant coefficients with asterisks (*, default) denoting the significance level
pdat <- data.frame(j0810sobj_116_1@meta.data)
mat <- pdat[,c("cloneCountSums","Pattern_9","iggscore", "Pattern_1","Pattern_2","Pattern_3","Pattern_5","Pattern_6","Pattern_4","Pattern_7","Pattern_8","Pattern_10")]
ggcorrplot.mixed(cor(mat), upper = "ellipse", lower = "number")















