# - - - - - - - - - - - - 
# Example preparation for CoGAPS Run on an HPC cluster 
# Author: Dimitri Sidiropoulos
# - - - - - - - - - - - - 

#directory to where the data files are
data_dir<- getwd()
#load in the sample names
samples<- list.files(data_dir)
#set the working directory

# Create separate matrices for each sample
for (i in 1:length(samples)){
 # Set up the path names for the samples
 barcode.path <- paste0(data_dir,"/",samples[i],'/filtered_feature_bc_matrix/barcodes.tsv.gz')
 features.path <- paste0(data_dir,"/", samples[i],"/filtered_feature_bc_matrix/features.tsv.gz")
 matrix.path <- paste0(data_dir,"/", samples[i],"/filtered_feature_bc_matrix/matrix.mtx.gz")
 #read in the matrix
 mat <- readMM(file = matrix.path)
 #read in the gene names
 feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
 #read in the spot ids
 barcode.names = read.delim(barcode.path,   
                           header = FALSE,
                           stringsAsFactors = FALSE)
 #set the column and row names for the matrix
 colnames(mat) = barcode.names$V1
 rownames(mat) = feature.names$V2

  # Filter matrix
  # Remove spots with no signal
  mat <- mat[,apply(mat,2,max)>0]
  # Remove genes with no signal/constant signal
  mat <- mat[apply(mat,1,sd)>0,]

  #transform matrix
  mat<-log2(mat+1)
  #convert matrix to a sparse matrix
  sparse.mat <- Matrix(mat, sparse = T)
  #save the matrix object
  writeMM(obj = sparse.mat, file=paste0(samples[i],"/",samples[i],"_Feature-BarcodeMatrix.mtx"))

  #create the variables for the gene names and spot names
  geneNames <- rownames(sparse.mat)
  cellNames <- colnames(sparse.mat)
  #save the gene names and row names
  save(geneNames,cellNames,file=paste0(samples[i],"/",samples[i],"_Names.rdata"))
}

# Create combined matrix to aggregate samples
cmat <- NULL
exgenes <- list()
for (i in 1:length(samples)){
  # Set up the path names for the samples
  barcode.path <- paste0(data_dir,"/",samples[i],'/filtered_feature_bc_matrix/barcodes.tsv.gz')
  features.path <- paste0(data_dir,"/", samples[i],"/filtered_feature_bc_matrix/features.tsv.gz")
  matrix.path <- paste0(data_dir,"/", samples[i],"/filtered_feature_bc_matrix/matrix.mtx.gz")
  #read in the matrix
  mat <- readMM(file = matrix.path)
  #read in the gene names
  feature.names = read.delim(features.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  #read in the spot ids
  barcode.names = read.delim(barcode.path,   
                             header = FALSE,
                             stringsAsFactors = FALSE)
  #set the column and row names for the matrix
  colnames(mat) = paste0(barcode.names$V1,":", samples[i])
  rownames(mat) = feature.names$V2

  #Filter matrix
  # Subset to the immune genes
  keep <- intersect(rownames(mat),geneVector)
  mat <- mat[keep,]
  # Remove spots with no signal
  mat <- mat[,apply(mat,2,max)>0]
  
  # Check which genes have signal/constant signal and save them in a list to remove downstream 
   exgenes<- c(exgenes, rownames(mat[apply(mat,1,sd)>0,]))
  
  #transform matrix
  mat<-log2(mat+1)
  
  #convert matrix to a sparse matrix
  sparse.mat <- Matrix(mat, sparse = T)
  cmat <- cbind(cmat,sparse.mat)
}

#save the matrix object
writeMM(obj = cmat, file="pseudoaggr_CoGAPSprep.mtx")
#create the variables for the gene names and spot names
geneNames <- rownames(cmat)
cellNames <- colnames(cmat)
save(geneNames,cellNames,file="pseudoaggrNames.rdata")

#create the params object for running cogaps
#load the gene names, cell names,and matrix for the combined matrix

load("pseudoaggrNames.rdata")

params <- new("CogapsParams")
params <- CogapsParams(
  sparseOptimization=TRUE,
  nPatterns=25,
  seed=123,
  singleCell=TRUE,
  geneNames=geneNames,
  sampleNames=cellNames,
  nIterations=50000,
  distributed='single-cell'
)
#set the number of sets
  # total genes / nsets >1000 and # spots/ nsets >1000
params <- setDistributedParams(params, nSets=2)
saveRDS(params, file="pseudoaggr_CogapsParams.rds")

#create params files for each sample separately:
samples <- samples[-c(1:3)]
sets = NA
#set the number of sets
# total genes / nsets >1000 and # spots/ nsets >1000
for (i in 1:length(samples)){
  sets=NA
  matrix.path <- paste0(samples[i],"/",samples[i],"_Feature-BarcodeMatrix.mtx")
  names.path <- paste0(samples[i],"/",samples[i],"_Names.rdata")
  load(names.path)
  print(samples[i])
  print(length(cellNames))
  if(length(cellNames)<2000){sets=1}
  if(length(cellNames)>2000 & length(cellNames)<3000){sets=2}
  if(length(cellNames)>3000){sets=3}
  print(sets)
  params <- new("CogapsParams")
  params <- CogapsParams(
    sparseOptimization=TRUE,
    nPatterns=10,
    seed=123,
    singleCell=TRUE,
    geneNames=geneNames,
    sampleNames=cellNames,
    nIterations=50000,
    distributed='single-cell'
  )
params <- setDistributedParams(params, nSets=sets)
saveRDS(params, file=paste0(samples[i],"/",samples[i],"CogapsParams.rds"))
}

