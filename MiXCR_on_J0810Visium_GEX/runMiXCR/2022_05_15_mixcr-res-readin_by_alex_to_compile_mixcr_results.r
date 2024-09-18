#Read-in and prepare mixcr data from visium data. 
#This mixcr output was generated from trust_toassemble.fa files, NOT  
#using mixcr from scratch on entire  fastq data (yet)
# 15 May 2022

#Set up workspace. ----
library(strex)
library(stringr)
library(dplyr)
library(ggplot2)
setwd("/media/aag7319/WDRed/ZZZ_PDAC-Visium/")
#setwd("E:/ZZZ_PDAC-Visium/")
#Functions----
readReport <- function(filedir, cropnames = FALSE){
  files <- dir(filedir)
  r.idx <- grep("clonotypes.ALL.txt", files)
  pts <- vector(mode = "list", length = length(r.idx))
  pnames <- files[r.idx]
  pnames <- unlist(strsplit(pnames, ".clonotypes.ALL.txt"))
  if(cropnames){pnames <- lapply(pnames, crop.name)}
  names(pts) <- pnames
  for(i in 1:length(r.idx)){
    pts[[i]] <- read.delim(paste0(filedir,"/", files[r.idx[i]]))
  }
  return(pts)
}
cleanScores <- function(input){
  #NASTY. Need to learn string manipulation. 
  #Removing numeric scores from gene call fields. 
  for(l in 1:length(input)){
    str <- input[l]
    ins <- sub('.*\\((.*)\\).*', '\\1', str)
    torm <- paste0("(",ins,")")
    str <- substr(str,1,nchar(str) - nchar(torm))
    input[l] <- str
  }
  return(input)
  
}
assignLocus <- function(airr){
  if(dim(airr)[1] == 0){return(NA)}
  
  ch <- c("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")
  locus <- vector(mode = "character", length = dim(airr)[1])
  for(c in 1:length(ch)){
    ch.idx <- unique(grep(ch[c], airr$v_call), 
                    # grep(ch[c], airr$d_call), 
                     grep(ch[c], airr$j_call)) 
                     #grep(ch[c], airr$c_call))
    locus[ch.idx] <- ch[c]
  }

  return(locus)
}
convertAirr <- function(input){
if(dim(input)[1] == 0){return(NA)}
  airr <- data.frame(
    sequence_id = NA,
    clone_id = input$cloneId,
    cloneCount = input$cloneCount,
    cloneFraction = input$cloneFraction,
    productive = input$aaSeqImputedCDR3!= "",
    v_call = cleanScores(input$allVHitsWithScore),
    d_call = cleanScores(input$allDHitsWithScore),
    j_call = cleanScores(input$allJHitsWithScore),
    c_call = cleanScores(input$allCHitsWithScore),
    cdr3 = input$nSeqImputedCDR3,
    cdr3_aa = input$aaSeqImputedCDR3,
    cdr3_length = nchar(input$nSeqImputedCDR3),
    cdr2 = input$nSeqImputedCDR2,
    cdr2_aa = input$aaSeqImputedCDR2,
    cdr1 = input$nSeqImputedCDR1,
    cdr1_aa = input$aaSeqImputedCDR1, 
    sequence = input$targetSequences
   
  )
  
  airr$locus = assignLocus(airr)
  return(airr)
}
returnEmpty <- function(byspot){
  #if empty return index
  if(dim(byspot)[1] == 0){
    return(1)
  }else{
    return(0)
  }
}
cleanCB <- function(barcodes){
  bc <- substr(barcodes,1, 16 )
  return(bc)
}
matchLoc <- function(airr, meta = spot.meta){
  tmp <- paste0(airr$barcode,"VI_",airr$samp)
  i <- match(tmp, paste0(spot.meta$aggr_barcodes, spot.meta$panel))
  airr$imagerow <- spot.meta$imagerow[i]
  airr$imagecol <- spot.meta$imagecol[i]
  airr$aggr_barcodes <- spot.meta$aggr_barcodes[i]
  return(airr)
}

#Load Metadata----
spot.meta <- read.delim("./data/J0810_data_orig/visium_metadata.csv", sep = ",")
spot.meta$aggr_barcodes <- substr(spot.meta$aggr_barcodes,1,16)
#Concat Dat----
snames <- c("113_1", "113_2", "113_3", "113_4", 
            "114_1", "114_2", "114_3", "114_4", 
            "115_1", "115_2", "115_3", "115_4", 
            "116_1", "116_2", "116_3", "116_4")
results <- vector(mode = "list", length = length(snames))
names(results) <- snames 
for(s in 1:length(snames)){
  samp <- snames[s]
  datadir <- paste0("/data/mixcr-res/split-bams-v1/", samp)
  filedir <- paste0(getwd(), datadir) 
  byspot <- readReport(filedir)
  #Double check all barcode names. 
 
    fbars <- substr(names(byspot),1,16)
    actualbars <- spot.meta$aggr_barcodes[which(spot.meta$panel == paste0("VI_",samp))]
    invals <- which(is.na(match(fbars,actualbars)))
    if(length(invals) > 0){
      print(paste0("Sample ", samp))
      print("Mismatch between directory files and the expected spot barcodes!")
      print(paste0("Removing ", length(invals), " barcodes from read-in."))
      byspot <- byspot[-invals]
    }
  
  isEmpty <- as.logical(unlist(lapply(byspot, returnEmpty)))
  notEmpty <- byspot[-which(isEmpty)]
  byspot_airr <- lapply(notEmpty, convertAirr)
  byspot_airr <- byspot_airr[!is.na(byspot_airr)]
  airr <- bind_rows(byspot_airr, .id = "barcode")
  airr$sample <- samp
  airr$pt <- substr(airr$sample,1,nchar(airr$sample) - 2)
  airr$barcode <- cleanCB(airr$barcode)
  airr$imagerow <- NA; airr$imagecol <- NA
  airr <- matchLoc(airr)
  
  airr$sequence_id <- paste0("S",airr$pt, "_", airr$clone_id)
  airr$clone_id <- paste0("MC",airr$sample,".",airr$clone_id)
  results[[s]] <- airr
}

 res <- bind_rows(results)
 results <- split(res, res$sample)
 # write.table(res, file = "./22_06_31_pdac-visium-v1.tsv", sep = "\t", 
 #             row.names = FALSE)
 # save(res, file =  "./22_06_31_pdac-visium-v1.rda")
#Plot locations :') ----
for(s in 1:length(results)){
  samp <- names(results)[s]
  rowlim <- c( min(spot.meta$imagerow[which(spot.meta$panel == paste0("VI_",samp))]), 
               max(spot.meta$imagerow[which(spot.meta$panel == paste0("VI_",samp))]) )
  
  collim <- c( min(spot.meta$imagecol[which(spot.meta$panel == paste0("VI_",samp))]), 
               max(spot.meta$imagecol[which(spot.meta$panel == paste0("VI_",samp))]) )
  
 # aspect_ratio <- (rowlim[2] - rowlim[1]) / (collim[2] - collim[1])
  aspect_ratio <- rowlim[2] / collim[2]
  png(filename = paste0("./results/bcr-maps/", samp,".png"),height = 3500 * aspect_ratio, 
      width = 3500 )
  p1 <- ggplot(results[[s]], aes(y = imagerow, x = imagecol,
                                 shape = locus, color = locus)) + 
    geom_point(size = 4) + 
    theme_void() + 
    theme(legend.position = "none") + 
    scale_color_manual(values = c("IGH" = "green", "IGK" = "red", "IGL" = "blue")) + 
    xlim(0,collim[2]) + 
    ylim(0,rowlim[2])
  
  print(p1)
  dev.off()
}

#Assess clonality----
library(alakazam)
library(shazam)
library(umap)
library(stats)
library(stringdist)
library(ggprism)
dl.igl <- adist(bych$IGL$cdr3, method = "dl",
                           costs = c(insertions = 0.25, deletions = 0.25, 
                                     substitutions = 1))
umap.igl <- umap(d = as.matrix(dl.igl), input = "dist")$layout
kclus <- cutree(hclust(dl.igl), h = 2 )

#IgL----
dfigl <- data.frame(
        seqid = bych$IGL$sequence_id,
        UMAP1 = umap.igl[,1],
        UMAP2 = umap.igl[,2],
        cluster = kclus,
        pt = bych$IGL$pt,
        sample = bych$IGL$sample
        
)



gl <- ggplot(dfigl, aes(x = UMAP1, y = UMAP2)) + 
  geom_point(aes(color = sample), size = 0.75) + 
  theme_prism() + labs(title = "IGL")




#IgK----
dl.igk <- adist(bych$IGK$cdr3, #method = "dl",
                costs = c(insertions = 2, deletions = 2, 
                          substitutions = 1))
umap.igk <- umap(d = as.matrix(dl.igk), input = "dist")$layout
#kclus <- cutree(hclust(dl.igk), h = 2 )

dfigk <- data.frame(
  seqid = bych$IGK$sequence_id,
  UMAP1 = umap.igk[,1],
  UMAP2 = umap.igk[,2],
  #cluster = kclus,
  pt = bych$IGK$pt,
  sample = bych$IGK$sample,
  cloneCount = log(bych$IGK$cloneCount, base = 10)
)



gk <- ggplot(dfigk, aes(x = UMAP1, y = UMAP2)) + 
  geom_point(aes(color = pt), size = 0.75) + 
  theme_prism() + labs(title = "IGK") + 
  scale_fill_viridis()


#IgH-----
dh.igh <- adist(bych$IGH$cdr3, #method = "dl",
                costs = c(insertions = 2, deletions = 2, 
                          substitutions = 1))
umap.igh <- umap(d = as.matrix(dh.igh), input = "dist")$layout
#kclus <- cutree(hclust(dl.igh), h = 2 )

dfigh <- data.frame(
  seqid = bych$IGH$sequence_id,
  UMAP1 = umap.igh[,1],
  UMAP2 = umap.igh[,2],
  #cluster = kclus,
  pt = bych$IGH$pt,
  sample = bych$IGH$sample
)



gh <- ggplot(dfigh, aes(x = UMAP1, y = UMAP2)) + 
  geom_point(aes(color = pt), size = 0.75) + 
  theme_prism() + labs(title = "IGH")

#Check out given mixcr freqs ----
source("./functions/immcantation-toolbox.r")
cll <- countClones(bych$IGL, group = "pt")


sampCloneNorm <- function(bych){
  samp <- names(table(bych$sample))
  nsamp <- length(samp)
  for(s in 1:nsamp){
    sidx <- which(bych$sample == samp[s])
    tot <- sum(bych$cloneFraction[sidx])
    bych$cloneFraction[sidx] <- bych$cloneFraction[sidx] / tot
  }
  
  df <- data.frame(
     clone_id = bych$clone_id,
     cloneFraction = bych$cloneFraction,
     sample = bych$sample,
     pt = bych$pt
  )
  
  df <- split(df, df$sample)
  out <- lapply(df, FUN = function(df){
    o = order(df$cloneFraction, decreasing = TRUE)
    out = df$cloneFraction[o]
    names(out) = df$clone_id[o]
    return(out)
  })
  return(out)
}
nClonesCalc <- function(bych){
  df <- split(bych, bych$sample)
  out <- vector(mode = "integer", length = length(df))
  for(i in 1:length(df)){
    out[i] <- dim(df[[i]])[1]
  }
  names(out) <- names(df)
  return(out)
}


nClonesCalc(bych$IGL)

