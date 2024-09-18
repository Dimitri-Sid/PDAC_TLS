
library(stringr)
library(dplyr)
library(ggplot2)
library(ggalluvial)

#### FUNCTIONS 

#Removeslist elements that contain all NA values
checkBlanks<- function(df, cloneCall) {
  count <- NULL
  for (i in seq_along(df)) {
    if (length(df[[i]][,cloneCall]) == length(which(is.na(df[[i]][,cloneCall]))) | 
        length(which(!is.na(df[[i]][,cloneCall]))) == 0 | 
        nrow(df[[i]]) == 0) {
      count <- c(i, count)
    } else {
      next()
    }
  }
  if (!is.null(count)) {
    df <- df[-count]
  }
  return(df)
}

off.the.chain <- function(dat, chain, cloneCall) {
  chain1 <- toupper(chain) #to just make it easier
  if (chain1 %in% c("TRA", "TRD", "IGH")) {
    x <- 1
  } else if (chain1 %in% c("TRB", "TRG", "IGL")) {
    x <- 2
  } else {
   warning("It looks like ", chain, " does not match the available options for `chain = `")
  }
  return(dat)
}

compareClonotypes_modified <- function(mat, cloneCall = "CDR3.aa", metadata = "Anns",
          samples = NULL, clonotypes = NULL, numbers = NULL, graph = "alluvial", 
          exportTable = FALSE) {
  print("* select cloneCall to be any column that has clonotype info etc.")
  print("** select one patient at a time, because there is the issue that that some combinations of group and x_axis will be duplicated, even though they may have different values of stratum. The alluvial plot will track each group (alluvium) across each x-axis through some stratum, and that stratum must be unique. (The alluvia are not allowed to reverse direction.) So, each group-x_axis pair must appear at most once.")
  df <- checkBlanks(mat, cloneCall)
  if (!is.null(numbers) & !is.null(clonotypes)) {
    stop("Make sure your inputs are either numbers or clonotype sequences.")
  }
  Con.df <- NULL
  for (i in seq_along(df)) {
    df[[i]][,cloneCall] <- paste0(df[[i]][,cloneCall] ,":",df[[i]][,metadata]) #combine CT info to separate downstream as metadata
    tbl <- as.data.frame(table(df[[i]][, cloneCall]))
    tbl[, 2] <- tbl[, 2]/sum(tbl[, 2])
    colnames(tbl) <- c("Clonotypes", "Proportion")
    tbl$Sample <- names(df[i])
    Con.df <- rbind.data.frame(Con.df, tbl) # adding cell type info here
  }
  Con.df[,metadata] <-  str_split(Con.df$Clonotypes, ":", simplify = TRUE)[,2]
  Con.df$Clonotypes <-  str_split(Con.df$Clonotypes, ":", simplify = TRUE)[,1]
  if (!is.null(samples)) {
    Con.df <- Con.df[Con.df$Sample %in% samples, ]
  }
  if (!is.null(clonotypes)) {
    Con.df <- Con.df[Con.df$Clonotypes %in% clonotypes, ]
  }
  if (!is.null(numbers)) {
    top <- Con.df %>% top_n(n = numbers, wt = Proportion)
    Con.df <- Con.df[Con.df$Clonotypes %in% top$Clonotypes, 
    ]
  }
  if (nrow(Con.df) < length(unique(Con.df$Sample))) {
    stop("Reasses the filtering strategies here, there is not \n            enough clonotypes to examine.")
  }
  
  Con.df[,"coldatvar"] <- Con.df[,metadata]
  
  if (exportTable == TRUE) {
    return(Con.df)
  }
  
  plot <- ggplot(Con.df, aes(x = coldatvar, fill = Clonotypes, 
                             group = Clonotypes, stratum = Clonotypes, alluvium = Clonotypes, 
                             y = Proportion, label = Clonotypes)) + xlab(metadata)+ theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust =1)) 
  
  if (graph == "alluvial") {
    plot = plot + geom_stratum() + geom_flow(stat = "alluvium")
  }
  
  if (graph == "area") {
    plot = plot + geom_area(aes(group = Clonotypes), color = "black")
  }
  return(plot)
}
