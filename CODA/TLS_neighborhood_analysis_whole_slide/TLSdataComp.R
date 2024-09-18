library(openxlsx)
library(pheatmap)

perms <- grep( "TLS_data_", list.files(), value = TRUE)
#list directories
slides <- list.files(paste0("./", perms[1]))
kept.slides <- c("S16-57877_2T2_HE",
                 "S16-57877_2T4_HE",
                 "S18-05538_FSC_HE",
                 "S19-27124_2O_HE",
                 "S19-28733_3O_HE",
                 "S21-27914_4O_HE",
                 "S21-27914_4N_HE",
                 "S19-05400_2T4_HE",
                 "S19-05400_2T2_HE",
                 "S20-21634_3N_HE",
                 "S17-38021_2T3_HE")

# combine all TLS from the slides above:

df <- data.frame()
for (slide in kept.slides){
  excel_path <- paste0("./", perms[4], "/", slide, "/",slide,".xlsx")
  mat <- read.xlsx(excel_path)
  df <- rbind(df, mat)
}

range(df$TLS)

df$nontissue <- df$nontissue + df$fiducials
df <- df[,-c(5)]
#df <- df[,colnames(df) != "fiducials" & colnames(df) != "TLS" & colnames(df) != "nontissue" ]
#pheatmap(df, show_rownames = FALSE)
 
# kmeans_result <- kmeans(scale(df), centers = 10)
df$TLSid <- rownames(df)
stacks <- rownames(df)
df <- df[order(df$ecm),]
df$stackorder <- stacks
df_melted <- melt(df, variable.name = "variable", value.name = "Value")
df_melted$stackorder <- as.numeric(df_melted$stackorder)
#df_melted[df_melted$variable == "fiducials",]$variable <- "nontissue"
# # Create stacked bar plot

# df_melted <- df %>%
#   gather(key = "variable", value = "value", -c(TLSid,stackorder))

df_melted$variable <- factor(df_melted$variable, levels = c("ecm", "PDAC","islet","duct","TLS","acini","fat", "nerve","vasculature","nontissue"))

ggplot(df_melted,aes(x=stackorder, y=Value, fill=variable)) + geom_col(stat="position") +
  scale_fill_manual(values = c("ecm" = "pink", 
                               "PDAC" = "orange", 
                               "islet" = "turquoise",
                               "duct" = "blue", 
                               "TLS" = "white",
                               "acini" = "purple", 
                               "fat" = "yellow",
                               "nerve" = "darkgreen", 
                               "vasculature" = "#33FF57", 
                               "nontissue" ="gray")) +  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5, margin = margin(t = 0, b = -10)),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.spacing = unit(0, "lines")
  ) 

# df <- df[df$ecm<50,]
# pheatmap(df, scale = "column")
# 
# df <- df[,colnames(df) != "ecm" & colnames(df) !="fat"]
# pheatmap(scale(df))



# Sample data matrix
set.seed(123) # For reproducibility

# Load the stats package for prcomp, though it should be available by default
library(stats)

# Running PCA
pca_result <- prcomp(df, center = TRUE, scale. = TRUE)

# center = TRUE means subtracting the mean of each variable (mean-centering)
# scale. = TRUE means scaling each variable to have unit variance before applying PCA

# Print summary of PCA results to get variance explained by each principal component
summary(pca_result)

# Access the principal components (rotation matrix)
pca_result$rotation

# Access the scores (coordinates of the data in the principal component space)
pca_result$x

# Plot variance explained by each principal component
plot(pca_result, type = "l") # 'l' for line plot

# Plotting the first two principal components
plot(pca_result$x[, 1], pca_result$x[, 2], xlab = "PC1", ylab = "PC2", main = "PCA Plot")
text(pca_result$x[, 1], pca_result$x[, 2], labels = row.names(pca_result$x), pos = 4)

library(ggplot2)
# Convert PCA scores to a data frame
pca_scores <- as.data.frame(pca_result$x)

# Add row names as a new column for labeling
pca_scores$ID = row.names(pca_scores)

# Create the plot
ggplot(pca_scores, aes(x = PC1, y = PC2, label = ID)) +
  geom_point() + 
  geom_text(vjust = -0.5, hjust = 0.5) +
  xlab("Principal Component 1") + 
  ylab("Principal Component 2") +
  ggtitle("PCA Plot") +
  theme_minimal()+  stat_ellipse()

#Another common way to visualize PCA results is to use a biplot, which shows both the scores (observations) and the loadings (variables) on the same plot, providing a more comprehensive view of how variables contribute to the principal components.
biplot(pca_result,cex =c(0.5,1.5))

