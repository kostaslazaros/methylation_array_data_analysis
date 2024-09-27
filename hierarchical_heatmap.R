library(ggplot2)
library(pheatmap)
library(minfi)
library(dplyr)

# Read the data
df <- read.csv(file="./diagenode_v2_results/normal_vs_indolent/data/normal_vs_indolent_heatmap_df.csv")

# Print the structure of df to ensure it is read correctly
str(df)

# Read the targets (sample sheet)
targets <- read.metharray.sheet("./diagenode_v2_normal_vs_indolent_data", pattern="SampleSheet.csv")

# Check if the targets dataframe has the expected structure and column
if (!"Prognosis" %in% colnames(targets)) {
  stop("Column 'Prognosis' not found in the targets dataframe")
}

# Create a new annotation data frame with column names and annotation values
annotation_df <- data.frame(
  Tissue = targets$Prognosis
)
rownames(annotation_df) <- colnames(df)

# Print the annotation_df to ensure it is created correctly
print(annotation_df)

# Perform hierarchical clustering
col_dist <- dist(t(as.matrix(df)), method = "euclidean")
col_hclust <- hclust(col_dist, method = "complete")

# Cut the dendrogram to obtain a specific number of clusters
num_clusters <- 3# Specify the desired number of clusters
col_clusters <- cutree(col_hclust, k = num_clusters)

# Add cluster assignments to the annotation data frame
annotation_df$Clusters <- as.factor(col_clusters)

# Print the updated annotation_df to ensure clusters are added correctly
print(annotation_df)

# Generate the heatmap with annotations, without the column name in the legend
pheatmap(as.matrix(df),
         color = colorRampPalette(c("green", "white", "red"))(50), # Color palette
         clustering_distance_cols = "euclidean", # Clustering distance for columns
         clustering_method = "complete", # Clustering method
         show_rownames = FALSE, # Show row names
         show_colnames = TRUE, # Show column names
         cluster_rows = FALSE, # Do not cluster rows
         annotation_col = annotation_df, # Add column annotations
         annotation_names_col = FALSE, # Do not show annotation column names in the legend
         cutree_cols = num_clusters) # Cut columns dendrogram into specified number of clusters

# Save the filtered data to a CSV file
write.table(df, file="./step01_heatmap_df_filtered.csv", sep=",", row.names=TRUE)