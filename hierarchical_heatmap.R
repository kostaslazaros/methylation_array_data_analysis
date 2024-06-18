library(ggplot2)
library(pheatmap)

df <- read.csv(file="./step01_epic_data_heatmap_df.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)

rows_to_keep <- apply(df, 1, function(row) all(row < 0.3 | row > 0.6))


df_filtered <- df[rows_to_keep, ]

# Create heatmap with hierarchical clustering only for columns
pheatmap(as.matrix(df),
         color = colorRampPalette(c("green", "white", "red"))(50), # Color palette
         clustering_distance_cols = "euclidean", # Clustering distance for columns
         clustering_method = "complete", # Clustering method
         show_rownames = FALSE, # Show row names
         show_colnames = TRUE, # Show column names
         cluster_rows = FALSE) # Do not cluster rows

write.table(df_filtered, file="./step01_heatmap_df_filtered.csv", sep=",", row.names=TRUE)
