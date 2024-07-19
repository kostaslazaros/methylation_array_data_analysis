library(ggplot2)
library(pheatmap)
library(minfi)
library(dplyr)

#df <- read.csv(file="./diagenode_data_results/data/diagenode_data_heatmap_df.csv")

df <- read.csv(file="./diagenode_public_combination_results/data/diagenode_public_combo_data_heatmap_df.csv")


#targets <- read.metharray.sheet("./diagenode_data", pattern="SampleSheet.csv")

targets <- read.metharray.sheet("./diagenode_public_combination", pattern="SampleSheet.csv")

#targets <- subset(targets, GSM_ID != "GSEDIAG01")
#targets <- subset(targets, GSM_ID != "GSEDIAG02")

# Create a new annotation data frame with column names and annotation values
annotation_df <- data.frame(
  Tissue = targets$Prognosis_simple
)
rownames(annotation_df) <- colnames(df)


# Generate the heatmap with annotations, without the column name in the legend
pheatmap(as.matrix(df),
         color = colorRampPalette(c("green", "white", "red"))(50), # Color palette
         clustering_distance_cols = "euclidean", # Clustering distance for columns
         clustering_method = "complete", # Clustering method
         show_rownames = FALSE, # Show row names
         show_colnames = TRUE, # Show column names
         cluster_rows = FALSE, # Do not cluster rows
         annotation_col = annotation_df, # Add column annotations
         annotation_names_col = FALSE) # Do not show annotation column names in the legend




write.table(df_filtered, file="./step01_heatmap_df_filtered.csv", sep=",", row.names=TRUE)
