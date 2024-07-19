# Diagenode bvals & mvals reading
bvals_diagenode <- read.csv(file="./sesame_results/data/diagenode_data_converted_bvals.csv", 
                            header=TRUE, row.names=1, stringsAsFactors=TRUE)


mvals_diagenode <- read.csv(file="./sesame_results/data/diagenode_data_converted_mvals.csv", 
                            header=TRUE, row.names=1, stringsAsFactors=TRUE)



# Public data bvals & mvals reading
bvals_public <- read.csv(file="./public_data_results/data/public_data_beta_vals.csv", 
                         header=TRUE, row.names=1, stringsAsFactors=TRUE)


mvals_public <- read.csv(file="./public_data_results/data/public_data_mvals.csv", 
                         header=TRUE, row.names=1, stringsAsFactors=TRUE)



# Find common row names (bvals)
bvals_common_rows <- intersect(rownames(bvals_diagenode), rownames(bvals_public))

# Subset data frames to keep only common rows (bvals)
bvals_diagenode <- bvals_diagenode[bvals_common_rows, ]
bvals_public <- bvals_public[bvals_common_rows, ]

# Concatenate the data frames by columns (bvals)
merged_bvals <- cbind(bvals_diagenode, bvals_public)

# Write the data frame to a CSV file (bvals)
write.csv(merged_bvals, file = "./diagenode_public_combination_results/data/merged_bvals.csv", 
          row.names = TRUE)


# Find common row names (mvals)
mvals_common_rows <- intersect(rownames(mvals_diagenode), rownames(mvals_public))

# Subset data frames to keep only common rows (mvals)
mvals_diagenode <- mvals_diagenode[mvals_common_rows, ]
mvals_public <- mvals_public[mvals_common_rows, ]

# Concatenate the data frames by columns (mvals)
merged_mvals <- cbind(mvals_diagenode, mvals_public)

# Write the data frame to a CSV file (mvals)
write.csv(merged_mvals, file = "./diagenode_public_combination_results/data/merged_mvals.csv", 
          row.names = TRUE)


