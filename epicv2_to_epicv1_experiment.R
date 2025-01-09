library(minfi)
library(sesame)
library(readxl)
library(openxlsx)

# Ensure sesame data is available
sesameDataCacheAll()

# getBetas is the default
#betas = openSesame("./diagenode_data", func = getBetas) 
#betas

betas <- read.csv(file="./diagenode_v2_bvals.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)

# Convert the data frame to a matrix
betas <- as.matrix(betas)

# use empirical evidence in mLiftOver
mapping = sesameDataGet("liftOver.EPICv2ToEPIC")
write.xlsx(mapping, "./python_part/data/epicv2_to_epic_mapping.xlsx", rowNames = FALSE)



betas_epicv1 <- mLiftOver(betas, "EPIC", impute=TRUE, mapping = mapping)
betas_epicv1

betas_epicv1_df <- as.data.frame(betas_epicv1)

mvals <- BetaValueToMValue(betas_epicv1)
mvals_df <- as.data.frame(mvals)

write.csv(betas_epicv1_df, file = "./sesame_results/data/diagenode_data_converted_bvals.csv", row.names = TRUE)
write.csv(mvals_df, file = "./sesame_results/data/diagenode_data_converted_mvals.csv", row.names = TRUE)
