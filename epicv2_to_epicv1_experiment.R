library(minfi)
library(sesame)

# Ensure sesame data is available
sesameDataCacheAll()

# getBetas is the default
#betas = openSesame("./diagenode_data", func = getBetas) 
#betas

betas <- read.csv(file="./diagenode_data_results/data/diagenode_data_beta_vals.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)

# Convert the data frame to a matrix
betas <- as.matrix(betas)

# use empirical evidence in mLiftOver
mapping = sesameDataGet("liftOver.EPICv2ToEPIC")

betas_epicv1 <- mLiftOver(betas, "EPIC", impute=TRUE, mapping = mapping)
betas_epicv1

betas_epicv1_df <- as.data.frame(betas_epicv1)

mvals <- BetaValueToMValue(betas_epicv1)
mvals_df <- as.data.frame(mvals)

write.csv(betas_epicv1_df, file = "./sesame_results/data/diagenode_data_converted_bvals.csv", row.names = TRUE)
write.csv(mvals_df, file = "./sesame_results/data/diagenode_data_converted_mvals.csv", row.names = TRUE)
