# Load all necessary packages
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(openxlsx)
library(Gviz)
library(DMRcate)
library(stringr)
library(reshape)
library(sesame)
library(lattice)
library(Matrix) 
library(data.table) 
library(dplyr)


# get the EPIC V2 annotation data
annEPICv2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
head(annEPICv2)

# read in the sample sheet for the experiment
targets <- read.metharray.sheet("./diagenode_data_v01_v02_combo", pattern="SampleSheet.csv")
targets

# read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets=targets, force = TRUE)
rgSet

# give the samples descriptive names
targets$ID <- paste(targets$Prognosis,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$GSM_ID
rgSet

# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)

pal <- brewer.pal(8,"Dark2")

# Generate the barplot and capture its return value for use with legend placement
par(mfrow=c(1,1))
bp <- barplot(colMeans(detP), col = pal[factor(targets$Prognosis)], las = 3,
              cex.names = 0.8, ylab = "Mean detection p-values")
legend("topleft", legend = levels(factor(targets$Prognosis)), fill = pal, bg = "white")

mean_detp_lst <- colMeans(detP)

mean_detp_df <- data.frame(mean_detp = mean_detp_lst)
mean_detp_df <- data.frame(Sample = rownames(mean_detp_df), mean_detp_df)
mean_detp_df
write.xlsx(mean_detp_df, "./mean_detp_diageonode_v01_v02_combo.xlsx")

# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)

qc <- getQC(mSetRaw)
median_intens_df <- as.data.frame(qc)
median_intens_df <- data.frame(Sample = rownames(median_intens_df), median_intens_df)
median_intens_df
write.xlsx(median_intens_df, "./median_intensity_diageonode_v01_v02_combo.xlsx")
