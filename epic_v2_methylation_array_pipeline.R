# Load all necessary packages
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(reshape)
library(sesame)
library(lattice)
library(Matrix) 
library(data.table) 
library(dplyr)
library(openxlsx)


# get the EPIC V2 annotation data
annEPICv2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
head(annEPICv2)

# Convert row names to a new dataframe with a single column named 'ID_source' 
#df_with_id_source <- data.frame(ID_source = rownames(annEPICv2))

# Select the specified columns from annEPICv2 and combine with ID_source
#selected_columns <- annEPICv2[, c("chr", "pos", "strand", "Relation_to_Island", 
#                                  "HMM_Island", "Regulatory_Feature_Group", 
#                                  "Phantom5_Enhancers", "Methyl450_Enhancer",
#                                  "Methyl450_Loci", "EPICv1_Loci")]

# Combine ID_source column with the selected columns
#final_df <- cbind(df_with_id_source, selected_columns)

# Save the new dataframe as an Excel file
#write.xlsx(final_df, "./EPICv2_source_data.xlsx", rowNames = TRUE)


# read in the sample sheet for the experiment
targets <- read.metharray.sheet("./diagenode_v2_data", pattern="SampleSheet.csv")
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

# Control Strip Plot (for Bisulfite Conversion I and II)
.isRGOrStop <- function(object) {
  if (!is(object, "RGChannelSet")) {
    stop("object is of class '", class(object), "', but needs to be of ",
         "class 'RGChannelSet' or 'RGChannelSetExtended'")
  }
}


mycontrolStripPlot <- function(rgSet,
                               controls = c("BISULFITE CONVERSION I",
                                            "BISULFITE CONVERSION II"),
                               sampNames = NULL, xlim = c(5, 17)) {
  .isRGOrStop(rgSet)
  
  r <- getRed(rgSet)
  g <- getGreen(rgSet)
  
  for (controlType in controls) {
    ctrlAddress <- getControlAddress(rgSet, controlType = controlType)
    
    # Red channel
    ctlWide <- as.matrix(log2(r[ctrlAddress, ,drop=FALSE]))
    if (!is.null(sampNames)) colnames(ctlWide) <- sampNames
    ctlR <- melt(ctlWide, varnames = c("address", "sample"))
    
    # Green channel
    ctlWide <- as.matrix(log2(g[ctrlAddress, ,drop=FALSE]))
    if (!is.null(sampNames)) colnames(ctlWide) <- sampNames
    ctlG <- melt(ctlWide, varnames = c("address", "sample"))
    
    # Plot
    ctl <- rbind(
      cbind(channel = "Red", ctlR),
      cbind(channel = "Green", ctlG))
    if (any((ctl$value < xlim[1]) | (ctl$value > xlim[2]))) {
      message("Warning: ", controlType, " probes outside plot range")
    }
    fig <- xyplot(
      x = sample ~ value | channel,
      groups = channel, horizontal = TRUE, pch = 19,
      col = c("darkgreen", "darkred"),
      xlab = "Log2 Intensity",
      xlim = xlim,
      main = paste("Control:", controlType),
      layout = c(2, 1),
      data = ctl,
      panel = function(x, y,...) {
        panel.stripplot(x, y,...)
        panel.abline(h = (as.numeric(y) - 0.5), lty = 3, col = "grey70")
      })
    print(fig)
  }
}

mycontrolStripPlot(rgSet, "BISULFITE CONVERSION II")

# remove poor quality samples
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
rgSet

# remove poor quality samples from targets data
targets <- targets[keep,]

# remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)

# normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessFunnorm(rgSet)


# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)

qc <- getQC(mSetRaw)
plotQC(qc, badSampleCutoff = 10)


# visualise what the data looks like before and after normalisation
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Prognosis,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Prognosis)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Prognosis,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Prognosis)), 
       text.col=brewer.pal(8,"Dark2"))



# MDS plots to look at largest sources of variation
par(mfrow=c(1,1))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Prognosis)])
legend("topleft", legend=levels(factor(targets$Prognosis)), text.col=pal,
       bg="white", cex=0.7)


# FILTERING

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)


mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt


# exclude cross reactive probes 
xReactiveProbes <- maxprobes::xreactive_probes(array_type = "EPIC")
xReactiveProbes <- unlist(xReactiveProbes)

keep2 <- !(featureNames(mSetSqFlt) %in% xReactiveProbes)
table(keep2)

mSetSqFlt <- mSetSqFlt[keep2,] 
mSetSqFlt

# MDS plots to look at largest sources of variation
par(mfrow=c(1,1))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Prognosis)])
legend("topleft", legend=levels(factor(targets$Prognosis)), text.col=pal,
       bg="white", cex=0.7)


bVals <- getBeta(mSetSqFlt)

# Convert the matrix to a data frame
bVals_df <- as.data.frame(bVals)

# Assuming your data frame is named df and the column of interest is "column_name"
#min_value <- min(bVals_df$GSEDIAG19, na.rm = TRUE)
#max_value <- max(bVals_df$GSEDIAG19, na.rm = TRUE)

#min_value
#max_value

# Assuming your data frame is named df and the column of interest is "column_name"
#zero_count <- sum(bVals_df$GSEDIAG19 < 0.1, na.rm = TRUE)
#zero_count

# Remove values that are either <= 0 or >=1 in 
# order to avoid inf error in PCA dim reduction
rows_to_keep <- apply(bVals_df, 1, function(row) all(row > 0 & row < 1))
bVals_df <- bVals_df[rows_to_keep, ]

rows_to_keep2 <- apply(bVals_df, 1, function(row) all(row < 0.3 | row > 0.6))
bVals_df <- bVals_df[rows_to_keep2, ]

min(min(bVals_df))

bVals <- as.matrix(bVals_df)

mVals <- BetaValueToMValue(bVals)

mvals_df <- as.data.frame(mVals)
mvals_df
# Write the data frame to a CSV file
write.csv(mvals_df, file = "./diagenode_v2_mvals.csv", row.names = TRUE)

head(bVals[,1:5])

# Differential CpG methylation analysis

# this is the factor of interest
prognosis <- factor(targets$Prognosis)

# use the above to create a design matrix
design <- model.matrix(~0+prognosis, data=targets)
colnames(design) <- c(levels(prognosis))

# fit the linear model 
fit1 <- lmFit(mVals, design)

# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(AVPC-Indolent, levels=design)

# fit the contrasts
fit2 <- contrasts.fit(fit1, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

# get the table of results for the first contrast (naive - rTreg)
annEPICv2Sub <- annEPICv2[match(rownames(mVals),annEPICv2$Name),
                      c(1:4,12:19,24:ncol(annEPICv2))]

DMPs <- topTable(fit2, num=Inf, coef=1, genelist=annEPICv2Sub)
head(DMPs)

write.table(DMPs, file="./avpc_vs_indolent_DMPs.csv", sep=",", row.names=FALSE)

mVals_filt1 <- rmSNPandCH(mVals)
mVals_filt2 <- rmPosReps(mVals_filt1, filter.strategy="mean")

# Differentially methylated regions
myAnnotation <- cpg.annotate(object = mVals_filt2, datatype = "array", what = "M",
                             arraytype = "EPICv2", epicv2Remap = TRUE, 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "AVPC - Group4", fdr=0.001)

DMRs <- dmrcate(myAnnotation, lambda=1000, pcutoff=0.001)
results.ranges <- extractRanges(DMRs, genome="hg38")



# Save the GRanges object
saveRDS(results.ranges, file = "./diagenode_v2_results/avpc_vs_group4/data/avpc_vs_group4_granges_object.rds")

# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$Prognosis_third_grouping))]

names(groups) <- levels(factor(targets$Prognosis_third_grouping))

cols <- groups[as.character(factor(targets$Prognosis_third_grouping))]

DMR.plot(ranges = results.ranges, dmr = 1, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "EPICv2", genome = "hg38")





