# Load all necessary packages
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(sesame)
library(stringr)
library(reshape)
library(lattice)
library(Matrix) 
library(data.table) 
library(dplyr)


# get the EPIC annotation data
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(annEPIC)

# Convert row names to a new dataframe with a single column named 'ID_source'
df_with_id_source <- data.frame(ID_source = rownames(annEPIC))

# Save the new dataframe as an Excel file
write.xlsx(df_with_id_source, "EPICv1_ID_source_data.xlsx", rowNames = FALSE)


# read in the sample sheet for the experiment
targets <- read.metharray.sheet("./step01_version02_public_data/data/step01_version02_epic_data", pattern="SampleSheet.csv")
targets

# read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets=targets, force = TRUE)
rgSet

# give the samples descriptive names
targets$ID <- paste(targets$Prognosis_simple,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$GSM_ID
rgSet


# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)



pal <- brewer.pal(8,"Dark2")

# Generate the barplot and capture its return value for use with legend placement
par(mfrow=c(1,1))
bp <- barplot(colMeans(detP), col = pal[factor(targets$Prognosis_simple)], las = 3,
              cex.names = 0.8, ylab = "Mean detection p-values")
legend("topleft", legend = levels(factor(targets$Prognosis_simple)), fill = pal, bg = "white")


# controlStripPlot(rgSet, controls="BISULFITE CONVERSION II")

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
plotQC(qc)

# visualise what the data looks like before and after normalisation
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Prognosis_simple,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Prognosis_simple)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Prognosis_simple,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Prognosis_simple)), 
       text.col=brewer.pal(8,"Dark2"))





# MDS plots to look at largest sources of variation
par(mfrow=c(1,1))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Prognosis_simple)])
legend("topleft", legend=levels(factor(targets$Prognosis_simple)), text.col=pal,
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
        col=pal[factor(targets$Prognosis_simple)])
legend("topleft", legend=levels(factor(targets$Prognosis_simple)), text.col=pal,
       bg="white", cex=0.7)



# calculate M-values for statistical analysis
#mVals <- getM(mSetSqFlt)
#head(mVals[,1:5])

#mVals_df <- as.data.frame(mVals)

# Write the data frame to a CSV file
#write.csv(mVals_df, file = "./public_data_mvals.csv", row.names = TRUE)

bVals <- getBeta(mSetSqFlt)

# Convert the matrix to a data frame
bVals_df <- as.data.frame(bVals)

rows_to_keep <- apply(bVals_df, 1, function(row) all(row < 0.3 | row > 0.6))
bVals_df <- bVals_df[rows_to_keep, ]
bVals <- as.matrix(bVals_df)
mVals <- BetaValueToMValue(bVals)

mVals_df <- as.data.frame(mVals)

write.csv(mVals_df, file = "./public_data_mvals.csv", row.names = TRUE)


# Write the data frame to a CSV file
write.csv(bVals_df, file = "./step01_version02_public_data/results/public_data_bvals.csv", row.names = TRUE)



head(bVals[,1:5])

par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Prognosis_simple, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("topright", legend = levels(factor(targets$Prognosis_simple)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Prognosis_simple, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Prognosis_simple)), 
       text.col=brewer.pal(8,"Dark2"))



# Differential CpG methylation analysis

# this is the factor of interest
prognosis <- factor(targets$Prognosis_simple)

# use the above to create a design matrix
design <- model.matrix(~0+prognosis, data=targets)
colnames(design) <- c(levels(prognosis))

# fit the linear model 
fit1 <- lmFit(mVals, design)

# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(indolent-benign, levels=design)


# fit the contrasts
fit2 <- contrasts.fit(fit1, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))


# get the table of results for the first contrast (naive - rTreg)
annEPICSub <- annEPIC[match(rownames(mVals),annEPIC$Name),
                      c(1:4,12:19,24:ncol(annEPIC))]

DMPs <- topTable(fit2, num=Inf, coef=1, genelist=annEPICSub)
head(DMPs)

write.table(DMPs, file="./step01_version02_public_data/results/indolent_vs_benign_dmps.csv", sep=",", row.names=FALSE)


# Differentially methylated regions
myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             epicv2Remap = FALSE, coef = "indolent - benign", 
                             arraytype = "EPICv1")

str(myAnnotation)

DMRs <- dmrcate(myAnnotation, lambda=1000, pcutoff = 0.001)
results.ranges <- extractRanges(DMRs)
results.ranges

# Save the GRanges object
saveRDS(results.ranges, file = "./step01_version02_public_data/results/indolent_vs_benign_granges.rds")

# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$Prognosis_simple))]
names(groups) <- levels(factor(targets$Prognosis_simple))
cols <- groups[as.character(factor(targets$Prognosis_simple))]

# draw the plot for the top DMR
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges, dmr = 1, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "EPICv1", genome = "hg19")

