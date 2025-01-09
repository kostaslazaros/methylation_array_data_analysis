# Load all necessary packages
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(reshape)
library(lattice)
library(openxlsx)


# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)

# Convert row names to a new dataframe with a single column named 'ID_source'
df_with_id_source <- data.frame(ID_source = rownames(ann450k))

# Save the new dataframe as an Excel file
write.xlsx(df_with_id_source, "450k_ID_source_data.xlsx", rowNames = FALSE)


# read in the sample sheet for the experiment
targets <- read.metharray.sheet("./GSE127985_data", pattern="SampleSheet.csv")
targets


# read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets=targets)
rgSet


# give the samples descriptive names
targets$ID <- paste(targets$Prognosis,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID
rgSet


# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)



pal <- brewer.pal(8,"Dark2")

# Generate the barplot and capture its return value for use with legend placement
par(mfrow=c(1,1))
bp <- barplot(colMeans(detP), col = pal[factor(targets$Prognosis)], las = 2,
              cex.names = 0.8, ylab = "Mean detection p-values")
legend("topleft", legend = levels(factor(targets$Prognosis)), fill = pal, bg = "white")


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
xReactiveProbes <- read.csv(file=paste("./",
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep2 <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep2)

mSetSqFlt <- mSetSqFlt[keep2,] 
mSetSqFlt


# MDS plots to look at largest sources of variation
par(mfrow=c(1,1))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Prognosis)])
legend("topleft", legend=levels(factor(targets$Prognosis)), text.col=pal,
       bg="white", cex=0.7)



# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
head(mVals[,1:5])


bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])


par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Prognosis, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("topright", legend = levels(factor(targets$Prognosis)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Prognosis, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Prognosis)), 
       text.col=brewer.pal(8,"Dark2"))



# Differential CpG methylation analysis

# this is the factor of interest
prognosis <- factor(targets$Prognosis)

# use the above to create a design matrix
design <- model.matrix(~0+prognosis, data=targets)
colnames(design) <- c(levels(prognosis))

# fit the linear model 
fit1 <- lmFit(mVals, design)

# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(bad-good, levels=design)


# fit the contrasts
fit2 <- contrasts.fit(fit1, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))


# get the table of results for the first contrast (naive - rTreg)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
head(DMPs)
write.table(DMPs, file="DMPs_v2.csv", sep=",", row.names=FALSE)


# Differentially methylated regions
myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "bad - good", arraytype = "450K")


str(myAnnotation)


DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges

# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$Prognosis))]
names(groups) <- levels(factor(targets$Prognosis))
cols <- groups[as.character(factor(targets$Prognosis))]

# draw the plot for the top DMR
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges, dmr = 55, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "450K", genome = "hg19")


