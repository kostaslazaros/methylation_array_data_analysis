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
library(stringr)
library(reshape)
library(lattice)
library(Matrix) 
library(data.table) 
library(dplyr)


# get the EPIC annotation data
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(annEPIC)



# read in the sample sheet for the experiment
targets <- read.metharray.sheet("./diagenode_public_combination", pattern="SampleSheet.csv")
targets


mvals <- read.csv(file="./diagenode_public_combination_results/data/merged_mvals.csv", 
                  header=TRUE, row.names=1, stringsAsFactors=TRUE)

bvals <- read.csv(file="./diagenode_public_combination_results/data/merged_bvals.csv", 
                  header=TRUE, row.names=1, stringsAsFactors=TRUE)


mvals <- as.matrix(mvals)

bvals <- as.matrix(bvals)




# Differential CpG methylation analysis

# this is the factor of interest
prognosis <- factor(targets$Prognosis_simple)

# use the above to create a design matrix
design <- model.matrix(~0+prognosis, data=targets)
colnames(design) <- c(levels(prognosis))

# fit the linear model 
fit1 <- lmFit(mvals, design)

# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(indolent-metastasis, levels=design)


# fit the contrasts
fit2 <- contrasts.fit(fit1, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))


# get the table of results for the first contrast (naive - rTreg)
annEPICSub <- annEPIC[match(rownames(mvals),annEPIC$Name),
                      c(1:4,12:19,24:ncol(annEPIC))]

DMPs <- topTable(fit2, num=Inf, coef=1, genelist=annEPICSub)
head(DMPs)

write.table(DMPs, file="./diagenode_public_combination_results/data/combination_indolent_vs_metastasis_DMPs.csv", 
            sep=",", row.names=FALSE)


# Differentially methylated regions
myAnnotation <- cpg.annotate(object = mvals, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "indolent - metastasis", arraytype = "EPICv1")


str(myAnnotation)



DMRs <- dmrcate(myAnnotation, lambda=1000, pcutoff = 0.001)
results.ranges <- extractRanges(DMRs)
results.ranges

# Save the GRanges object
saveRDS(results.ranges, file = "./diagenode_public_combination_results/data/diagenode_avpc_vs_high_grade_granges_object.rds")

pal <- brewer.pal(8,"Dark2")

# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$Prognosis_simple))]

names(groups) <- levels(factor(targets$Prognosis_simple))

cols <- groups[as.character(factor(targets$Prognosis_simple))]

DMR.plot(ranges = results.ranges, dmr = 3, CpGs = bvals, phen.col = cols, 
         what = "Beta", arraytype = "EPICv1", genome = "hg19")
