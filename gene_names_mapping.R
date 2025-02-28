library(dplyr)
library(tibble)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
library(minfi)
library(tidyr)

# get the EPIC annotation data
annEPIC_v2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
head(annEPIC_v2)
annEPIC_v2$GencodeV41_Name


#annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
#head(annEPIC)


df <- read.csv(file="./diagenode_v2_mvals.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)

# get the table of results for the first contrast
annEPICSub_v2 <- annEPIC_v2[match(rownames(df),annEPIC_v2$Name),
                     c(1:4,12:19,24:ncol(annEPIC_v2))]

# get the table of results for the first contrast
#annEPICSub <- annEPIC[match(rownames(df),annEPIC$Name),
                            #c(1:4,12:19,24:ncol(annEPIC))]

df$Genes <- annEPICSub_v2$GencodeV41_Name

df <- separate(df, Genes, into = c("GeneName", "Rest"), sep=";", extra="drop")

df <- subset(df, select=-Rest)

df <- df[df$GeneName != "", ]


df_filtered <- df

# Remove unecessary 
rm(df)
rm(annEPIC_v2)
rm(annEPICSub_v2)

#rm(annEPIC)
#rm(annEPICSub)

# Turn cpg rownames into a column

df_filtered <- df_filtered %>% rownames_to_column(var="CpG")

grouped_df <- df_filtered %>%
  group_by(GeneName) %>%
  summarise(across(where(is.numeric), mean))

final_df <- grouped_df %>% column_to_rownames(var = "GeneName")

write.table(final_df, file="./diagenode_v2_mvals_genes.csv", sep=",", row.names=TRUE)
