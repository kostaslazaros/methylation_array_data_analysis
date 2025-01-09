library(dplyr)
library(tibble)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
library(minfi)
library(tidyr)
library(readxl)
library(openxlsx)

# get the EPIC annotation data
annEPIC_v2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
head(annEPIC_v2)
annEPIC_v2$GencodeV41_Name

df <- read_excel("./intersection_lists/high_grade_intersection.xlsx")

# Convert to a data frame
df <- as.data.frame(df)

rownames(df) <- df[[1]]
df

# get the table of results for the first contrast
annEPICSub_v2 <- annEPIC_v2[match(rownames(df),annEPIC_v2$Name),
                            c(1:4,12:19,24:ncol(annEPIC_v2))]

df$Genes <- annEPICSub_v2$GencodeV41_Name

df <- separate(df, Genes, into = c("GeneName", "Rest"), sep=";", extra="drop")

df <- subset(df, select=-Rest)

df <- df[df$GeneName != "", ]

df_filtered <- df
# Remove unecessary 
rm(df)
rm(annEPIC_v2)
rm(annEPICSub_v2)

grouped_df <- df_filtered %>%
  group_by(GeneName) %>%
  summarise(
    CpGs = toString(unique(Name)),
    diff_meth = toString(unique(diff_meth))
  )

final_df <- grouped_df %>%
  select(GeneName, diff_meth) %>%
  as.data.frame()


write.xlsx(final_df, "./intersection_lists/genes_from_dmps/genes_high_grade_dmp_list.xlsx", rowNames = FALSE)
