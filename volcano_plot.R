library(tidyverse)
library(RColorBrewer)
library(ggrepel)

df <- read.csv("./diagenode_public_combination_results/data/combination_indolent_vs_metastasis_DMPs.csv")
df$neg_log10_pval <- -log10(df$P.Value)
df

# Biostatsquid theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))


# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
df$diffexpressed <- "NO"
  df$diffexpressed[df$logFC > 1 & df$P.Value < 0.001] <- "UP"
  df$diffexpressed[df$logFC < -1 & df$P.Value < 0.001] <- "DOWN"
  head(df[order(df$P.Value) & df$diffexpressed == 'DOWN', ])


  
# df_annotated <- df[df$genes %in% genes_to_annotate,]
  # For Downregulated genes: prioritize by -log10(pval) first (high significance), then by logfoldchanges (most negative)
  top_downregulated <- df[df$diffexpressed == "DOWN", ][order(-df$neg_log10_pval[df$diffexpressed == "DOWN"], df$logFC[df$diffexpressed == "DOWN"]), ][1:5,]
  
  top_upregulated <- df[df$diffexpressed == "UP", ][order(-df$neg_log10_pval[df$diffexpressed == "UP"], -df$logFC[df$diffexpressed == "UP"]), ][1:5,]

  #)
  
  
  # top_up <- c("cg23688334", "cg11059561", "cg02968883", "cg04656757", "cg12366974", "cg02695267")
  
  # top_down <- c("cg06788803", "cg27294950", "cg26687120", "cg21511817", "cg25339705", "cg11682124")
  
  
  # Combine into a list
  top_dmps_combined <- c(top_downregulated$Name, top_upregulated$Name)
  top_dmps_combined
  
  # Convert the vector to a data frame
  # dmps_df <- data.frame(GeneNames=top_genes_combined)
  # write.csv(genes_df, "yang_fwd_pericytes_disease_vs_control_degs.csv", row.names=FALSE)
  
  df_annotated <- df[df$Name %in% top_dmps_combined,]
  
  
  
  # Note. with coord_cartesian() even if we have genes with p-values or log2FC ourside our limits, they will still be plotted.
  # Your adjusted plot code
p <- ggplot(data=df, aes(x=logFC, y= -log10(P.Value), col=diffexpressed)) +
    geom_vline(xintercept=c(-2.5, 2.5), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept=-log10(0.001), col = "gray", linetype = 'dashed') + 
    geom_point(size=2) + 
    scale_color_manual(values=c("#00AFBB", "grey", "#bb0c00"),
                       labels=c("Hypomethylated", "Not significant", "Hypermethylated")) +
    coord_cartesian(ylim=c(0, 11), xlim=c(-6, 6)) +
    labs(color='', x=expression("log"[2]*"FC"), y=expression("-log"[10]*"p-value")) +
    scale_x_continuous(breaks=seq(-6, 6, 1)) +
    ggtitle("DMPs (Indolent vs Metastasis)")
  
  # Adding labels with ggrepel for better visibility and avoiding overlaps
 p + geom_label_repel(data=df_annotated, aes(label=Name, x=logFC, y= -log10(P.Value)),
                       box.padding=0.35, point.padding = 0.5, 
                       size=3, segment.color='grey50', show.legend=FALSE)
  
  ggsave("./step01_volcano_indolent_benign.png", dpi = 600)

    