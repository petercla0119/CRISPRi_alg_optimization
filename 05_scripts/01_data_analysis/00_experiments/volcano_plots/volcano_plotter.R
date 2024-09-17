### VOLCANO PLOTTER ###
### Takes 7 inputs: (1) hits file, (2) allgenes file both from MAGECK. 
### (3) a .csv file with gene name in one column and gene color in the other. 
### (4) the FDR you used in your screen (5) a prefix for the output 
### (6) if you want an overlay or not, "yes" or "no" 
### (7) if you want to annotate genes "yes" or "no"

library(ggplot2)
library(dplyr)

# Get arguments and make dataframes
args <- commandArgs(trailingOnly = TRUE)

hitsfile <- args[1]
allgenes_file <- args[2]
gene_list_path <- args[3]
fdr_input <- as.numeric(args[4])
output_name <- args[5]
overlay <- args[6]
annotate <- args[7]

# df_hits <- read.csv(hitsfile) # sgRNA summary file
# df_allgenes <- read.csv(allgenes_file) # Gene summary file
df_hits <- read.delim("/Users/Claire/Downloads/corrected_tbl.sgrna_summary.txt")
df_allgenes <- read.delim("/Users/Claire/Downloads/corrected_tbl.gene_summary.txt")
df_ntc <- df_allgenes %>% dplyr::filter(grepl('negative_control', Gene))

print("plotting....")

# Setup axes
df_pos <- df_hits %>% filter(epsilon > 0)
df_neg <- df_hits %>% filter(epsilon < 0)

merged_df <- df_allgenes %>% 
  left_join(df_hits, by = "gene", suffix = c("", ".y")) %>% 
  filter(is.na(epsilon.y))

dfnothits <- merged_df %>% select(-ends_with(".y"))

# Define threshold function based on FDR
product_threshold_fdr <- function(df, fdr = fdr_input) {
  maxi <- max(abs(df$product), na.rm = TRUE)
  for (pro in seq(0, maxi, by = 1)) {
    df_thres <- df %>% filter(abs(product) > pro)
    if ((sum(grepl('NTC', df_thres$gene)) / nrow(df_thres)) < fdr) {
      break
    }
  }
  return(list(thres = pro, df_hits = df_thres))
}

result <- product_threshold_fdr(df_allgenes, fdr_input)
thres <- result$thres
df_hits <- result$df_hits

# Create plot
ggplot() +
  geom_point(data = df_allgenes, aes(x = epsilon, y = -log10(pvalue)), color = "grey", alpha = 0.3, size = 1) +
  geom_point(data = df_ntc, aes(x = epsilon, y = -log10(pvalue)), color = "black", alpha = 0.3, size = 1) +
  geom_point(data = df_pos, aes(x = epsilon, y = -log10(pvalue)), color = "thistle", alpha = 0.5, size = 1) +
  geom_point(data = df_neg, aes(x = epsilon, y = -log10(pvalue)), color = "lightblue", alpha = 0.5, size = 1) +
  geom_line(aes(x = seq(0.01, 50, 0.01), y = thres / seq(0.01, 50, 0.01)), linetype = "dotted", color = "black") +
  geom_line(aes(x = -seq(0.01, 50, 0.01), y = thres / seq(0.01, 50, 0.01)), linetype = "dotted", color = "black") +
  xlim(-max(abs(min(df_allgenes$epsilon) - 1), max(df_allgenes$epsilon) + 2), 
       max(abs(min(df_allgenes$epsilon) - 1), max(df_allgenes$epsilon) + 2)) +
  ylim(0, -log10(min(df_allgenes$pvalue)) + 0.5) +
  labs(x = "Normalized Phenotype", y = "-log10 P", title = "Plot") +
  theme_classic() +
  theme(text = element_text(family = "Arial", size = 10))

# Annotate genes
if (gene_list_path != "none") {
  gene_list <- read.csv(gene_list_path, header = FALSE, stringsAsFactors = FALSE)
  gene_dict <- setNames(gene_list$V2, gene_list$V1)
  
  for (i in 1:nrow(df_allgenes)) {
    gene <- df_allgenes$gene[i]
    if (gene %in% names(gene_dict)) {
      phenotype <- df_allgenes$epsilon[i]
      p_value <- -log10(df_allgenes$pvalue[i])
      if (annotate == "yes") {
        ggplot2::annotate("text", x = phenotype, y = p_value, label = gene, size = 3)
      }
      ggplot2::geom_point(aes(x = phenotype, y = p_value), color = gene_dict[gene], size = 2)
    }
  }
}

# Save plot
ggsave(paste0(output_name, "_volcano_plot.png"), width = 6.6, height = 5.28, dpi = 600, bg = "transparent")
