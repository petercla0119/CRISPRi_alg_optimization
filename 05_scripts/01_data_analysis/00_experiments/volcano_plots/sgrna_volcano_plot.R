library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggrepel)
source('clean_preprocess_module.R')
# Functions -----------------------

# The rank_test function is similar to the CRISPR-iNC feature (originally in python) from the Kampmann lab.
rank_test <- function(df) {
  # Filter rows where 'treat_mean' is greater than 20
  df <- df %>% filter(treat_mean > 300)
  
  # Split the dataframe into negative controls and targeting genes
  df_ntc <- df %>% filter(grepl('negative_control', Gene))
  df_targeting <- df %>% filter(!grepl('negative_control', Gene))
  
  # Extract relevant columns from the negative control dataframe
  ntc_sgRNA_p <- df_ntc$p.twosided
  ntc_sgRNA_p_lfc <- data.frame(p.twosided = df_ntc$p.twosided, LFC = df_ntc$LFC)
  
  # Get unique genes in the targeting dataframe
  genes <- unique(df_targeting$Gene)
  num_of_genes <- length(genes)
  
  # Initialize an empty list to store LFC and p-values for each gene
  gene_lfc_p <- list()
  
  # Loop through each gene and perform calculations
  for (gene in genes) {
    df_gene <- df_targeting %>% filter(Gene == gene) %>% arrange(p.twosided)
    
    # Calculate the mean LFC of up to the top 3 rows
    lfc <- mean(df_gene$LFC[1:min(3, nrow(df_gene))])
    
    # Perform the Wilcoxon rank-sum test
    pvalue <- wilcox.test(df_gene$p.twosided[1:min(3, nrow(df_gene))], ntc_sgRNA_p, alternative = "two.sided")$p.value
    
    # Store the results in the list
    gene_lfc_p[[gene]] <- list(lfc = lfc, pvalue = pvalue)
  }
  
  # Set the seed for reproducibility
  set.seed(10)
  # Multiple Testing Correction - sample 5 NTCs and calculate KD phenotype scores and pvalues
  # Perform randomization for the negative controls
  for (j in seq_len(num_of_genes)) {
    # Shuffle the ntc_sgRNA_p_lfc list
    ntc_sgRNA_p_lfc <- ntc_sgRNA_p_lfc[sample(nrow(ntc_sgRNA_p_lfc)), ]
    
    # Select the first 5 elements after shuffling
    ntc_selected <- ntc_sgRNA_p_lfc[1:5, ]
    
    # Extract the p-values and calculate the LFC
    ntc_selected_p <- ntc_selected$p.twosided
    ntc_lfc <- mean(ntc_selected$LFC[order(ntc_selected_p)][1:min(3, length(ntc_selected_p))])
    
    # Perform the Wilcoxon rank-sum test
    ntc_pvalue <- wilcox.test(ntc_selected_p, ntc_sgRNA_p, alternative = "two.sided")$p.value
    # Store the results in the list
    gene_lfc_p[[paste0('negative_control_', j)]] <- list(lfc = ntc_lfc, pvalue = ntc_pvalue)
  }
  return(gene_lfc_p)
}

product_threshold_fdr <- function(df, fdr) {
  # Find the maximum absolute value in the 'product' column
  maxi <- max(abs(df$product), na.rm = FALSE)
  
  # Iterate through the sequence from 0 to maxi with a step of 0.1
  for (pro in seq(0, maxi, by = 0.1)) {
    # Filter the dataframe based on the current threshold
    df_thres <- df %>% filter(abs(product) > pro)
    
    # Check if the proportion of 'NTC' entries is below the FDR threshold
    if ((sum(grepl('negative_control', df_thres$gene)) / nrow(df_thres)) < fdr) {
      break
    }
  }
  # Return the threshold and the filtered dataframe
  return(list(thres = pro, df_thres = df_thres))
}

# Import data ----------------------------
# c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#999999")
raw_counts <- read.delim("/Users/Claire/Library/CloudStorage/Box-Box/MSTP Research Rotations/Summer_2024_Might_Ward_Lab/CRISPRi_alg_optimization/00_references_and_code/files_from_others/00_JH/01_data/00_processed_data/00_mageck_counts/counts/all.count.txt", sep = "\t", header = TRUE)
# stmn2_sgrna_summary <- read.delim("/Users/Claire/Downloads/corrected_tbl.sgrna_summary.txt")
stmn2_sgrna_summary <- separate_sgrna_summary("/Users/Claire/Downloads/thresh_300_.sgrna_summary.txt")
stmn2_sgrna_summary <- separate_sgrna_summary("/Users/Claire/Downloads/corrected_tbl.sgrna_summary.txt")

stmn2_sgrna_summary <- calculate_avg_sd_cv_sgrna_summary(stmn2_sgrna_summary)
stmn2_gene_summary <- read.delim("/Users/Claire/Downloads/corrected_tbl.gene_summary.txt")
stmn2_gene_summary <- stmn2_gene_summary %>% dplyr::rename(Gene = id)
# df_ntc <- stmn2_sgrna_summary %>% dplyr::filter(Gene =="negative_control") %>% 
#   dplyr::select(sgrna, Gene, LFC,
#                 pvalue = p.twosided,
#                 fdr = FDR) # Not needed here but too afraid to delete. 'df_ntc' assigned later

# Setup axes ---------------------
# sgrna_pos <- stmn2_sgrna_summary %>% 
#   dplyr::filter(LFC > 0) %>% 
#   dplyr::select(sgrna, Gene, LFC,
#                 pvalue = p.high,
#                 fdr = FDR)
# sgrna_neg <- stmn2_sgrna_summary %>% 
#   dplyr::filter(LFC < 0) %>% 
#   dplyr::select(sgrna, Gene, LFC,
#                 pvalue = p.low,
#                 fdr = FDR)

# Reformat Gene Summary ---------------
# df_pos <- stmn2_gene_summary %>%
#   dplyr::filter(pos.lfc > 0) %>%
#   dplyr::select(Gene,lfc = pos.lfc,
#                 score = pos.score,
#                 pvalue = pos.p.value,
#                 fdr = pos.fdr)
# df_neg <- stmn2_gene_summary %>%
#   dplyr::filter(neg.lfc < 0) %>%
#   dplyr::select(Gene, lfc = neg.lfc,
#                 score = neg.score,
#                 pvalue = neg.p.value,
#                 fdr = neg.fdr)
# stmn2_gene_summary <- rbind(df_neg, df_pos) %>%
#   dplyr::mutate(
#     qvalue = -log(pvalue),
#     hit_direction = ifelse(lfc < 0, "Negative", "Positive")
#   )

# Merge sgRNA summary with gene summary --------------------
# merged_df <- stmn2_gene_summary %>% 
#   left_join(stmn2_sgrna_summary, by = "Gene")

# Visualize distribution of columns with boxplots ------------------------------------
# numeric_columns <- stmn2_sgrna_summary %>% dplyr::select_if(is.numeric)
# metric_columns <- numeric_columns %>% dplyr::select(FDR, LFC, p.high, p.low, p.twosided)
# plot_all_boxplots(numeric_columns)
# plot_all_boxplots(metric_columns)


# Thresholding ---------------------------
# Filter 'raw_counts' by read count threshold
## Setup threshold parameters -----------------------
control_groups <- c("Low_1", "Low_2", "Low_3")
treatment_groups <- c("High_1", "High_2", "High_3")
control_thres <- 300
treatment_thres <- 300
# Assuming control_groups and treatment_groups are vectors of column names.
#     ('raw_counts[, control_groups]') Subset control_groups from treatment_groups,  
#     ('raw_counts[, control_groups] > control_thres') create a boolean matrix if element is TRUE, 
#     ('rowSums(raw_counts[, control_groups] > control_thres)') calculate the sum of control samples exceeding threshold for each row
#     ('rowSums(raw_counts[, control_groups] > control_thres) == length(control_groups)') Check that all 3 replicate samples exceed the threshold
df_thres <- raw_counts[rowSums(raw_counts[, control_groups] > control_thres) == length(control_groups) &
rowSums(raw_counts[, treatment_groups] > treatment_thres) == length(treatment_groups), ]

# Noise modeling ---------------
# TODO: consider methods of empirical noise modeling
# observed_data <- df_thres  # Replace with actual observed data
# estimated_noise = sd(observed_data$LFC - mean(observed_data$LFC))

# Process sgRNA data -----------------------------
ranked_results <- rank_test(stmn2_sgrna_summary)
df <- do.call(rbind, lapply(ranked_results, as.data.frame))
# Rename columns if necessary
colnames(df) <- c('epsilon', 'pvalue')
df$gene <- row.names(df)
row.names(df) <- NULL
df <- df %>% dplyr::select(gene, epsilon, pvalue)
# Extract the gene name from the 'index' column
df$gene_red <- sapply(df$gene, function(x) strsplit(x, '_')[[1]][1])
# Negate the 'epsilon' values
df$epsilon <- -df$epsilon
# Filter the rows where 'gene' is 'negative'
df_ntc <- subset(df, gene_red == 'negative') 
# Calculate the median of 'epsilon' in the negative control group
ntc_median_epsilon <- median(df_ntc$epsilon, na.rm = TRUE)

# p_reduced_df_ntc <- df_ntc %>% dplyr::filter(pvalue <= ntc_pvalue)
# p_reduced_df_target <- df %>% dplyr::filter(pvalue <= 0.083872 & !gene_red == "negative")
# Subtract the median epsilon from the epsilon values
df$sub_ntc_med_epsilon <- df$epsilon - ntc_median_epsilon

# Calculate the 'product' column
df$product <- df$sub_ntc_med_epsilon * (-log10(df$pvalue))
# df$product_untransformed <- df$sub_ntc_med_epsilon * (df$pvalue)

# Sort the data frame by the 'product' column in descending order
df <- df[order(-df$product), ]

# FDR --------------------- 

df_hits <- product_threshold_fdr(df, fdr = 0.05)
thres <- df_hits[[1]]
df_hits <- df_hits[[2]]
df_hits$qvalue <- -log10(df_hits$pvalue)
# Sorting and saving all genes
df_sorted <- df[order(-df$product), ]

# Sorting and saving hits
df_hits_sorted <- df_hits[order(-df_hits$product), ]

# Filtering for rows containing 'NTC' in the 'index' column
# df_ntc <- df[grepl("NTC", df$index), ]

npg <- c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2")

# Separate positive and negative hits
df_pos <- df_hits[df_hits$epsilon > 0, ]
df_neg <- df_hits[df_hits$epsilon < 0, ]

# Create volcano plot ---------
# Create a new column in df to indicate the category of each point
df$category <- "Other genes"
df$category[df$gene %in% df_pos$gene] <- "Positive hits"
df$category[df$gene %in% df_neg$gene] <- "Negative hits"
df$category[df$gene_red %in% 'negative'] <- "Negative control"
# Top genes from "Low" group
low_top_genes <- c('GET4','MAP2K4','CLASRP','POU3F1','ECEL1','NOTCH1','ATF2','CHORDC1','CALM2','UBE3C')
high_top_genes <- c('RANBP1','KIAA1919','FADS2','CCNC','TARDBP','STMN2','ZC3H13','LARP6','RBL2','PAPOLG')
# Use outliers for labeling points on graph
# outliers <- df %>% dplyr::filter(gene %in% high_top_genes)
# outliers <- df %>% dplyr::filter(-log10(pvalue)>1.3 & gene %in% df_pos$gene & !gene %in% df_ntc$gene)

## Plot attempt 1 ----------
# Create scatter plot with legend
p <- ggplot() +
  geom_jitter(data = df, aes(x = sub_ntc_med_epsilon, y = -log10(pvalue), color = category), size = 1.5, alpha = 0.5) +
  geom_jitter(data = df_hits %>% dplyr::filter(gene == "STMN"), aes(x = sub_ntc_med_epsilon, y = -log10(pvalue)), size = 1.5, alpha = 0.7) +
  # Add labels for outliers using geom_label_repel
  # geom_label_repel(data = outliers, aes(x = sub_ntc_med_epsilon, y = -log10(pvalue), label = gene), 
  #                  size = 3, color = "black", box.padding = 0.2, point.padding = 0.2, 
  #                  segment.color = "grey50") +
  # geom_text(data = outliers, aes(x = sub_ntc_med_epsilon, y = -log10(pvalue), label = gene),
  #           vjust = -0.5, hjust = -0.5, angle = 60, size = 3, color = "black") +
  scale_color_manual(values = c("Positive hits" = "red", 
                                "Negative hits" = "#3C5488FF", 
                                "Other genes" = '#F39B7FFF', 
                                "Negative control" = 'grey')) +
  labs(title = "STMN2 Screen sgRNA (FDR=0.05), thresh_25",
       x = "Phenotype", y = "-log10 P", color = "Category") +
  theme_linedraw() +
  theme(legend.position = "right", legend.text = element_text(size = 10))
# Add the FDR threshold lines
x <- seq(0.01, 10, by = 0.01)
y <- sapply(x, function(i) thres / i)

# Adding FDR lines
p <- p + geom_line(aes(x = x, y = y), linetype = "dashed", color = "black") +
  geom_line(aes(x = -x, y = y), linetype = "dashed", color = "black") +
  annotate("text", x = max(x), y = y[1], label = "FDR = 0.05", hjust = 1, size = 5)

# Define limits
lim <- max(abs(min(df$epsilon) - 1), (max(df$epsilon) + 2))
p <- p + xlim(-lim, lim) + ylim(0, -log10(min(df$pvalue)) + 0.5)
print(p)

## Plot attempt 2 ----------
p <- ggplot() +
  # geom_point(data = df %>% dplyr::filter(gene == "STMN2" | gene == "TARDBP"), size = 1.5, alpha = 0.7) +
  geom_jitter(data = df_pos, aes(x = sub_ntc_med_epsilon, y = -log10(pvalue)), color = "red", size = 1.5, alpha = 0.7) +
  geom_jitter(data = df_neg, aes(x = sub_ntc_med_epsilon, y = -log10(pvalue)), color = "#3C5488FF", size = 1.5, alpha = 0.7) +
  # geom_jitter(data = df %>% dplyr::filter(!gene %in% df_pos$gene & !gene %in% df_neg$gene), aes(x = sub_ntc_med_epsilon, y = -log10(pvalue)), color = '#F39B7FFF', size = 1.5, alpha = 0.3) +
  geom_point(data = df_ntc, aes(x = epsilon, y = -log10(pvalue)), color = 'grey', size = 1.5, alpha = 0.1) +
  labs(title = "STMN2 Screen sgRNA (FDR=0.05), thresh_25",
    x = "Phenotype", y = "-log10 P") +
  theme(legend.position = "top", legend.text = element_text(size = 12)) +
  theme_minimal()

# Add the FDR threshold lines
x <- seq(0.01, 10, by = 0.01)
y <- sapply(x, function(i) thres / i)

# Adding FDR lines
p <- p + geom_line(aes(x = x, y = y), linetype = "dashed", color = "black") +
  geom_line(aes(x = -x, y = y), linetype = "dashed", color = "black") +
annotate("text", x = max(x), y = y[1], label = "FDR = 0.05", hjust = 1, size = 5)

# Define limits
lim <- max(abs(min(df$epsilon) - 1), (max(df$epsilon) + 2))
p <- p + xlim(-lim, lim) + ylim(0, -log10(min(df$pvalue)) + 0.5)
print(p)






















###############################################################
###############################################################
###############################################################

# TODO run mageck-iNC steps with normalized counts csv. Currenlty running it on raw counts file


# z <- df_targeting %>% 
#   dplyr::group_by(Gene) %>%
#   filter(n() > 1) %>%  arrange(p.twosided)
# 
# # There is a subset of genes that with multiple sgRNAs. Should test to eval if wilcox output is same as python for this subset
# df <- stmn2_sgrna_summary %>% filter(treat_mean > 300)
# 
# # Split the dataframe into negative controls and targeting genes
# df_ntc <- df %>% filter(grepl('negative_control', Gene))
# df_targeting <- df %>% filter(!grepl('negative_control', Gene))
# 
# # Extract relevant columns from the negative control dataframe
# ntc_sgRNA_p <- df_ntc$p.twosided
# ntc_sgRNA_p_lfc <- data.frame(p.twosided = df_ntc$p.twosided, LFC = df_ntc$LFC)
# 
# # Commented out two lines below: 
# # genes <- unique(df_targeting$Gene)
# # num_of_genes <- length(genes)
# df_targeting <- df_targeting %>% dplyr::group_by(Gene) %>% filter(n() > 1) # added line to subset genes with multiple sgRNAs
# genes <- unique(df_targeting$Gene) # added line to subset genes with multiple sgRNAs
# num_of_genes <- length(genes) # added line to subset genes with multiple sgRNAs
# 
# # Initialize an empty list to store LFC and p-values for each gene
# gene_lfc_p <- list()
# 
# # Loop through each gene and calculate the mean LFC from the top 3 sgRNAs and  calculations
# for (gene in genes) {
#   df_gene <- df_targeting %>% filter(Gene == gene) %>% arrange(p.twosided)
#   
#   # Calculate the mean LFC of up to the top 3 rows
#   lfc <- mean(df_gene$LFC[1:min(3, nrow(df_gene))])
#   
#   # Perform the Wilcoxon rank-sum test
#   # TODO: try with wilcox signed rank test... wilcox rank-sum/mann whitney needs independent samples... but are samples really independent?
#   pvalue <- wilcox.test(df_gene$p.twosided[1:min(3, nrow(df_gene))], ntc_sgRNA_p, alternative = "two.sided", correct = TRUE)$p.value
#   
#   # Store the results in the list
#   gene_lfc_p[[gene]] <- list(lfc = lfc, pvalue = pvalue)
# }

# TODO rerun entire thing with thresh_300_.sgrna_summary




# ATTEMPT TO PLOT FILTERED PLOT EVERYTHING BELOW WORKS. DO NOT CHANGE ---------------------------------------
## rank_test function implemented with single case. trying to recreate plot previously generated -----------
df <- stmn2_sgrna_summary %>% filter(treat_mean > 300)

# Split the dataframe into negative controls and targeting genes
df_ntc <- df %>% filter(grepl('negative_control', Gene))
df_targeting <- df %>% filter(!grepl('negative_control', Gene))

# Extract relevant columns from the negative control dataframe
ntc_sgRNA_p <- df_ntc$p.twosided
ntc_sgRNA_p_lfc <- data.frame(p.twosided = df_ntc$p.twosided, LFC = df_ntc$LFC)

# Get unique genes in the targeting dataframe
genes <- unique(df_targeting$Gene)
num_of_genes <- length(genes)

# Initialize an empty list to store LFC and p-values for each gene
gene_lfc_p <- list()

# Loop through each gene and calculate the mean LFC from the top 3 sgRNAs and  calculations
for (gene in genes) {
  df_gene <- df_targeting %>% filter(Gene == gene) %>% arrange(p.twosided)
  
  # Calculate the mean LFC of up to the top 3 rows
  lfc <- mean(df_gene$LFC[1:min(3, nrow(df_gene))])
  
  # Perform the Wilcoxon rank-sum test
  # TODO: try with wilcox signed rank test... wilcox rank-sum/mann whitney needs independent samples... but are samples really independent?
  pvalue <- wilcox.test(df_gene$p.twosided[1:min(3, nrow(df_gene))], ntc_sgRNA_p, alternative = "two.sided", correct = FALSE, exact = TRUE)$p.value
  
  # Store the results in the list
  gene_lfc_p[[gene]] <- list(lfc = lfc, pvalue = pvalue)
}

# Set the seed for reproducibility
set.seed(10)
# Perform randomization for the negative controls to create pseudo-genes
for (j in seq_len(num_of_genes)) {
  # Shuffle the ntc_sgRNA_p_lfc list
  ntc_sgRNA_p_lfc <- ntc_sgRNA_p_lfc[sample(nrow(ntc_sgRNA_p_lfc)), ]
  ntc_selected <- ntc_sgRNA_p_lfc[1:5, ]
  ntc_selected_p <- ntc_selected$p.twosided
  ntc_lfc <- mean(ntc_selected$LFC[order(ntc_selected_p)][1:min(3, length(ntc_selected_p))])
  ntc_pvalue <- wilcox.test(ntc_selected_p, ntc_sgRNA_p, alternative = "two.sided")$p.value
  gene_lfc_p[[paste0('negative_control_', j)]] <- list(lfc = ntc_lfc, pvalue = ntc_pvalue)
}

## Post-process rank_test output --------------------------------
ranked_results_df <- do.call(rbind, lapply(gene_lfc_p, as.data.frame))
# Rename columns if necessary
colnames(ranked_results_df) <- c('epsilon', 'pvalue')
ranked_results_df$gene <- row.names(ranked_results_df)
row.names(ranked_results_df) <- NULL
ranked_results_df <- ranked_results_df %>% dplyr::select(gene, epsilon, pvalue)
# Extract the gene name from the 'index' column
ranked_results_df$gene_red <- sapply(ranked_results_df$gene, function(x) strsplit(x, '_')[[1]][1])
# Negate the 'epsilon' values
ranked_results_df$epsilon <- -ranked_results_df$epsilon
# Filter the rows where 'gene' is 'negative'
ranked_results_df_ntc <- subset(ranked_results_df, gene_red == 'negative')
# Calculate the median of 'epsilon' in the negative control group
ntc_median_epsilon <- median(ranked_results_df_ntc$epsilon, na.rm = TRUE)
# Filter by pvalue
# p_reduced_df_ntc <- ranked_results_df_ntc %>% dplyr::filter(pvalue <= ntc_pvalue)
# p_reduced_df_target <- ranked_results_df %>% dplyr::filter(pvalue <= 0.08387257 & !gene %in% df_ntc$gene)
# Subtract the median epsilon from the epsilon values
ranked_results_df <- ranked_results_df %>% dplyr::mutate(sub_ntc_med_epsilon = epsilon - ntc_median_epsilon)

# Calculate the 'product' column
ranked_results_df$product <- ranked_results_df$sub_ntc_med_epsilon * (-log10(ranked_results_df$pvalue))

# Sort the data frame by the 'product' column in descending order
ranked_results_df <- ranked_results_df[order(-ranked_results_df$product), ]

## Product threshold function applied to single case. Trying to recreate previouly generated plot --------------------- 
# Find the maximum absolute value in the 'product' column
maxi <- max(abs(df$product), na.rm = FALSE)

# Iterate through the sequence from 0 to maxi with a step of 0.1
for (pro in seq(0, maxi, by = 0.1)) {
  # Filter the dataframe based on the current threshold
  df_thres <- df %>% filter(abs(product) > pro)
  
  # Check if the proportion of 'NTC' entries is below the FDR threshold
  if ((sum(grepl('negative_control', df_thres$gene)) / nrow(df_thres)) < fdr) {
    break
  }
}
# Return the threshold and the filtered dataframe
return(list(thres = pro, df_thres = df_thres))

### FDR ------------------
df_hits <- product_threshold_fdr(ranked_results_df, fdr = 0.05)
thres <- df_hits[[1]]
df_hits <- df_hits[[2]]
df_hits$qvalue <- -log10(df_hits$pvalue)

## Create volcano plot -------------
## Attempt 1 ----------------------
# Create a new column in df to indicate the category of each point
df_hits$category <- "Other genes"
df_hits$category[df_hits$gene %in% df_pos$gene] <- "Positive hits"
df_hits$category[df_hits$gene %in% df_neg$gene] <- "Negative hits"
df_hits$category[df_hits$gene_red %in% 'negative'] <- "Negative control"

df_pos <- df_hits[df_hits$epsilon > 0, ]
df_neg <- df_hits[df_hits$epsilon < 0, ]
p <- ggplot() +
  # geom_point(data = df %>% dplyr::filter(gene == "STMN2" | gene == "TARDBP"), size = 1.5, alpha = 0.7) +
  geom_jitter(data = df_pos, aes(x = sub_ntc_med_epsilon, y = -log10(pvalue)), color = "red", size = 1.5, alpha = 0.7) +
  geom_jitter(data = df_neg, aes(x = sub_ntc_med_epsilon, y = -log10(pvalue)), color = "#3C5488FF", size = 1.5, alpha = 0.7) +
  geom_jitter(data = df_hits %>% dplyr::filter(!gene %in% df_pos$gene & !gene %in% df_neg$gene & !gene %in% df_ntc$gene), aes(x = sub_ntc_med_epsilon, y = -log10(pvalue)), color = '#F39B7FFF', size = 1.5, alpha = 0.5) +
  geom_point(data = ranked_results_df_ntc, aes(x = epsilon, y = -log10(pvalue)), color = 'grey', size = 1.5, alpha = 0.1) +
  labs(title = "STMN2 Screen sgRNA (FDR=0.05)",
       x = "Phenotype", y = "-log10 P") +
  theme(legend.position = "right", legend.text = element_text(size = 12)) +
  theme_minimal()

# Add the FDR threshold lines
x <- seq(0.01, 10, by = 0.01)
y <- sapply(x, function(i) thres / i)
p <- p + geom_line(aes(x = x, y = y), linetype = "dashed", color = "black") +
geom_line(aes(x = -x, y = y), linetype = "dashed", color = "black") +
annotate("text", x = max(x), y = y[1], label = "FDR = 0.05", hjust = 1, size = 5)

# Define limits
lim <- max(abs(min(df_hits$epsilon) - 1), (max(df_hits$epsilon) + 2))
p <- p + xlim(-lim, lim) + ylim(0, -log10(min(df_hits$pvalue)) + 0.5)

# Display the plot
print(p)



















# SCRATCH ------------------------------------------------------
# # Assuming df contains the data and is already created as per the previous code
# 
# # Create a factor to differentiate categories for coloring
# df$category <- "Other genes"
# df$category[df$epsilon > 0] <- "Positive hits"
# df$category[df$epsilon < 0] <- "Negative hits"
# 
# # Assign colors to each category
# colors <- c("Positive hits" = "#DC0000FF", "Negative hits" = "#3C5488FF", "Other genes" = "#F39B7FFF")
# 
# # Create the plot
# p <- ggplot(df, aes(x = epsilon, y = -log10(pvalue), color = category)) +
#   geom_point(size = 1.5, alpha = 0.7) +
#   scale_color_manual(values = colors) +
#   labs(x = "Phenotype", y = "-log10 P", color = "Category") +
#   theme_minimal() +
#   theme(legend.position = "top", legend.text = element_text(size = 12))
# 
# # Add the FDR threshold lines
# x <- seq(0.01, 10, by = 0.01)
# y <- sapply(x, function(i) thres / i)
# # p <- p + geom_line(aes(x = x, y = y), linetype = "dashed", color = "black") +
#   # geom_line(aes(x = -x, y = y), linetype = "dashed", color = "black") +
#   # annotate("text", x = max(x), y = y[1], label = "FDR = 0.05", hjust = 1, size = 5)
# 
# # Define limits
# lim <- max(abs(min(df$epsilon) - 1), (max(df$epsilon) + 2))
# p <- p + xlim(-lim, lim) + ylim(0, -log10(min(df$pvalue)) + 0.5)
# 
# # Display the plot with legend
# print(p)



## Alternative rank_test function. Avoids creating a for loop and uses dplyr/magrittr. Have not tested if it works -----
# rank_test <- function(df) {
#   # Filter rows where 'treat_mean' is greater than 20
#   df <- df %>% filter(treat_mean > 20)
#   
#   # Split the dataframe into negative controls and targeting genes
#   df_ntc <- df %>% filter(grepl('negative_control', Gene))
#   df_targeting <- df %>% filter(!grepl('negative_control', Gene))
#   
#   # Extract relevant columns from the negative control dataframe
#   ntc_sgRNA_p <- df_ntc$p.twosided
#   
#   # Define a function to calculate LFC and perform the Wilcoxon test for each gene
#   process_gene <- function(df_gene) {
#     df_gene <- df_gene %>% arrange(p.twosided)
#     top_lfc <- df_gene$LFC[1:min(3, nrow(df_gene))]
#     lfc <- mean(top_lfc)
#     pvalue <- wilcox.test(top_lfc, ntc_sgRNA_p, alternative = "two.sided")$p.value
#     return(c(lfc = lfc, pvalue = pvalue))
#   }
#   
#   # Apply the process_gene function to each gene group
#   results <- df_targeting %>% 
#     group_by(Gene) %>% 
#     summarize(lfc_pvalue = list(process_gene(cur_data()))) %>% 
#     unnest_wider(lfc_pvalue)
#   
#   # Convert the results to a named list, similar to the original output
#   gene_lfc_p <- split(results[-1], results$Gene)
#   
#   return(gene_lfc_p)
# }



# # Plotting all genes ------------------------
# ggplot() +
#   geom_point(data = both %>% dplyr::filter(id == "STMN2" | id == "TARDBP"), color = "red", size = 2)
#   geom_point(data = stmn2_gene_summary, aes(x = lfc, y = -log10(pvalue)), color = "grey", alpha = 0.3, size = 1) +
#   geom_point(data = df_pos, aes(x = lfc, y = -log10(pvalue)), color = "thistle", alpha = 0.5, size = 1) +
#   geom_point(data = df_neg, aes(x = lfc, y = -log10(pvalue)), color = "lightblue", alpha = 0.5, size = 1) +
#   geom_line(aes(x = seq(0.01, 50, 0.01), y = 0.5 / seq(0.01, 50, 0.01)), linetype = "dotted", color = "black") +
#   geom_line(aes(x = -seq(0.01, 50, 0.01), y = 0.5 / seq(0.01, 50, 0.01)), linetype = "dotted", color = "black") +
#   geom_point(data = df_ntc, aes(x = LFC, y = -log10(pvalue)), color = "black", alpha = 0.3, size = 1) +
#   xlim(-max(abs(min(stmn2_gene_summary$lfc) - 1), max(stmn2_gene_summary$lfc) + 2), 
#        max(abs(min(stmn2_gene_summary$lfc) - 1), max(stmn2_gene_summary$lfc) + 2)) +
#   ylim(0, -log10(min(stmn2_gene_summary$pvalue)) + 0.5) +
#   labs(x = "Normalized Phenotype", y = "-log10 P", title = "STMN2") +
#   theme_classic() +
#   theme(text = element_text(family = "Arial", size = 10))
# 
# 
# 
# 
# 


