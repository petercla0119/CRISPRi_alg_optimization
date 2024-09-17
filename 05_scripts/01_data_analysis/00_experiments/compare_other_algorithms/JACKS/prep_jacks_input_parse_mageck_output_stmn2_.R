# Parse rra-High-Low-STMN2.sgrna_summary.txt file into a counts file for JACKS pipeline

library(tidyverse)
library(dplyr)
library(readxl)

rra.High.Low.STMN2.sgrna_summary <- read.delim("/Users/Claire/Library/CloudStorage/Box-Box/MSTP Research Rotations/Summer_2024_Might_Ward_Lab/CRISPRi_screenfiles/rra-High-Low-STMN2.sgrna_summary.txt")

# Separate 
# Separate the 'combined_values' column into three new columns
counts_sep <- rra.High.Low.STMN2.sgrna_summary %>%
  dplyr::select(sgrna, Gene, control_count, treatment_count) %>% 
  separate(control_count, into = c("low_rep1", "low_rep2", "low_rep3"), sep = "/") %>% 
  separate(treatment_count, into = c("high_rep1", "high_rep2", "high_rep3"), sep = "/")

# Replace negative_control with CTRL
counts_sep[counts_sep=="negative_control"] <- "CONTROL"

# Print the modified data frame
print("Data Frame after Separation:")
print(head(counts_sep))

output_file <- "counts_sep_CONTROL.txt"
# Save the modified data frame to a .txt file
write.table(counts_sep, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

###########################################
# Create replicatemapfile for JACKS pipeline
column_names <- colnames(counts_sep)
-t(colnames(counts_sep))
# Create a new data frame with column names as rows
column_names_df <- data.frame(Column_Names = column_names, stringsAsFactors = FALSE)
############################################

#create sgrnamappingfile
sgrnamapping <- counts_sep %>% 
  dplyr::select(sgrna, Gene)
sgrnamapping_file <- "sgrnamapping.txt"
# Save the modified data frame to a .txt file
write.table(sgrnamapping, file = sgrnamapping_file, sep = "\t", row.names = FALSE, quote = FALSE)

# create file for ctrl_gene - used as negative control genes

# Import sgRNA library file
sgrna_lib <- read_excel("~/Downloads/elife-81856-supp4-v2.xlsx", sheet = "library 1+2")
# gene and essentiality columns are coded as characters
sgrna_lib$gene_source <- as.factor(sgrna_lib$gene_source) # recode essentiality col as factor to easily view options
# view levels - should be four levels
levels(sgrna_lib$gene_source)
# Create negative control file with non-essential genes and negative_controls
sgrna_lib %>% 
  dplyr::filter(gene_source=="nonessential")
