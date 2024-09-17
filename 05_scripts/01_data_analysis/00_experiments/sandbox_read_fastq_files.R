#######################################
# IGNORE LINES BELOW - USED AS SANDBOX
#######################################

library(ShortRead)
get_read_counts <- function(fastq_file) {
  # Read the .fastq file
  fq <- readFastq(fastq_file)
  
  # Extract sequences
  seqs <- sread(fq)
  
  # Convert to character vector
  seqs_char <- as.character(seqs)
  
  # Count unique sequences
  read_counts <- table(seqs_char)
  
  # Convert to data frame
  read_counts_df <- as.data.frame(read_counts)
  
  # Rename columns for clarity
  colnames(read_counts_df) <- c("Sequence", "Count")
  
  return(read_counts_df)
}

# Path to your .fastq file
fastq_file <- "~/Downloads/hits.JH8105_1_S1_L001_R1_001.fastq"

# Get read counts
read_counts <- get_read_counts(fastq_file)

# Print the read counts
print(read_counts)

# Optionally, save to a CSV file
write.csv(read_counts, "read_counts.csv", row.names = FALSE)


ShortRead::readFastq(fastq_file)
df=ShortRead::readFastq(fastq_file)

# df2=ShortRead::readFastq(fastq_file) %>% mutate(Length=str_length(Sequence)) %>% filter(Length>200) %>% ShortRead::writeFasta(out.file="long_reads.fa")
library(dplyr)
library(ggplot2)
library(Biostrings)
summary(df)fastq_file <- "~/Downloads/hits.JH8105_1_S1_L001_R1_001.fastq"
# Get read counts

df2=ShortRead::readFastq(fastq_file) %>% mutate(Length=str_length(Sequence)) %>% filter(Length>200) %>% ShortRead::writeFasta(out.file="long_reads.fa")
head(df)
reads=ShortRead::sread(df)
reads@ranges@width
widths=reads@ranges@width
widths=as.data.frame(reads@ranges@width)
ggplot(widths) +
  geom_histogram(aes(x=reads@ranges@width))
quals=Biostrings::quality(df)
