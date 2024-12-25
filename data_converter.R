# Script to fetch series matrix from GEO (GSE64810) and convert to CSV.
# Also, "cleans up" data for input in app.
# I also converted the normalized counts and DESeq results here.

library(tidyverse)
library(readr)
library(GEOquery)
library(stringr)

# getting series matrix
gse <- getGEO('GSE64810',GSEMatrix=TRUE)

# saving matrix as csv
write.csv(gse, file = "GSE64810_series_matrix.csv")

# reading in matrix csv
gse_csv <- read_delim("GSE64810_series_matrix.csv")

# converting data to tibble
gse_tib <- as_tibble(gse_csv) %>%
  # selecting for relevant columns
  select(2:16) %>%
  # removing GSE64810 prefix from column headers
  rename_with(~ gsub("GSE64810_series_matrix.txt.gz.", "", .)) %>%
  # renaming column headers with cleaner names
  rename(
    source = 8,
    organism = 9,
    tissue = 10,
    diagnosis = 11,
    post_mortem_interval = 12,
    age_of_death = 13,
    RNA_integrity_number = 14,
    mRNAseq_reads = 15
  ) %>%
  # removing redundant strings from rows in specific columns
  mutate(
    tissue = str_replace(tissue, "tissue: ", ""),
    diagnosis = str_replace(diagnosis, "diagnosis: ", ""),
    post_mortem_interval = str_replace(post_mortem_interval, "pmi: ", ""),
    age_of_death = str_replace(age_of_death, "age of death: ", ""),
    RNA_integrity_number = str_replace(RNA_integrity_number, "rin: ", ""),
    mRNAseq_reads = str_replace(mRNAseq_reads, "mrna-seq reads: ", "")
  )

write_csv(gse_tib, "data/GSE64810_samples_input.csv")


# converting normalized counts txt file to csv
counts <- read_delim("data/GSE64810_mlhd_DESeq2_norm_counts_adjust.txt")
write_csv(counts, "data/GSE64810_counts_input.csv")

# converting DESeq results txt file to csv
deseq <- read_delim("data/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt")
write_csv(deseq,"data/GSE64810_DESeq2_input.csv")

max(deseq$padj) # 0.9999222
min(deseq$padj) # 2.687232e-34
