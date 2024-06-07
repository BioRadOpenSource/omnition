#!/usr/bin/env Rscript
# Bio-Rad Laboratories, Inc.

# Purpose: Count the number of features in the final expression matrix

# Setting environment -----------------------------------
# --------------------------------------

# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
sampleId <- args[1]
filteredMatrix <- args[2]

# Checking if inputs are set
if (is.na(sampleId)) {
  stop("Sample ID must be provided.")
} else if (!file.exists(filteredMatrix)) {
  stop(paste(filteredMatrix, "does not exist."))
}



# Loading dependencies -----------------------------------
# -------------------------------------

library(Matrix)
library(tidyverse)



# Analysis ------------------------------------------------
# ------------------------------------

cat("[PROGRESS]: Importing matrix.\n")

# Reading in sparse matrix
matrix <-  readMM(filteredMatrix)


cat("[PROGRESS]: Counting features.\n")

# Total number of cells in sample (should match numcells results)
estimated_cells <- ncol(matrix)

# Total number of unique features sample
total_unique_features <- sum(rowSums(matrix) > 0)

# Number of unique features per cell
#nolint start
median_unique_features_per_cell <- median(colSums(matrix > 0))
#nolint end

# Total UMI counts for sample
total_umi_counts <- sum(matrix)

# Number of UMI counts per cell
median_umi_counts_per_cell <- median(colSums(matrix))

# Combining into tidy df
feature_counts <- tibble(
  sample = sampleId,
  estimated_cells,
  total_unique_features,
  median_unique_features_per_cell,
  total_umi_counts,
  median_umi_counts_per_cell
) %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value")


cat("[PROGRESS]: Writing output.\n")

# Outputting as CSV
write_csv(feature_counts, paste0(sampleId, "_matrix_features.csv"))
