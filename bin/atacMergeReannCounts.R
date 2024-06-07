#!/usr/bin/env Rscript
# Bio-Rad Laboratories, Inc.

# Purpose: Merge

# Setting environment -------------------------------------
# ------------------------------------


# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
inputDir <- args[1]


# Checking if inputs are set
if (!dir.exists(inputDir)) {
  stop(paste(inputDir, "does not exist."))
}

library(tidyverse)

file_pattern <- "_reannotate_bam_read_counts\\.csv"
counts_files <-
  list.files(path = inputDir,
             pattern = file_pattern,
             full.names = TRUE)
summarized_counts_df <-
  map_dfr(counts_files, read_csv, col_types = cols()) %>%
  group_by(sample, process, metric) %>%
  summarize(count = sum(count))

write.csv(
  summarized_counts_df,
  paste0(
    inputDir,
    "/",
    summarized_counts_df$sample[1],
    "_merged_reannotate_read_counts.csv"
  ),
  row.names = FALSE
)
