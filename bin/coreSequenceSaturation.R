#!/usr/bin/env Rscript
# atacSequenceSaturation.R
# Bio-Rad Laboratories, Inc.

# Purpose: Perform sequence saturation sampling
# Usage: sequenceSaturation.R  sample_id fragments_file summary_file assay



# Dependencies ------------------------------------------------------------

library(tidyverse)
library(data.table)



# Environment -------------------------------------------------------------

# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
sample_id <- args[1]
fragments_file <- args[2]
summary_file <- args[3]
assay <- args[4]


# Functions ---------------------------------------------------------------

# Function for calculating the unique observations at a subsampling depth
rarefy <- function(objects_int_vec, choose_int) {

    objects_int_vec <- objects_int_vec[objects_int_vec > 0]
    n <- sum(objects_int_vec)
    p1 <- ifelse(n - objects_int_vec < choose_int,
                 0,
                 exp(lchoose(n - objects_int_vec, choose_int)
                     - lchoose(n, choose_int)))
    obs <- sum(1 - p1)
    names(obs) <- choose_int

    return(obs)

}

# Function for calculating saturation metrics over multiple data frames
map_rarefy <- function(data, choose_int_vector) {

    map(choose_int_vector, rarefy, objects_int_vec = data$support) %>%
        unlist() %>%
        enframe(name = "mean_depth_per_cell", value = "obs") %>%
        mutate(mean_depth_per_cell = as.numeric(mean_depth_per_cell))

}



# Analysis ----------------------------------------------------------------

if (assay == "ATAC") {
    # Import sample-level fragments file
    fragment_df <- fread(fragments_file,
                        col.names = c("barcode", "support"),
                        select = c(4:5))

    # Import sample-level alignment summary file
    summary_df <- fread(summary_file,
                        skip = 6,
                        nrows = 3,
                        select = c("CATEGORY", "TOTAL_READS")) %>%
        rename_all(tolower) %>%
        mutate(category = tolower(category))

    # Calculate mean reads per cell and subsampling intervals
    read_pairs <- summary_df[category == "pair", total_reads] / 2
    mean_read_pairs_per_cell <- read_pairs / length(unique(fragment_df$barcode))

} else {
    # Import cell-barcode level counts file and format to match ATAC frag file
    fragment_df <- fread(fragments_file,
                         select = c(1, 3),
                        col.names = c("barcode", "support"))

    # Import total read count
    metric_file <- read_csv(summary_file, col_types = cols())
    total_reads <- as.integer(metric_file[which(metric_file$sample == sample_id &&
                                     metric_file$metric == "tot_seq"), 7])

    # Calculate mean reads per cell
    read_pairs <- total_reads / 2
    mean_read_pairs_per_cell <- read_pairs / length(unique(fragment_df$barcode))
}

# Calculate subsampling intervals
subsample_intervals <- round(seq(0, mean_read_pairs_per_cell, length.out = 100))

# Calculate saturation curve
# x-axis: mean reads per cell
# y-axis: median number of unique fragments per cell
saturation_df <- fragment_df %>%
    group_by(barcode) %>%
    nest() %>%
    ungroup() %>%
    mutate(barcode = log(row_number())) %>%
    mutate(saturation = map(data,
                            map_rarefy,
                            choose_int_vector = subsample_intervals)) %>%
    select(saturation) %>%
    unnest(saturation) %>%
    group_by(mean_depth_per_cell) %>%
    summarize(median_obs_per_cell = median(obs, na.rm = TRUE))

# Write out saturation curve coordinates
output_file_name <- paste0(sample_id, "_sequence_saturation.csv")
write_csv(saturation_df, file = output_file_name)
