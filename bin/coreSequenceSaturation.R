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

set.seed(56)

# Functions ---------------------------------------------------------------

# Function form calculating saturation metrics over multiple data frames
subsampler <- function(amount_sample, data, total_gene_counts) {
    # this grabs a number of random indices (gene counts) between 0 and
    # the total gene counts equal to the current level of downsampling
    draw_indices <- sample.int(total_gene_counts, amount_sample)

    # this takes the indices we previously had and maps them back
    # onto the buckets for the genes/cell
    items <- findInterval(draw_indices - 1,
        c(0, cumsum(data$support)),
        rightmost.closed = TRUE)

    # This creates a list of the sorted unique buckets of all
    # genes per cell where at least one count was seen
    unique_gene_per_cell <- sort(unique(items))

    # initialize new empty column which will tell us if a gene is observed
    # in a cell after downsampling
    data[, observed := 0]

    # if the gene was observed during downsampling set its value to 1
    data$observed[unique_gene_per_cell] <- 1

    # groups by barcode and counts the total number of unique genes observed per barcode
    data %>%
        select(barcode, observed) %>%
        group_by(barcode) %>%
        summarize(total = sum(observed)) %>%
        pull(total) %>%
        median()

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
    total_reads <- metric_file %>%
        filter(sample == sample_id
               & process == "fastqc"
               & metric == "tot_seq") %>%
        summarize(total_reads = sum(value)) %>%
        pull(total_reads) %>%
        as.integer()

    # Calculate mean reads per cell
    read_pairs <- total_reads / 2
    mean_read_pairs_per_cell <- read_pairs / length(unique(fragment_df$barcode))
}

# Calculate subsampling intervals for plotting (mean reads per cell)
mean_depth_per_cell <- round(seq(0, mean_read_pairs_per_cell, length.out = 25))

# Get the total size of the pool of gene counts we are sampling from
total_gene_counts <- sum(fragment_df$support)

# Calculate subsampling intervals for downsampling gene counts
subsample_intervals_downsample <- round(seq(0, total_gene_counts, length.out = 25))

# Calculate saturation curve
# x-axis: mean reads per cell
# y-axis: median number of unique fragments per cell
saturation_df <- as.data.frame(subsample_intervals_downsample) %>%
    rowwise %>%
    mutate(median_obs_per_cell = subsampler(
        amount_sample = subsample_intervals_downsample,
        data = fragment_df,
        total_gene_counts = total_gene_counts)) %>%
    add_column(mean_depth_per_cell = mean_depth_per_cell) %>%
    select(mean_depth_per_cell, median_obs_per_cell)

# Write out saturation curve coordinates
output_file_name <- paste0(sample_id, "_sequence_saturation.csv")
write_csv(saturation_df, file = output_file_name)
