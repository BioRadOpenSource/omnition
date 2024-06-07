#!/usr/bin/env Rscript
# coreSequenceSaturation.R
# Bio-Rad Laboratories, Inc.

# Purpose: Perform sequence saturation sampling
# Usage: coreSequenceSaturation.R  \
#   sample_id feature_counts_file metrics_file allowlist_file symbols_file \
#   species_id cDNA_counts assay

# Dependencies ------------------------------------------------------------

packages <- c("data.table", "magrittr", "dqrng", "parallel", "argparse")
invisible(lapply(packages, library, character.only = TRUE))

# Environment -------------------------------------------------------------

# Parsing command line arguments
parser <- ArgumentParser()

# nolint start
parser$add_argument("-id", "--sample_id", type = "character", default = NULL, help = "Sample name")
parser$add_argument("-c", "--feature_counts_file", type = "character", default = NULL, help = "Cell by gene counts file")
parser$add_argument("-mf", "--metrics_file", type = "character", default = NULL, help = "Metric summary file")
parser$add_argument("-af", "--allowlist_file", type = "character", default = NULL, help = "Allow list file for droplets above knee")
parser$add_argument("-sf", "--symbols_file", type = "character", default = NULL, help = "File with gene symbols")
parser$add_argument("-si", "--species_id", type = "character", default = NULL, help = "Species id (eg. Homo_sapiens)")
parser$add_argument("-cc", "--cDNA_counts", type = "character", default = FALSE, help = "Value to use for calculating read depth. If false, total reads from metric summary used.")
parser$add_argument("-a", "--assay", type = "character", required = TRUE, help = "Must be one of the following: ATAC, cATAC, RNA")
# nolint end
parser <- parser$parse_args()

sample_id <- parser$sample_id
feature_counts_file <- parser$feature_counts_file
metrics_file <- parser$metrics_file
allowlist_file <- parser$allowlist_file
symbols_file <- parser$symbols_file
species_id <- parser$species_id
assay <- parser$assay

cDNA_counts <- ifelse(
    is.null(parser$cDNA_counts) | parser$cDNA_counts == "false", FALSE,
    as.integer(parser$cDNA_counts)
)

set.seed(56)

# Set data.table threading
setDTthreads(0)

# Functions ---------------------------------------------------------------

subsampler <- function(amount_sample, feature_dt, total_gene_counts, interval_form) {
    # Function form calculating saturation metrics over multiple data frames

    # We must sample without replacement and keep track of the cell, gene, & count
    # We could paste the cell barcode & gene name together,
    #   and give every cell-gene combination as many rows as it has counts
    # Then we could sample a single vector/column.
    # However, it is very inefficient to have a row for every count
    # So instead of a row for every single count,
    #   we can have a single range to represent the number of rows
    # If a cell-gene has 10 counts, it would have a 1:10 range,
    #   if the next cell has 5 counts, it would have a 11:15 range

    # To achieve this, cumulative sum is used on the gene count column
    #   so that every row becomes an interval along the range of total gene counts
    # Now if we randomly select along the range of total gene counts,
    #    we can find the interval (cell-gene combination) that contains that number
    # This is equivalent to giving every count a row
    # This also retains the order of the orginal DF
    # So that when we sample from the range of total gene counts, and find the interval
    #   that contains that #, we can use it as the row # of the DT with cell-gene info

    # The steps are as follows:
    # 1. Randomly sample from the range of total gene counts
    # 2. Find the interval that contains the sampled number
    # 3. Use the interval row # as the row # of the table with cell-gene info
    # 4. Count the number of unique genes per cell
    # 5. Return the median number of unique genes per cell

    # This creates a vector of random numbers between 1 and the total number of counts
    # So each draw will be a single count for a particular cell-gene
    draw_indices <- dqsample.int(total_gene_counts, amount_sample)

    # This finds the intervals that contains the random number draws
    # where each interval represents a cell-gene and each draw a single count
    items <- findInterval(draw_indices - 1,
        interval_form,
        rightmost.closed = TRUE
    )

    # We only care if the gene was observed in a cell, not how many times
    gene_indicies_observed <- unique(items)

    # items contain row indicies of feature_dt
    # I want to count how many times each row index appears
    # and store that count in a new column called observed_count
    # I can then use that column to filter out rows that were not observed

    # initialize new empty column which will tell us if a gene is observed
    # in a cell after downsampling
    feature_dt[, observed := FALSE]

    # Mark cell-genes that were observed
    feature_dt[gene_indicies_observed, observed := TRUE]

    # Groups by barcode (cell), can counts how many genes were observed per cell
    genes_per_cell <- feature_dt[,
        .(gene_count = sum(observed)),
        by = .(barcode)
    ]

    stats_table <- data.table(
        amount_sampled = amount_sample,
        median_obs_per_cell = median(genes_per_cell$gene_count)
    )

    return(stats_table)
}

# Analysis ----------------------------------------------------------------

if (assay == "ATAC" || assay == "cATAC") {
    # Import feature counts file
    feature_dt <- fread(feature_counts_file,
        select = c(4:5),
        col.names = c("barcode", "count")
    )

    # Import sample-level alignment summary file
    summary_dt <- fread(metrics_file,
        skip = 6,
        nrows = 3,
        select = c("CATEGORY", "TOTAL_READS")
    )
    setnames(summary_dt, names(summary_dt), tolower(names(summary_dt)))
    summary_dt[, category := tolower(category)]

    # Calculate mean reads per cell and subsampling intervals
    read_pairs <- summary_dt[category == "pair", total_reads] / 2
} else {
    # Import cell level feature counts file
    feature_dt <- fread(feature_counts_file,
        select = c(1, 2, 3),
        col.names = c("barcode", "gene", "count")
    )
    # Import metric file to get total read count
    summary_dt <- fread(metrics_file)

    # Select total read count for sample unless cDNA_counts is provided
    if (cDNA_counts == FALSE) {
        read_pairs <- summary_dt[
            sample == sample_id &
                process == "fastqc" &
                metric == "read_count_r1", value
        ]
    } else {
        read_pairs <- cDNA_counts
    }
}

# Calculate mean reads per cell
mean_read_pairs_per_cell <- read_pairs / uniqueN(feature_dt$barcode)

# Calculate subsampling intervals for plotting (mean reads per cell)
mean_depth_per_cell <- seq(0, mean_read_pairs_per_cell, length.out = 25) %>%
    round(.)

# Remove barcodes & gene symbols not found in allowlist & gene symbols file for mixed
if (!is.null(species_id)) {
    allowlist_df <- fread(allowlist_file, select = c(1), col.names = c("barcode"))
    symbols_df <- fread(symbols_file, select = c(1), col.names = c("gene"))

    feature_dt <- feature_dt[
        barcode %in% allowlist_df$barcode &
            gene %in% symbols_df$gene,
    ]
}

# Get the total size of the pool of gene counts we are sampling from
total_gene_counts <- sum(feature_dt$count)

# Calculate subsampling intervals for downsampling gene counts
subsample_intervals_downsample <- seq(0, total_gene_counts, length.out = 25) %>%
    round(.)

# This makes intervals for each cell-gene count along the range of total gene count
# So the range 1:total_gene_counts can now represents 'rows' within these intervals
# FROM: cell gene count
# Cell-A    Gene-X  276
# Cell-A    Gene-Z  73
# Cell-B    Gene-X  253
interval_form <- c(0, cumsum(feature_dt$count)) %>% as.integer()
# TO: interval form
# Cell-A    Gene-X  276
# Cell-A    Gene-Z  349
# Cell-B    Gene-X  602

# Calculate median genes per cell for each subsampling
subsample_stat_table <- lapply(
    subsample_intervals_downsample,
    FUN = subsampler,
    feature_dt = feature_dt,
    total_gene_counts = total_gene_counts,
    interval_form = interval_form
) %>% rbindlist(.)

# Output csv file with saturation curve coordinates
# x-axis: mean reads per cell
# y-axis: median number of unique fragments per cell
output_file <- cbind(
    mean_depth_per_cell,
    subsample_stat_table[, c("median_obs_per_cell")]
)

if (is.null(species_id)) {
    output_file_name <- paste0(sample_id, "_sequence_saturation.csv")
} else {
    output_file_name <- paste0(sample_id, ".", species_id, "_sequence_saturation.csv")
}

# Write out saturation curve coordinates
fwrite(output_file, output_file_name)
