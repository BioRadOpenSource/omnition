#!/usr/bin/env Rscript

library(tidyverse)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-bc", "--bead_counts", help = "Path to bead counts file.")
parser$add_argument("-a", "--allowlist", help = "Path to allowlist.")
parser$add_argument("-c", "--cell_expression", help = "Path to cell expression file.")
parser$add_argument("-b", "--bead_expression", help = "Path to bead expression file.")
parser$add_argument("-bt", "--barcode_translate",
                    help = "Path to barcode translation file.")
parser$add_argument("-s", "--sample", help = "Sample name.")

parser <- parser$parse_args()

bead_counts_file <- parser$bead_counts
allowlist_file <- parser$allowlist
cell_expr_file <- parser$cell_expression
bead_expr_file <- parser$bead_expression
barcode_translate_file <- parser$barcode_translate
sample_name <- parser$sample

# Function calls
assemble_bead_to_cell_table <- function(
    bead_counts,
    allowlist,
    cell_expr,
    bead_expr,
    barcode_translate
) {
    # assemble bead-to-cell table
    bead_cell_data <-
        # add cell barcode for each bead
        left_join(bead_expr, barcode_translate, by = c("barcode" = "BeadBarcode")) %>%
        # reorganize columns for readability
        relocate(DropBarcode, .before = "barcode") %>%
        # add in cell barcode data; this creates duplicate values for bead multiplets
        left_join(cell_expr, by = c("DropBarcode" = "barcode")) %>%
        # rename columns to denote bead or cell as origin of the count
        rename_with(~ gsub(".x", ".bead", .x, fixed = TRUE), ends_with(".x")) %>%
        rename_with(~ gsub(".y", ".cell", .x, fixed = TRUE), ends_with(".y")) %>%
        # add in total aligned (with duplicates) and unaligned read counts
        left_join(., bead_counts, by = "barcode") %>%
        relocate(unalignedReads, .before = "input_reads.bead") %>%
        # move columns around
        relocate(alignedReads, .before = "unalignedReads") %>%
        group_by(DropBarcode) %>%
        arrange(-umi.cell) %>%
        relocate(nBeads, .after = last_col()) %>%
        # add boolean for whether the cell barcode was above the knee
        mutate(cell_above_knee = case_when(
            is.na(DropBarcode) ~ barcode %in% allowlist$CellBarcode,
            TRUE ~ DropBarcode %in% allowlist$CellBarcode)) %>%
        # rename columns for readability and consistency
        rename(bead_barcode = barcode) %>%
        rename(cell_barcode = DropBarcode) %>%
        rename(total_aligned_reads = alignedReads) %>%
        rename(unaligned_reads = unalignedReads) %>%
        rename(deduplicated_reads.bead = input_reads.bead) %>%
        rename(deduplicated_reads.cell = input_reads.cell) %>%
        rename(n_beads = nBeads)
    return(bead_cell_data)
}

calc_bead_merging_stats <- function(
    bead_cell_data,
    sample_name
) {
    # Calculate bead merging metrics
    n_cells <- nrow(bead_cell_data %>%
        filter(cell_above_knee) %>%
        distinct(cell_barcode))
    cells_with_merges <- nrow(bead_cell_data %>%
        filter(cell_above_knee) %>%
        filter(n_beads > 1) %>%
        distinct(cell_barcode))
    pct_cells_with_multiple_beads <- round(cells_with_merges / n_cells * 100, 2)

    bead_merging_metrics <- data.frame(
        sample = c(sample_name),
        process = c("bead_merging"),
        metric = c(
            "n_cells",
            "cells_with_merges",
            "percent_cells_with_multiple_beads"
        ),
        value = c(
            n_cells,
            cells_with_merges,
            pct_cells_with_multiple_beads
        )
    )
    return(bead_merging_metrics)
}

calc_sum_expr_metrics <- function(
    cell_expr,
    allowlist,
    sample_name
) {
    # Summarize expression metrics
    cell_expr <- cell_expr %>%
        rename(deduplicated_reads = input_reads)
    cell_expr_sum <- cell_expr %>%
        select(-barcode) %>%
        summarise(across(everything(), sum, na.rm = TRUE)) %>%
        select(-genes)

    # If we have columns starting with "input", that means it was a mixed run
    # Need to calculate the number of reads for each species in cells
    if (any(grepl("^input", names(cell_expr)))) {
        mixed_in_cells_sum <- cell_expr %>%
            filter(barcode %in% allowlist$CellBarcode) %>%
            # if the column starts with input, add it
            select(starts_with("input")) %>%
            summarise(across(everything(), sum, na.rm = TRUE)) %>%
            # rename to add _in_cell to the end of the existing name
            rename_with(~ paste0(.x, "_in_cells"), everything()) %>%
            # Drop the input_ prefix from the column names
            rename_with(~ gsub("^input_", "", .x), everything())
        mixed_total_sum <- cell_expr %>%
            # if the column starts with input, add it
            select(starts_with("input")) %>%
            summarise(across(everything(), sum, na.rm = TRUE))
        species <- gsub("_reads_in_cells$", "", names(mixed_in_cells_sum))
        # Calculate fraction of reads in cells by species
        mixed_fraction <- data.frame(
            sample = c(paste0(sample_name, "_", species[1]),
                       paste0(sample_name, "_", species[2])),
            process = c("summarize_expression"),
            metric = c("fraction_reads_in_cells"),
            value = c(
                mixed_in_cells_sum[[1, 1]] / mixed_total_sum[[1, 1]] * 100,
                mixed_in_cells_sum[[1, 2]] / mixed_total_sum[[1, 2]] * 100
            )
        )
        # Capture the reads in cells by species
        cell_expr_sum <- cell_expr_sum %>%
            bind_cols(mixed_in_cells_sum)
    }
    # Calculate percentages for relevant metrics
    pct_df <- cell_expr_sum %>%
        select(-deduplicated_reads, -umi, -transcripts) %>%
        # If they exist, drop columns that start with "input"
        select(-starts_with("input")) %>%
        mutate(across(everything(), ~ .x / cell_expr_sum$deduplicated_reads * 100)) %>%
        rename_with(~ paste0("pct_", .x), everything())
    # Combine the two data frames and calculate fraction of reads in cells
    # Reads mapped transcriptome named that way for the report
    # Redundant, but for consistency with the report
    cell_expr_sum <- cell_expr_sum %>%
        bind_cols(pct_df) %>%
        mutate(fraction_reads_in_cells = sum(
            cell_expr$deduplicated_reads[cell_expr$barcode %in% allowlist$CellBarcode]
        ) / cell_expr_sum$deduplicated_reads * 100) %>%
        mutate(reads_mapped_transcriptome = pct_genic_reads)
    # Pivot and format for metric summary
    cell_expr_sum_long <- cell_expr_sum %>%
        pivot_longer(cols = everything(), names_to = "metric", values_to = "value") %>%
        mutate(sample = sample_name, process = "summarize_expression") %>%
        select(sample, process, metric, value)
    if (exists("mixed_fraction")) {
        cell_expr_sum_long <- cell_expr_sum_long %>%
            bind_rows(mixed_fraction)
    }
    return(cell_expr_sum_long)
}

# read input data frames
allowlist <- read_csv(allowlist_file, col_names = "CellBarcode")
bead_counts <- read_csv(bead_counts_file)
cell_expr <- read_csv(cell_expr_file)
bead_expr <- read_csv(bead_expr_file)
barcode_translate <- read_csv(barcode_translate_file)

bead_cell_data <- assemble_bead_to_cell_table(
    bead_counts,
    allowlist,
    cell_expr,
    bead_expr,
    barcode_translate
)

write_csv(bead_cell_data, paste0(sample_name, ".bead_summary.csv"))

bead_merging_metrics <- calc_bead_merging_stats(
    bead_cell_data,
    sample_name
)

sum_exp_metrics <- calc_sum_expr_metrics(
    cell_expr,
    allowlist,
    sample_name
)
all_metrics_df <- bind_rows(bead_merging_metrics, sum_exp_metrics)

write_csv(
    all_metrics_df, paste0(sample_name, ".bead_merge_sum_exp_metrics.csv")
)
