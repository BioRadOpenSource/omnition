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

parser <- parser$parse_args()

bead_counts_file <- parser$bead_counts
allowlist_file <- parser$allowlist
cell_expr_file <- parser$cell_expression
bead_expr_file <- parser$bead_expression
barcode_translate_file <- parser$barcode_translate

# read input data frames
allowlist <- read_csv(allowlist_file, col_names = "CellBarcode")
bead_counts <- read_csv(bead_counts_file)
cell_expr <- read_csv(cell_expr_file)
bead_expr <- read_csv(bead_expr_file)
barcode_translate <- read_csv(barcode_translate_file)

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

# write output
write_csv(bead_cell_data,
          gsub(".barcodeTranslate.tsv", ".bead_summary.csv", barcode_translate_file)
         )
