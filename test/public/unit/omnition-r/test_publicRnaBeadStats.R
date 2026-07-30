# unit tests to run for publicRnaFalseConnections.R

### Test data ###
# bead_counts_file = counts of aligned and unaligned reads by bead barcode
# allowlist_file = list of cell barcodes above the knee cutoff
# cell_expr_file = counts of reads by cell barcode with categorizations
# bead_expr_file = counts of reads by bead barcode with categorizations
# barcode_translate_file = mapping of bead barcodes to cell barcodes
# sample_name = name of sample

#load input test files
bead_counts_file <- "files/test_publicRnaBeadStats_files/DemoRnaMixed_S1.counts_per_barcode.csv"
allowlist_file <- "files/test_publicRnaBeadStats_files/DemoRnaMixed_S1_barcode_allowlist.csv"
cell_expr_file <- "files/test_publicRnaBeadStats_files/DemoRnaMixed_S1.final.scrnaseq_counts.csv"
bead_expr_file <- "files/test_publicRnaBeadStats_files/DemoRnaMixed_S1.bead_scrnaseq_counts.csv"
barcode_translate_file <- "files/test_publicRnaBeadStats_files/DemoRnaMixed_S1_barcodeTranslate.tsv"
sample_name <- "DemoRnaMixed_S1"

# load comparison output files
bead_summary_file <- "files/test_publicRnaBeadStats_files/DemoRnaMixed_S1.bead_summary.csv"
metrics_file <- "files/test_publicRnaBeadStats_files/DemoRnaMixed_S1.bead_merge_sum_exp_metrics.csv"

# Read lines and include only dependency and functions
lines <- scan("../../../../bin/publicRnaBeadStats.R", what = character(), sep = "\n")

# find start of function calls
function_call_start <- grep("# read input data frames", lines)

# trim down sourced file excluding function calls
lines <- lines[-c(function_call_start:length(lines))]

# source the trimmed file
source(textConnection(paste(lines, collapse = "\n")))

### Test assemble_bead_to_cell_table ###
# ===================================
# Merges bead counts, cell expression, bead expression, and barcode translate data
# Adds a column to designate which cells are above the knee cutoff
# ===================================

# Read in test data files
allowlist <- read_csv(allowlist_file, col_names = "CellBarcode")
bead_counts <- read_csv(bead_counts_file)
cell_expr <- read_csv(cell_expr_file)
bead_expr <- read_csv(bead_expr_file)
barcode_translate <- read_csv(barcode_translate_file)

# Set test case variables
row_count <- 952
column_count <- 42
above_knee_count <- 805

# Read in comparison output file
bead_summary <- read_csv(bead_summary_file)

# Get function output
bead_cell_data <- assemble_bead_to_cell_table(
    bead_counts,
    allowlist,
    cell_expr,
    bead_expr,
    barcode_translate
)

# Test
test_that("Pass - assemble_bead_to_cell_table", {
    # Test that input is loaded in correctly
    bead_cell_data %>%
        is.tibble(.) %>%
        expect_true(.)
    expect_equal(nrow(bead_cell_data), row_count)
    expect_equal(ncol(bead_cell_data), column_count)

    # Check cells
    bead_cell_data %>%
        filter(cell_above_knee == TRUE) %>%
        nrow(.) %>%
        expect_equal(., above_knee_count)

    # Test that the df matches the bead summary file
    expect_equal(
        bead_cell_data,
        # group by cell_barcode
        bead_summary %>% group_by(cell_barcode)
    )
})

### Test calc_bead_merging_stats ###
# ===================================
# Calculates the number of cells, number of cells with merges,
# and the percentage of cells with multiple beads
# ===================================

# Set test case variables
n_cells <- 800
cells_with_merges <- 31
percent_cells_with_multiple_beads <- round(cells_with_merges / n_cells * 100, 2)

# Read in comparison output file
all_metrics <- read_csv(metrics_file)

# Get function output
merging_metrics <- calc_bead_merging_stats(
    bead_cell_data,
    sample_name
)

# Test
test_that("Pass - calc_bead_merging_stats", {
    # Test that input is loaded in correctly
    merging_metrics %>%
        is.data.frame(.) %>%
        expect_true(.)

    # Test that the metrics are correct
    expect_equal(
        (merging_metrics %>% filter(metric == "n_cells") %>% select(value))[[1,1]],
        n_cells
    )
    expect_equal(
        (merging_metrics %>% filter(metric == "cells_with_merges") %>% select(value))[[1,1]],
        cells_with_merges
    )
    expect_equal(
        (merging_metrics %>% filter(metric == "percent_cells_with_multiple_beads") %>% select(value))[[1,1]],
        percent_cells_with_multiple_beads
    )

    # Test that the df matches the metric summary file
    expect_equal(
        merging_metrics,
        as.data.frame(all_metrics %>% filter(process == "bead_merging"))
    )
})

### Test calc_sum_expr_metrics ###
# ===================================
# Sums columns in cell_expr
# Calculates percentages
# Calculates fraction of reads in cells for single and mixed species data
# ===================================

# Set test case variables
total_fraction <- 94.1
human_fraction <- 92.2
mouse_fraction <- 96.0

# Get function output
summarized_metrics <- calc_sum_expr_metrics(
    cell_expr,
    allowlist,
    sample_name
)

# Test
test_that("Pass - calc_sum_expr_metrics", {
    # Test that input is loaded in correctly
    summarized_metrics %>%
        is.tibble(.) %>%
        expect_true(.)

    # Test that the metrics are correct
    expect_equal(
        round((summarized_metrics %>% filter(
            sample == sample_name,
            metric == "fraction_reads_in_cells"
        ) %>% select(value))[[1,1]], 1),
        total_fraction
    )
    expect_equal(
        round((summarized_metrics %>% filter(
            sample == paste0(sample_name, "_DemoReferenceRnaMixedHomoSapiens_S1"),
            metric == "fraction_reads_in_cells"
        ) %>% select(value))[[1,1]], 1),
        human_fraction
    )
    expect_equal(
        round((summarized_metrics %>% filter(
            sample == paste0(sample_name, "_DemoReferenceRnaMixedMusMusculus_S2"),
            metric == "fraction_reads_in_cells"
        ) %>% select(value))[[1,1]], 1),
        mouse_fraction
    )

    # Test that the df matches the metric summary file
    expect_equal(
        summarized_metrics,
        (all_metrics %>% filter(process == "summarize_expression"))[]
    )
})