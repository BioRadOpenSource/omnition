# unit tests to run for rnaMergeBeads.R

### Test data ###
# droplet_2 = No false connection when transcript cutoff applied
# droplet_2 inclusion of beads w/ <200 will result it in 4 false connections
# droplet_3 has 2 false connections
# CCTCCTTAATGGCCATGGTGT = Cell barocde = NA, >200, DO=NA, human
# CGGCCAGCACGCCAATGGTGT = Cell barocde = NA, >200, DO=NA, mixed
# AATTCCATCAGCAACGCAATC = Cell barocde = NA, >200 DO=NA, mouse
# CATCCCATCACCAACGGAATC = Cell barocde = NA, >200, above knee = false, DO=NA, Mixed

# load test file
test_input <- "files/rnaWorkflow.bead_summary.csv"

# Read lines and include only dependency and functions
lines <- scan("../../../bin/rnaCalcBeadPlots.R", what = character(), sep = "\n")

# find start of function calls
function_call_start <- grep("# Function calls", lines)

# trim down sourced file exclusing function calls
lines <- lines[-c(function_call_start:length(lines))]

# source the trimmed file
source(textConnection(paste(lines, collapse = "\n")))

### Test prep_input ###
# ===================================
# read in bead summary csv with columns of interest
# If NA cell barcodes, fix with bead barcode, print warning
# Fixed cell parcode made unique
# Fix n_beads with NA, this should be 1. Same beads with no cell barcode.
# Make separate above knee only data.table
# Output list with above knee and all bead
# ===================================

# Prep test cases
above_knee_count <- 10
all_bead_count <- 13
expected_cols_count <- 3
expected_all_nbeads_counts <- c(1, 1, 1, 1)
expected_abv_knee_nbeads_counts <- c(1, 1, 1)

missing_cellbarcode_all <- c(
    "CCTCCTTAATGGCCATGGTGT",
    "CGGCCAGCACGCCAATGGTGT",
    "AATTCCATCAGCAACGCAATC",
    "CATCCCATCACCAACGGAATC"
)

missing_cellbarcode_abv_knee <- c(
    "CCTCCTTAATGGCCATGGTGT",
    "CGGCCAGCACGCCAATGGTGT",
    "AATTCCATCAGCAACGCAATC"
)

warn_msg <- paste0(
    "Missing cell barcodes detected.",
    " Using bead barcode in place"
)

# Get function output
prep_output <- prep_input(test_input)

# Test
test_that("Pass - prep_input", {
    # Test that input is loaded in correctly
    # Class checks
    prep_output %>%
        is.list(.) %>%
        expect_true(.)

    # Class checks
    prep_output$all_beads %>%
        is.data.table(.) %>%
        expect_true(.)
    prep_output$abv_knee_beads %>%
        is.data.table(.) %>%
        expect_true(.)

    # bead and above knee bead counts
    expect_equal(nrow(prep_output$all_beads), all_bead_count)
    expect_equal(ncol(prep_output$all_beads), expected_cols_count)
    expect_equal(nrow(prep_output$abv_knee_beads), above_knee_count)
    expect_equal(ncol(prep_output$abv_knee_beads), expected_cols_count)

    # Test that NAs detected in cell barcode
    expect_warning(
        prep_input(test_input),
        warn_msg
    )

    # Test that cell barcode corrected with bead barcode
    prep_output$all_beads[cell_barcode == bead_barcode, bead_barcode] %>%
        expect_setequal(., missing_cellbarcode_all)
    prep_output$abv_knee_beads[cell_barcode == bead_barcode, bead_barcode] %>%
        expect_setequal(., missing_cellbarcode_abv_knee)
    # Test that when bead barcode not merged, n_beads = 1
    prep_output$all_beads[cell_barcode == bead_barcode, n_beads] %>%
        expect_setequal(., expected_all_nbeads_counts)
    prep_output$abv_knee_beads[cell_barcode == bead_barcode, n_beads] %>%
        expect_setequal(., expected_abv_knee_nbeads_counts)
})

### Test calc_poiss ###
# ===================================
# Take in data.table from prep_input
# Correct NA deconvolution oligos
# Reformat events per cell barcode
# Events either DOs or n_beads in drop
# remove duplicate entries, as cell level stats are repeated
# If NA detected, post warning
# Output expected distribution (pois_df)
# Output observed (per drop)
# Outputs is list with 2x data.tables inside
# ===================================
# Prep test cases
expected_abv_knee_total <- 11
expected_all_total <- 12
warn_msg2 <- paste0(
    "NAs detects in events. Verify input is correct"
)

poiss_beads_ouput <- calc_poiss(
    comps = prep_output$all_beads,
    events = "n_beads",
    lambda = 3
)

poiss_abv_knee_beads_ouput <- calc_poiss(
    comps = prep_output$abv_knee_beads,
    events = "n_beads",
    lambda = 3
)


test_that("Pass - calc_poiss", {
        # Class checks
        poiss_beads_ouput %>%
            is.data.table(.) %>%
            expect_true(.)
        poiss_abv_knee_beads_ouput %>%
            is.data.table(.) %>%
            expect_true(.)

        # Check above knee beads add up
        poiss_abv_knee_beads_ouput[Observed > 0, totals := (Count * Observed)] %>%
            .[is.na(totals), totals := 0] %>%
            .[, totals] %>%
            sum(.) %>%
            expect_equal(., expected_abv_knee_total)

        # Check all beads add up
        poiss_beads_ouput[Observed > 0, totals := (Count * Observed)] %>%
            .[is.na(totals), totals := 0] %>%
            .[, totals] %>%
            sum(.) %>%
            expect_equal(., expected_all_total)
})

### Test fmt_dist_tbl ###
# ===================================
# This functon converts Observed & Expected counts to percentages &
# groups all counts >= 8 into single 8+ bin
# 1: Take in data.table from function calc_poiss()
# 2: Group all counts >= 8 into single bin (8+)
# 3: Reformat to always include count bins 1-7 & everything 8+
# 4: Calculate percentage of total for each bin
# 5: Output percentages with original naming for report
# ===================================

# Prep test cases
# This is for defining the different ranges of bins that will be tested
bin_ranges <- list(
    small = c(1:5),
    exact = c(1:8),
    big = c(1:21),
    gap = c(1:3, 6:12)
)

# This table is for testing when we have exact number of desired bins
exact_test <- data.table(
    Count = as.factor(bin_ranges$exact),
    Observed = bin_ranges$exact,
    Expected = bin_ranges$exact
)

# This table is for testing when  there is less than desired number of bins
small_test <- data.table(
    Count = as.factor(bin_ranges$small),
    Observed = bin_ranges$small,
    Expected = bin_ranges$small
)

# This table is for testing when there is more than desired number of bins
big_test <- data.table(
    Count = as.factor(bin_ranges$big),
    Observed = bin_ranges$big,
    Expected = bin_ranges$big
)

# This table is for testing when there is missing bins inbetween
gap_test <- data.table(
    Count = as.factor(bin_ranges$gap),
    Observed = bin_ranges$gap,
    Expected = bin_ranges$gap
)

test_tbls <- list(
    exact = exact_test,
    small = small_test,
    big = big_test,
    gap = gap_test
)

# Define expectations
# We always only want 8 bins, with the 8th including all counts >=8
fmt_dist_tbl_out_nrow <- 8
fmt_dist_tbl_out_ncol <- 3
fmt_dist_tbl_out_colnames <- c("Count", "Observed", "Expected")
fmt_dist_tbl_out_fct_lvls <- as.factor(c(1:7, "8+"))

expected_pct <- list(
    small = c(0.0667, 0.1333, 0.2000, 0.2667, 0.3333, 0, 0, 0) * 100,
    exact = c(0.0278, 0.0556, 0.0833, 0.1111, 0.1389, 0.1667, 0.1944, 0.2222) * 100,
    big = c(0.0043, 0.0087, 0.0130, 0.0173, 0.0216, 0.0260, 0.0303, 0.8788) * 100,
    gap = c(0.0145, 0.0290, 0.0435, 0, 0, 0.0870, 0.1014, 0.7246) * 100
)

# Set up and call function
count_bins <- as.factor(c(1:7, "8+"))

fmt_dist_tbl_out <- lapply(test_tbls, FUN = fmt_dist_tbl, count_bins = count_bins)

# Func for running multiple checks on outputs for each of 4 test objects
check_fmt_dist_tbl <- function(x) {
    DT <- fmt_dist_tbl_out[x][[1]]
    expected_pct_values <- expected_pct[x][[1]] %>% round(., 2)

    # The all() returns a single true is all values in vec are true
    # This makes downstream test eval easier
    results <- c(
        class_check = is.data.table(DT),
        nrow_check = nrow(DT) == fmt_dist_tbl_out_nrow,
        colname_check = (colnames(DT) == fmt_dist_tbl_out_colnames) %>% all(.),
        ncol_check = ncol(DT) == fmt_dist_tbl_out_ncol,
        factor_lvl_check = (DT$Count == fmt_dist_tbl_out_fct_lvls) %>% all(.),
        observed_check = (round(DT$Observed, 2) == expected_pct_values) %>% all(.),
        observed_sum_check = round(sum(DT$Observed)) == 100,
        expected_sum_check = round(sum(DT$Expected)) == 100
    )

    return(results)
}

fmt_dist_tbl_test_results <- lapply(
    names(fmt_dist_tbl_out),
    FUN = check_fmt_dist_tbl) %>%
    setNames(., names(fmt_dist_tbl_out))

# Check that each object's output passed all tests
fmt_dist_tbl_test_summary <- fmt_dist_tbl_test_results %>%
    sapply(., FUN = all)

test_that("Pass - fmt_dist_tbl", {
    expect_true(fmt_dist_tbl_test_summary["small"],
        info = "small test object failed, see fmt_dist_tbl_test_results['small']"
    )

    expect_true(fmt_dist_tbl_test_summary["exact"],
        info = "exact test object failed, see fmt_dist_tbl_test_results['exact']"
    )

    expect_true(fmt_dist_tbl_test_summary["big"],
        info = "big test object failed, see fmt_dist_tbl_test_results['big']"
    )

    expect_true(fmt_dist_tbl_test_summary["gap"],
        info = "gap test object failed, see fmt_dist_tbl_test_results['gap']"
    )
})
