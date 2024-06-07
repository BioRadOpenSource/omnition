#!/usr/bin/env Rscript
packages <- c("data.table", "dplyr", "actuar", "argparse")
invisible(lapply(packages, library, character.only = TRUE))

options(stringsAsFactors = FALSE)
options(warn = -1)

parser <- ArgumentParser()
#nolint start
parser$add_argument("-o", "--output_directory", type = "character", default = getwd(), help = "Directory to write output files.")
parser$add_argument("-p", "--prefix", type = "character", default = "alignments", help = "Sample id prefix.")
parser$add_argument("-bs", "--bead_summary", type = "character", default = NULL, help = "Path to bead summary csv.")
#nolint end
parser <- parser$parse_args()

out_dir <- parser$output_directory
prefix <- parser$prefix
bead_lambda <- 1.7
bead_summary <- parser$bead_summary

# Set data.table threading
setDTthreads(0)

# #####################################################################################
# FUNCTIONS
# #####################################################################################
prep_input <- function(bead_summary) {
    # ===================================
    # read in bead summary csv with columns of interst
    # If NA cell barcodes, fix with bead barcode, print warning
    # Fixed cell parcode made unique
    # Fix n_beads with NA, this should be 1. Same beads with no cell barcode.
    # Get above knee only
    # Output list with above knee and all bead
    # ===================================
    ### All beads ###
    # Read barcode summary file and select revelant columns
    all_beads_comps <- fread(bead_summary,
        select = c(
            "bead_barcode", "cell_barcode", "n_beads",
            "cell_above_knee"
        )
    )

    # Fix bead bardcodes that dont have cell barcode, print warning
    if (anyNA(all_beads_comps[, cell_barcode])) {
        warn_msg <- paste0(
            "Missing cell barcodes detected.",
            " Using bead barcode in place and setting n_beads = 1"
        )
        warning(warn_msg)
        all_beads_comps[
            is.na(cell_barcode) == TRUE,
            cell_barcode := bead_barcode
        ]
        all_beads_comps[is.na(n_beads) == TRUE, n_beads := 1]
    }

    ### Above knee beads ###
    abv_knee_beads_comps <- all_beads_comps %>%
        .[cell_above_knee == TRUE] %>%
        .[, cell_above_knee := NULL]

    all_beads_comps <- all_beads_comps[, cell_above_knee := NULL]

    list(all_beads = all_beads_comps, abv_knee_beads = abv_knee_beads_comps) %>%
        return(.)
}

calc_poiss <- function(comps, events, lambda) {
    # ===================================
    # Take in data.table from prep_input
    # Correct NA Deconvolution oligos
    # Reformat events per cell barcode
    # Events either DOs or n_beads in drop
    # remove duplicate entries, as cell level stats are repeated
    # If NA detected, post warning
    # Outputs poiss_dt, contain observed and expected
    # ===================================
    # Generate zero-truncated Poisson distribution at indicated lambda
    if (anyNA(comps[, ..events])) {
        warn_msg <- paste0(
            "NAs detects in events. Verify input is correct."
        )
        warning(warn_msg)
    }

    cols2pull <- c("cell_barcode", events)
    tbl <- comps[, ..cols2pull] %>%
        distinct(.)
    per_drop <- as.data.table(table(tbl[, ..events]))
    colnames(per_drop) <- c("drop_size", "Frequency")
    per_drop$drop_size <- as.numeric(as.character(per_drop$drop_size))

    # Poisson distribution
    obs <- rep(0, max(per_drop$drop_size))
    obs[per_drop$drop_size] <- per_drop$Frequency
    # Same poisson function, but zero-truncated
    # This is for n_beads, which should always be a minimum of 1
    expected_poiss <- round(
        sum(per_drop$Frequency) *
            dztpois(c(1:max(per_drop$drop_size)), lambda = lambda),
        5
    )
    counts <- c(1:max(per_drop$drop_size))

    pois_dt <-
        data.table(
            Count = counts,
            Observed = obs,
            Expected = expected_poiss
        )
    return(pois_dt = pois_dt)
}

fmt_dist_tbl <- function(dist_dt, count_bins) {
    # ===================================
    # This functon converts Observed & Expected counts to percentages &
    # groups all counts >= 8 into single 8+ bin
    # 1: Take in data.table from function calc_poiss()
    # 2: Group all counts >= 8 into single bin (8+)
    # 3: Reformat to always include count bins 1-7 & everything 8+
    # 4: Calculate percentage of total for each bin
    # 5: Output percentages with original naming for report
    # ===================================

    # changing to factor to prep fo
    dist_dt[, Count := as.factor(Count)]

    grouped_bin_observed <- dist_dt[!Count %in% count_bins, Observed] %>% sum(.)
    grouped_bin_expected <- dist_dt[!Count %in% count_bins, Expected] %>% sum(.)

    # Create dt with bins that should always exists
    blank_dt <- data.table(Count = count_bins)

    # This creates a table with counting bins we defined in blank_dt
    merged <- dist_dt[blank_dt, on = c("Count")]

    # This will fill in missing bins with 0
    setnafill(merged, fill = 0, cols = c("Observed", "Expected"))

    merged[Count == "8+", Observed := ..grouped_bin_observed]
    merged[Count == "8+", Expected := ..grouped_bin_expected]

    merged[, pct_observed := ((Observed / sum(Observed)) * 100)]
    merged[, pct_expected := ((Expected / sum(Expected)) * 100)]

    merged[, `:=`(Observed = NULL, Expected = NULL)]

    # Changing percent columns to original naming to fit whats expects by reports
    setnames(merged, c("Count", "Observed", "Expected"))
    return(merged)
}

# ####################################################################################
# Function calls
# ####################################################################################
inp_list <- prep_input(bead_summary)

### Calculate Posson distributions ###
poiss_dist_list <- parallel::mclapply(inp_list,
        FUN = calc_poiss,
        events = "n_beads",
        lambda = bead_lambda
)

# ####################################################################################
# Specifying outputs
# ####################################################################################
out_prefix <- paste0(out_dir, "/", prefix)
out_file1 <- paste0(out_prefix, ".above_knee_beads_per_partition_poisson.csv")
out_file2 <- paste0(out_prefix, ".above_knee_beads_per_partition.csv")
out_file3 <- paste0(out_prefix, ".beads_per_partition_poisson.csv")
out_file4 <- paste0(out_prefix, ".beads_per_partition.csv")

# raw, unformatted output
out_file5 <- paste0(out_prefix, ".above_knee_beads_per_partition_raw_poisson.csv")

# ####################################################################################
# Export
# ####################################################################################

# Defining bins we always want to keep or group together
count_bins <- as.factor(c(1:7, "8+"))

### Above knee beads only ###
# Expected distribution
poiss_dist_list$abv_knee_beads %>%
    fmt_dist_tbl(., count_bins = count_bins) %>%
    fwrite(., file = out_file1)

poiss_dist_list$abv_knee_beads %>%
    fwrite(., file = out_file5)

# Observed distribution
poiss_dist_list$abv_knee_beads[, c("Count", "Observed")] %>%
    setNames(., c("drop_size", "Frequency")) %>%
    fwrite(., file = out_file2)

### All beads ###
# Expected distribution
poiss_dist_list$all_beads %>%
    fmt_dist_tbl(., count_bins = count_bins) %>%
    fwrite(., file = out_file3)

# Observed distribution
poiss_dist_list$all_beads[, c("Count", "Observed")] %>%
    setNames(., c("drop_size", "Frequency")) %>%
    fwrite(., file = out_file4)
