#!/usr/bin/env Rscript

library(argparse)
library(dplyr)
library(Matrix)
library(data.table)


cleanFrame <-
  function(data,
           x_colname,
           y_colname,
           x_cutoff = 0,
           y_cutoff = 0) {
    x <- data[[x_colname]]
    y <- data[[y_colname]]
    negative_infinities <- y %in% -Inf
    positive_infinities <- y %in% Inf
    all_NA <- y %in% NA
    all_NA <- all_NA | (x %in% NA)
    all_infinities <-
      negative_infinities | positive_infinities | all_NA
    data <- data[!all_infinities, ]
    x_lt_cutoff <- x < x_cutoff
    data <- data[!x_lt_cutoff, ]
    y_lt_cutoff <- y < y_cutoff
    data <- data[!y_lt_cutoff, ]
    removed_by_cutoff <-
      sum(x_lt_cutoff, na.rm = TRUE) + sum(y_lt_cutoff, na.rm = TRUE)
    return(data)
  }

generateCumFracTable <- function(raw_data) {
  # Generates a UMI table with fractional counts and cumulative fraction
  # Input is a datatable with CELL_BARCODE and NUM_GENIC_READS columns, ie a ed0 file
  # Cut off everything except the columns of interest
  numi_data <- select(raw_data, CELL_BARCODE, NUM_GENIC_READS)
  # Remove all barcodes with zero UMI
  numi_data <- subset(numi_data, NUM_GENIC_READS > 0)
  # Ensure entries are in descending order by NUM_GENIC_READS
  numi_data <- arrange(numi_data, desc(NUM_GENIC_READS))
  # Add column with fractional counts
  fractional_counts <-
    numi_data$NUM_GENIC_READS / sum(numi_data$NUM_GENIC_READS)
  numi_data <-
    cbind(numi_data, FRACTIONAL_COUNTS = fractional_counts)
  numi_data <-
    cbind(numi_data, CUM_FRAC = cumsum(numi_data$FRACTIONAL_COUNTS))
  return(numi_data)
}

distanceFromDiag <- function(curve_frame, x_limit = 0) {
  # Implements the Kneedle algorithm for concave-down, increasing curves
  # Returns the index in the given frame where the knee occurs
  colnames(curve_frame) <- c("x", "y")
  # If an x limit is given, subset the curve frame
  if (x_limit > 0) {
    if (nrow(curve_frame) > x_limit) {
      curve_frame <- curve_frame[1:x_limit, ]
    }
  }
  # Normalize using Kneedle's algorithm
  normalized_x <-
    (curve_frame[, "x"] - min(curve_frame[, "x"])) /
    (max(curve_frame[, "x"]) - min(curve_frame[, "x"]))
  normalized_y <-
    (curve_frame[, "y"] - min(curve_frame[, "y"])) /
    (max(curve_frame[, "y"]) - min(curve_frame[, "y"]))
  normalized_curve <-
    data.frame(cbind(x = normalized_x, y = normalized_y))
  difference_values <-
    apply(normalized_curve, 1, function(row)
      row["y"] - row["x"])
  difference_frame <-
    data.frame(cbind(x = normalized_x, difference = difference_values))
  knee_index <- which.max(difference_frame[, "difference"])
  return(knee_index)
}

#
# Main
#

parser <-
#nolint start
  ArgumentParser(description =
                   "Calculate the number of cells present in an RNA experiment")
parser$add_argument("--count_matrix", help = "path to the UMI count matrix")
parser$add_argument("--results_folder", help = "folder to write output to")
parser$add_argument("--sample_name", help = "the name of the sample being analyzed")
parser$add_argument(
  "--retain_objects",
  action = "store_true",
  default = FALSE,
  help = "keep objects in memory instead of removing them"
)
parser$add_argument("--override", type = "integer",
                    help = "instead of analyzing the sample, use the given value as the number of cells")
parser$add_argument("--loaded", type = "character",
                    help = "The number of cells loaded in the experiment")
#nolint end
args <- parser$parse_args()
args$loaded <- as.integer(args$loaded)

# Set data.table threading
setDTthreads(0)

count_matrix <- readMM(args$count_matrix)
barcodes <-
  read.table(gsub(".mtx.gz", ".barcodes.tsv", args$count_matrix),
              header = F)
genes <-
  read.table(gsub(".mtx.gz", ".genes.tsv", args$count_matrix),
              header = F)
dimnames(count_matrix) <-
  #store barcodes and genes as dimnames for downstream compatibility
  list(genes$V1, barcodes$V1)
umi_sums <- as.numeric(colSums(count_matrix))
raw_numi_data <-
  data.frame(CELL_BARCODE = colnames(count_matrix),
              NUM_GENIC_READS = umi_sums)
# Generate a cumulative fraction table, which can be used the
# same way one generated from a
# ed0 file would be
numi_data <- generateCumFracTable(raw_numi_data)

#
# Now that we have the raw number of umi per cell data, make the modified
# log10 and cumulative fraction data frames
#

# Make graphing-ready frame from the cumulative fraction
kneedle_cumfrac_frame <-
  #nolint start
  data.frame(1:length(numi_data$CUM_FRAC), numi_data$CUM_FRAC)
#nolint end
# Change the column names from x and y to something descriptive
colnames(kneedle_cumfrac_frame) <-
  c("Cell Index", "Cumulative Fraction")
kneedle_cumfrac_frame <-
  cleanFrame(kneedle_cumfrac_frame, "Cell Index", "Cumulative Fraction")

kneedle_cumfrac_knee <-
  tryCatch(
    distanceFromDiag(kneedle_cumfrac_frame, 50000),
    error = function(e)
      e
  )

new_xlimit <- (kneedle_cumfrac_knee * 12.5)

# Recall knee w/ new x_limit for low input samples
if (new_xlimit < 50000) {
  kneedle_cumfrac_knee <-
    tryCatch(
      distanceFromDiag(kneedle_cumfrac_frame, new_xlimit),
      error = function(e) {
        e
      }
    )
}

# This is the "Catch" of try catch
if (inherits(kneedle_cumfrac_knee, "error") ||
    kneedle_cumfrac_knee == 0 ||
    kneedle_cumfrac_knee == nrow(kneedle_cumfrac_frame)) {
  kneedle_cumfrac_knee <- NA
}

if (!args$retain_objects) {
  if (exists("e.out")) {
    rm(e.out)
  }
  if (exists("count_matrix")) {
    rm(count_matrix)
  }
}

# Make a reverse-order x column for the numi data, then subtract
# returned knee from length of x vector
# This removes postive and negative infinities, as well as "NA" from the lists
log10numi_frame <-
  data.frame(
    #nolint start
    1:length(numi_data$NUM_GENIC_READS),
    numi_data$NUM_GENIC_READS,
    log10(1:length(numi_data$NUM_GENIC_READS)),
    #nolint end
    log10(numi_data$NUM_GENIC_READS)
  )
# Change the column names from x and y to something descriptive
colnames(log10numi_frame) <-
  c("Cell Index", "nUMI", "Log10 Cell Index", "Log10 nUMI")
log10numi_frame <-
  cleanFrame(log10numi_frame, "Log10 Cell Index", "Log10 nUMI")

# Garbage collection
gc(verbose = FALSE)

# Set the safety threshold for the number of cells
kneedle_safety <- 250
total_barcodes <- nrow(numi_data)

# Message file to record which knee condition was used
message_file <- paste(args$sample_name, "CELL_CALLING_messages.txt", sep = "_")
message_prefix <- paste0("Warning: [3' RNA Droplet] Sample ", args$sample_name)
message_cells_called <- paste0(
  " had ",
  kneedle_cumfrac_knee,
  " cells called. "
)
message_suffix <- paste0(
  "You can rerun the sample with an override for the knee call,",
  " but it may lead to failures in software."
)

# Choose knee or override
# user_override is set only when kneedle_cumfrac_knee was not used
# Check for the override param value
override <- is.null(args$override)
override_value <- args$override

if (!override && (override_value <= total_barcodes)) { # nolint: cyclocomp_linter.
  writeLines(paste0(
    message_prefix,
    message_cells_called,
    "Overridding called cells with user defined override value of ",
    override_value
  ), message_file)
  # Create the allowlist
  user_override <- override_value
  allowlist <- as.character(numi_data$CELL_BARCODE[1:user_override])
} else if (!override && (override_value > total_barcodes)) {
  writeLines(paste0(
    message_prefix,
    message_cells_called,
    "User defined override of ",
    override_value,
    " is greater than total cells. Using all cells instead."
  ), message_file)
  # If the override is larger than the number of CBC,
  # set the override to the number of CBC
  user_override <- total_barcodes
  allowlist <- as.character(numi_data$CELL_BARCODE[1:user_override])
} else if (is.na(kneedle_cumfrac_knee) && total_barcodes < 10000) {
  # If there's no override, check the kneedle_cumfrac_knee
  # If knee call failed (NA) and there is less that 10k total barcodes, use all
  writeLines(
    paste0(
      message_prefix,
      " cell calling failed with less than 10k total unfiltered cells, ",
      "using all available cells instead. ",
      message_suffix
    ),
    message_file
  )
  allowlist <- as.character(numi_data$CELL_BARCODE[1:total_barcodes])
  user_override <- total_barcodes
} else if (is.na(kneedle_cumfrac_knee) && total_barcodes >= 10000) {
  # If knee call failed (NA) and there is more than 10k barcode, use 10k
  writeLines(paste0(
    message_prefix,
    " cell calling failed, using 10k cells instead. ",
    message_suffix
  ), message_file)
  allowlist <- as.character(numi_data$CELL_BARCODE[1:10000])
  user_override <- 10000
} else if ((kneedle_cumfrac_knee < kneedle_safety) &&
  (total_barcodes >= kneedle_safety)) {
  # If knee call is < safety & total barcodes greater than saftey, use safety
  writeLines(
    paste0(
      message_prefix,
      message_cells_called,
      "Cell call is below safety value, using 250 cells instead. ",
      message_suffix
    ),
    message_file
  )
  allowlist <- as.character(numi_data$CELL_BARCODE[1:kneedle_safety])
  allowlist <- na.omit(allowlist)
  user_override <- kneedle_safety
} else if ((kneedle_cumfrac_knee < kneedle_safety) &&
  (total_barcodes < kneedle_safety)) {
  # If knee call is < safety & total barcodes less than saftey, use total barcodes
  writeLines(
    paste0(
      message_prefix,
      message_cells_called,
      "Cell call is below saftey threshold of 250 ",
      "& total unfiltered cells is below 250. ",
      "Using all cells instead. ",
      message_suffix
    ),
    message_file
  )
  allowlist <- as.character(numi_data$CELL_BARCODE[1:total_barcodes])
  allowlist <- na.omit(allowlist)
  user_override <- total_barcodes
} else if (kneedle_cumfrac_knee > 20000) {
  # If knee is > 20k (2 * max num cells supported), use 10k instead
  writeLines(
    paste0(
      message_prefix,
      message_cells_called,
      ". Cell call is above max supported cells (20k), ",
      "using 10k cells instead. ",
      message_suffix
    ),
    message_file
  )
  allowlist <- as.character(numi_data$CELL_BARCODE[1:10000])
  user_override <- 10000
} else {
  # Use kneedle_cumfrac_knee if all above checks pass
  allowlist <- as.character(numi_data$CELL_BARCODE[1:kneedle_cumfrac_knee])
  user_override <- NA
}

write.table(
  sort(allowlist),
  file = file.path(
    args$results_folder,
    paste(args$sample_name, "_barcode_allowlist.csv", sep = "")
  ),
  row.names = FALSE,
  sep = ",",
  col.names = FALSE,
  quote = F
)

results_table <- rbind(
    c("User_override", user_override),
    c("kneedle_cumfrac_knee", kneedle_cumfrac_knee)
)

# This block is for reporting purposes
if (!is.na(user_override))  {
  # If forced cell overide, final called cells set to user overide
  results_table <- rbind(
    c("final", user_override),
    results_table
  )
} else {
  results_table <- rbind(
    c("final", kneedle_cumfrac_knee),
    results_table
  )
}

results_table <-
  cbind(rep(args$sample_name, nrow(results_table)), results_table)
colnames(results_table) <- c("sample", "variable", "value")
results_table <- data.frame(results_table)

write.table(
  results_table,
  file = file.path(
    args$results_folder,
    paste(args$sample_name, "_numcell_analysis.csv", sep = "")
  ),
  row.names = FALSE,
  sep = ","
)

write.table(
  kneedle_cumfrac_frame,
  file = file.path(
    args$results_folder,
    paste(args$sample_name, "_allalgos_cumfrac.csv", sep = "")
  ),
  row.names = FALSE,
  col.names = TRUE,
  sep = ","
)

write.table(
  log10numi_frame,
  file = file.path(
    args$results_folder,
    paste(args$sample_name, "_allalgos_loglog.csv", sep = "")
  ),
  row.names = FALSE,
  col.names = TRUE,
  sep = ","
)

# Remove the unnecessary columns for the report file
barcode_rank_df <- select(log10numi_frame, Barcode = `Cell Index`, nUMI = `nUMI`)
write.table(
  barcode_rank_df,
  file = file.path(
    args$results_folder,
    paste(args$sample_name, "_barcode_rank.csv", sep = "")
  ),
  row.names = FALSE,
  col.names = TRUE,
  sep = ",",
  quote = FALSE
)

sessionInfo()

# Can release genic umi data now
if (!interactive()) {
  rm(numi_data)
  gc(verbose = FALSE)
}
