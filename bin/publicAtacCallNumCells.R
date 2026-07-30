#!/usr/bin/env Rscript

library(ggplot2)
library(argparse)
library(dplyr)
library(DropletUtils)
library(tools)
library(Matrix)
library(BUSpaRse)
library(Seurat)
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

secondDerivMin <- function(curve_frame, spar_val, tol = 0) {
  # Use a fitting spline and the 'minimum of second derivative'
  # algorithm to call a knee point
  # The function behaves differently if you give it any tol
  # value, can't find a suitable default
  if (tol > 0) {
    spl <-
      smooth.spline(curve_frame[, 1],
                    y = curve_frame[, 2],
                    spar = spar_val,
                    tol = tol)
  } else{
    spl <-
      smooth.spline(curve_frame[, 1], y = curve_frame[, 2], spar = spar_val)
  }
  pred <- predict(spl)
  pred.prime <- predict(spl, deriv = 1)
  pred.second_deriv <- predict(spl, deriv = 2)
  # The indices may not be the same in the orig data and in pred
  # (if a value for tol is given), but the values are
  # So we have to take the values at the knee point, then
  # figure out at which indices in the original
  # data these values lie
  spline_knee_x <- pred$x[which.min(pred.second_deriv$y)]
  # Find the closest x value in the original data to the knee x value
  knee_index <- which.min(abs(curve_frame[, 1] - spline_knee_x))
  knee_point <- curve_frame[, 1][knee_index]
  knee_y_value <- curve_frame[, 2][knee_index]
  pred_frame <- data.frame(pred)
  orig_plot <- ggplot() +
    geom_point(
      data = curve_frame,
      aes(
        x = colnames(curve_frame)[1],
        y = colnames(curve_frame)[2],
        color = "original"
      ),
      size = 1.5
    ) +
    geom_line(data = pred_frame, aes(x = x, y = y, color = "fit_line")) +
    scale_color_manual(
      "",
      breaks = c("original", "fit_line", "knee_point"),
      values = c(
        "original" = "black",
        "fit_line" = "green",
        "knee_point" = "red"
      )
    ) +
    expand_limits(x = 0, y = 0) +
    geom_vline(aes(xintercept = knee_point, color = "knee_point")) +
    geom_hline(aes(yintercept = knee_y_value, color = "knee_point")) +
    labs(title = "Original Data with Fit Line") +
    xlab("Barcode Rank") +
    ylab("Num UMI")
  pred_prime_frame <- data.frame(pred.prime)
  deriv_plot <-
    ggplot(data = pred_prime_frame, aes(x = x, y = y, color = "deriv")) +
    geom_line() +
    geom_vline(aes(xintercept = knee_point, color = "knee_point")) +
    scale_color_manual(
      "",
      breaks = c("deriv", "knee_point"),
      values = c("deriv" = "black", "knee_point" = "red")
    ) +
    labs(title = "First Derivative of Fit Line") +
    xlab("Barcode Rank") +
    ylab("First Derivative Value")
  pred_doubleprime_frame <- data.frame(pred.second_deriv)
  second_deriv_plot <-
    ggplot(data = pred_doubleprime_frame, aes(x = x, y = y, color = "second_deriv")) +
    geom_line() +
    geom_vline(aes(xintercept = knee_point, color = "knee_point")) +
    scale_color_manual(
      "",
      breaks = c("second_deriv", "knee_point"),
      values = c("second_deriv" = "black", "knee_point" = "red")
    ) +
    labs(title = "Second Derivative of Fit Line") +
    xlab("Barcode Rank") +
    ylab("Second Derivative Value")
  return(knee_index)
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
  normalized_plot <-
    ggplot(data = normalized_curve, aes(x = x, y = y, color = "normalized")) +
    geom_line() +
    geom_line(data = difference_frame, aes(x = x, y = difference, color =
                                             "difference")) +
    geom_vline(aes(xintercept = normalized_curve[, "x"][knee_index])) +
    scale_color_manual(
      "",
      breaks = c("normalized", "difference"),
      values = c("normalized" = "black", "difference" = "red")
    ) +
    labs(title = "Normalized Curve")
  orig_plot <- ggplot(data = curve_frame, aes(x = x, y = y)) +
    geom_line() +
    geom_vline(aes(xintercept = curve_frame[, "x"][knee_index])) +
    labs(title = paste("Original Curve ", args$sample_name, sep = ""))
  return(knee_index)
}

#
# Main
#

# The number of UMI's below which dropletuitls barcodrRanks
# and emptyDrops assume cells are empty droplets
dropletutils_lower <- 50

parser <-
#nolint start
  ArgumentParser(description =
                   "Calculate the number of cells present in an experiment")
parser$add_argument("--ed0_file",
                    help = "path to .ed0min0readsAll.synthesis_stats.txt file for given sample")
parser$add_argument("--count_matrix", help = "path to the UMI count matrix")
parser$add_argument("--results_folder", help = "folder to write output to")
parser$add_argument("--sample_name", help = "the name of the sample being analyzed")
parser$add_argument("--counter_output",
                    help = "output folder for the sample from kallisto/bustools")
parser$add_argument(
  "--retain_objects",
  action = "store_true",
  default = FALSE,
  help = "keep objects in memory instead of removing them"
)
parser$add_argument("--override", type = "integer",
                    help = "instead of analyzing the sample, use the given value as the number of cells")
parser$add_argument("--loaded", type = "integer",
                    help = "The number of cells loaded in the experiment")
#nolint end
args <- parser$parse_args()

# Set data.table threading
setDTthreads()

# Read in data
if (!is.null(args$count_matrix)) {
  count_matrix <-
    fread(file = args$count_matrix,
          sep = ",",
          header = FALSE)
  raw_numi_data <-
    data.frame(CELL_BARCODE = count_matrix$V1,
               NUM_GENIC_READS = count_matrix$V2)
  # Generate a cumulative fraction table
  numi_data <- generateCumFracTable(raw_numi_data)
}

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
reverse_x <-
  abs(log10numi_frame[, "Log10 Cell Index"] -
        max(log10numi_frame[, "Log10 Cell Index"]))

rev_log10numi_frame <-
  data.frame(reverse_x, log10numi_frame[, "Log10 nUMI"])
colnames(rev_log10numi_frame) <-
  c("Reverse Cell Index", "Log10 nUMI")
rev_log10numi_frame <-
  cleanFrame(rev_log10numi_frame, "Reverse Cell Index", "Log10 nUMI")

second_deriv_min_frame <-
  data.frame(
    #nolint start
    1:length(numi_data$NUM_GENIC_READS),
    numi_data$NUM_GENIC_READS,
    log10(1:length(numi_data$NUM_GENIC_READS)),
    #nolint end
    log10(numi_data$NUM_GENIC_READS)
  )
colnames(second_deriv_min_frame) <-
  c("Cell Index", "nUMI", "Log10 Cell Index", "Log10 nUMI")
second_deriv_min_frame <-
  cleanFrame(second_deriv_min_frame, "Log10 Cell Index", "Log10 nUMI", 0, 0.5)

# Do emptydrops first so that we can free the memory from count_matrix sooner
skip_emptydrops <- TRUE
if ((!is.null(args$count_matrix)) ||
    (!is.null(args$counter_output))) {
  skip_emptydrops <- FALSE
}


if (!skip_emptydrops) {
  # Emptydrops only works on whole count matrices, not num genic reads per barcode
  e.out <- tryCatch(
    emptyDrops(
      count_matrix,
      lower = dropletutils_lower,
      test.ambient = TRUE,
      retain = 500
    ),
    error = function(e)
      e
  )
  if (inherits(e.out, "error") || !(any(e.out$FDR <= 0.01))) {
    emptydrops_num_cells <- NA
  } else {
    # Sort even though it's probably not necessary
    e.out <- e.out[order(e.out$Total, decreasing = T), ]
    # Select everything with false discovery rate score less than 1%
    is.cell <- e.out$FDR <= 0.01
    # Replace any "NA" with FALSE so we can use this as a logical vector for subsetting
    is.cell[is.na(is.cell)] <- FALSE
    # Subset rows of emptydrops output by whether or not it is called a cell
    called_cells <- e.out[is.cell, ]
    # Write a file with all the called cell cbcs
    write.table(
      called_cells,
      file = file.path(
        args$results_folder,
        paste(args$sample_name, "_called_cells.csv", sep = "")
      ),
      row.names = TRUE,
      sep = ","
    )
    # Count the number of called cells
    emptydrops_num_cells <- sum(is.cell, na.rm = TRUE)
    # Make a tibble of FDR vs FDR barcode rank
    FDR_tibble <-
      #nolint start
      tibble(rank = 1:length(sort(e.out$FDR)),
             #nolint end
             FDR = sort(e.out$FDR))
    # Only graph the first 10k, it automatically truncates if there's less than 10k CBCs
    emptydrops_FDR_plot <-
      ggplot(data = FDR_tibble[1:10000, ], aes(rank, FDR)) +
      geom_point() +
      labs(
        title = paste(
          args$sample_name,
          " EmptyDrops False Discovery Ratio, first 10K CBCs"
        ),
        sep = ""
      )
    if (!args$retain_objects) {
      rm(is.cell)
    }
  }
}
if (!args$retain_objects) {
  rm(e.out)
  rm(count_matrix)
}
gc(verbose = FALSE)


kneedle_cumfrac_knee <-
  tryCatch(
    distanceFromDiag(kneedle_cumfrac_frame, 50000),
    error = function(e)
      e
  )
# This is the "Catch" of try catch
if (inherits(kneedle_cumfrac_knee, "error") ||
    kneedle_cumfrac_knee == 0 ||
    kneedle_cumfrac_knee == nrow(kneedle_cumfrac_frame)) {
  kneedle_cumfrac_knee <- NA
}

# Use the kneedle algorithm to find the knee of the log-log
# numi graph, selects the "top" of the knee
rev_log10numi_index <-
  tryCatch(
    distanceFromDiag(rev_log10numi_frame),
    error = function(e)
      e
  )
rev_log10numi_frame_cbcs <- nrow(rev_log10numi_frame)
# Remove the data frame and garbage collect the memory
# if not running in interactive mode
if (!interactive()) {
  rm(rev_log10numi_frame)
  gc(verbose = FALSE)
}
if (inherits(rev_log10numi_index, "error") ||
    rev_log10numi_index == 0 ||
    rev_log10numi_index == rev_log10numi_frame_cbcs) {
  rev_log10numi_index <- NA
}


# Estimate a knee point by finding the minimum of the second derivative
second_deriv_min_knee <-
  tryCatch(
    secondDerivMin(second_deriv_min_frame[c("Log10 Cell Index", "Log10 nUMI")], 0.7, 0),
    error = function(e)
      e
  )
second_deriv_min_frame_cbcs <- nrow(second_deriv_min_frame)
# Remove the data frame and garbage collect the memory
# if not running in interactive mode
if (!interactive()) {
  rm(second_deriv_min_frame)
  gc(verbose = FALSE)
}
if (inherits(second_deriv_min_knee, "error") ||
    second_deriv_min_knee == 0 ||
    second_deriv_min_knee == second_deriv_min_frame_cbcs) {
  second_deriv_min_knee <- NA
}


# This tricks DropletUtils into thinking you've given it a count
#  matrix with only one gene
# It returns the UMI count at the knee, not a number of cells
# This can fail unexpectedly, likely due to barcodeRanks smooth.
# spline default parameters
DU_ranks <- tryCatch(
  barcodeRanks(t(numi_data[, 2]), lower = dropletutils_lower),
  error = function(e)
    e
)
if (inherits(DU_ranks, "error")) {
  DU_knee_index <- NA
  DU_inflection_index <- NA
} else{
  DU_knee_yval <- metadata(DU_ranks)$knee
  DU_inflection_yval <- metadata(DU_ranks)$inflection
  DU_knee_index <-
    which.min(abs(numi_data$NUM_GENIC_READS - DU_knee_yval))
  DU_inflection_index <-
    which.min(abs(numi_data$NUM_GENIC_READS - DU_inflection_yval))
}

# Create a datatable to store the results
prelim_results_table <- rbind(
  c("DropletUtils_inflection", DU_inflection_index),
  c("kneedle_cumfrac_knee", kneedle_cumfrac_knee),
  c("kneedle_log10_knee", rev_log10numi_index),
  c("second_deriv_min_knee", second_deriv_min_knee),
  c("DropletUtils_knee", DU_knee_index)
)

if (!skip_emptydrops) {
  prelim_results_table <-
    rbind(prelim_results_table,
          c("emptydrops_num_cells", emptydrops_num_cells))
}

# Create the allowlist
# Get the override value to add to the results table
if (!is.null(args$override)) {
  # If the override is larger than the number of CBC,
  # set the override to the number of CBC
  if (args$override <= nrow(numi_data)) {
    user_override <- args$override
  } else{
    user_override <- nrow(numi_data)
  }
  allowlist <- as.character(numi_data$CELL_BARCODE[1:user_override])
# If there's no override, check the knee call against the safety
} else if (!all(is.na(prelim_results_table[, 2]))) {
  best_knee <-
    as.numeric(prelim_results_table[, 2][min(which(!is.na(prelim_results_table[, 2])))])
  if (!is.null(args$loaded)) {
    # The knee can't be > 2x the loaded value
    if (best_knee > args$loaded * 2) {
      best_knee <- args$loaded * 2
      user_override <- args$loaded * 2
    } else {
      user_override <- NA
    }
  } else {
    user_override <- NA
  }
  allowlist <- as.character(numi_data$CELL_BARCODE[1:best_knee])
# If all algorithms fail, set to the smaller of all barcodes or 2x loaded
} else {
  if (!is.null(args$loaded)) {
    if (nrow(numi_data) < args$loaded * 2) {
      allowlist <- as.character(numi_data$CELL_BARCODE)
    } else {
      allowlist <- as.character(numi_data$CELL_BARCODE[1:(args$loaded * 2)])
    }
  } else {
    allowlist <- as.character(numi_data$CELL_BARCODE)
  }
}

write.table(
  allowlist,
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
    prelim_results_table)

results_table <-
  cbind(rep(args$sample_name, nrow(results_table)), results_table)
colnames(results_table) <- c("sample", "variable", "value")
results_table <- data.frame(results_table)
# Coerce values to numeric so they can be used by ggplot2
results_table$value <- as.numeric(as.character(results_table$value))

write.table(
  results_table,
  file = file.path(
    args$results_folder,
    paste(args$sample_name, "_numcell_analysis.csv", sep = "")
  ),
  row.names = FALSE,
  sep = ","
)

# Show the calculated knee points on a plot
cumfrac_plot <-
  ggplot(data = kneedle_cumfrac_frame,
         aes(x = `Cell Index`, y = `Cumulative Fraction`)) +
  geom_line() +
  geom_vline(data = results_table, mapping = aes(xintercept = value)) +
  geom_text(
    data = results_table,
    mapping = aes(x = value, y = 0, label = variable),
    angle = 90,
    size = 3,
    vjust = -0.5,
    hjust = 0
  ) +
  labs(title = paste(args$sample_name, " Cumulative Fraction", sep = ""))
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


# Show the calculated knee points on a plot
log10_plot <-
  ggplot(data = log10numi_frame, aes(x = `Log10 Cell Index`, y = `Log10 nUMI`)) +
  geom_line() +
  geom_vline(data = results_table, mapping = aes(xintercept = log10(value))) +
  geom_text(
    data = results_table,
    mapping = aes(
      x = log10(value),
      y = 0,
      label = variable
    ),
    angle = 90,
    size = 3,
    vjust = -0.5,
    hjust = 0
  ) +
  labs(title = paste(args$sample_name, " Log10 nUMI", sep = ""))

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



# Create a plot with just kneedle_cumfrac line on the log-log umi data
kneedle_cumfrac_log10_plot <-
  ggplot(data = log10numi_frame, aes(x = `Log10 Cell Index`, y = `Log10 nUMI`)) +
  geom_line() +
  geom_vline(aes(xintercept = log10(kneedle_cumfrac_knee))) +
  labs(title = paste(args$sample_name, " kneedle_cumfrac_knee on log10 nUMI", sep =
                       ""))

# Create a plot with just kneedle_cumfrac line on the cumfrac data
kneedle_cumfrac_plot <-
  ggplot(kneedle_cumfrac_frame,
         aes(x = `Cell Index`, y = `Cumulative Fraction`)) +
  geom_line() +
  geom_vline(aes(xintercept = kneedle_cumfrac_knee)) +
  labs(title = paste(args$sample_name, " kneedle_cumfrac_knee on cumfrac nUMI", sep =
                       ""))

# Can release genic umi data now
if (!interactive()) {
  rm(numi_data)
  gc(verbose = FALSE)
}
