#!/usr/bin/env Rscript
# Caleb Lareau
# Bio-Rad Laboratories, Inc.

options(warn = -1)

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tools)))

"%ni%" <- Negate("%in%")

options(scipen = 999)

# This function / script for processed chromosome
# files to produce consensus barcode doublets
# as well as summary and QC metrics and visuals

substrRight <- function(x, n) {
  substr(x, nchar(x) - n + 1, nchar(x))
}

# Input Output
args <- commandArgs(trailingOnly = FALSE)

# Grab knee call functions from figuring out where this script lives
#  and then changing the file name
needle <- "--file="
match <- grep(needle, args)

# Simple function to determine the mode of a vector
# See here: https://stackoverflow.com/questions/
# 2547402/is-there-a-built-in-function-for-finding-the-mode
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


get_local_minima_CL <- function(x) {
  # Find the indicies where we go from negative diff's to positive
  # if this returns less then 0 then we are descending
  y <- diff(c(Inf, x)) < 0L
  # rle to get the number of trues (downwards steps) and falses (upwards)
  # cumsum will then the index at which each of these inflections happens
  y <- cumsum(rle(y)$lengths)
  # If we only get three values, then this becomes a problem--
  # the every/other TRUE,FALSE,TRUE will only keep the extreme values
  # and then remove the only true minimum
  if (length(y) > 3) {
    y <- y[seq.int(1L, length(y), 2L)]
  }
  # Seth's modification for removing duplicated elements at the beginning
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

# Complicated function converted by CL
get_density_threshold_CL <-
  function(count_vector = NULL,
           type,
           logTransform = TRUE) {
    # Initialize using some reasonable value and filter anything below
    threshold <- min(get_mode(count_vector)) * 0.001
    filtered_tbl <-
      data.frame(count = count_vector[count_vector > threshold])
    # Parameterize the log transformation to work with non-count data
    # May get confusing evtually but for now, seemingly a decent hack
    if (logTransform) {
      filtered_tbl$log_counts <- log10(filtered_tbl$count)
    } else{
      filtered_tbl$log_counts <- filtered_tbl$count
    }
    # Calculate the density using a gaussian kernel
    xx_values <- 10000
    vector_density <-
      density(
        filtered_tbl$log_counts,
        bw = 0.1,
        kernel = "gaussian",
        n = xx_values,
        from = min(filtered_tbl$log_counts),
        to = max(filtered_tbl$log_counts)
      )
    local_mins <- get_local_minima_CL(vector_density$y)
    # If a minima was called at the very start or end of the distribution, remove it
    local_mins <-
      local_mins[!(local_mins == 1 |
                     local_mins == length(vector_density$y))]
    # Find the first local min that meets the following criteria, starting right to left
    local_min <- Find(function(x) {
      # Make sure that the selected min includes at least 20% of barcodes
      # and that the difference between the min and the max differences
      # by some appreciable amount
      # both in terms of absolute difference (Caleb changed 0.5
      # to 0.05) AND relative difference
      abs_difference <- ifelse(logTransform, 0.5, 0.05)
      return(x >= (0.2 * xx_values) &
                ((max(filtered_tbl$log_counts) - vector_density$x[x]) > abs_difference |
                   (vector_density$x[x] < max(filtered_tbl$log_counts) / 2)
                ))
    }, rev(local_mins))
    if (!is.null(local_min)) {
      if (logTransform) {
        threshold <- 10 ^ vector_density$x[local_min]
      } else {
        threshold <- vector_density$x[local_min]
      }
      message("Setting knee threshold to: ", threshold)
    } else {
      # Null local min-- take a best guess
      if (logTransform) {
        message("No reliable knee found-- setting threshold to 0")
        threshold <- 0
      } else {
        message("No reliable knee found-- setting threshold to 0")
        threshold <- 0
      }
      local_min <- 1
      local_mins <- 1
    }
    safety <- 0
    # Safe guard for Jaccard Index failure
    if (type == "jaccard" & (threshold > 0.5 | threshold < 0.000001)) {
      message("No reliable knee found-- setting threshold to 0.005")
      safety <- 0.005
    }
    # Safe guard for knee counts failure
    if (type == "bead" & (threshold > 100000 | threshold < 100)) {
      message("No reliable knee found-- setting threshold to 500")
      safety <- 500
    }
    # Safety is with the guard rails; threshold is what the knee calls
    safety <- ifelse(safety > 0, safety, threshold)
    return(c(safety, threshold))
  }

# Import parameters using logic from the end3
nn <- length(args)
csvDir <- args[nn - 10] # directory of csv (previously .rds) files
n_bc_file <- args[nn - 9] # file for the number of reads supporting each barcode
hq_bc_file <- args[nn - 8] # file for the HQ barcodes that were nominated
params_file <- args[nn - 7] # deconvolution params file
tblOut <- args[nn - 6] # filepath to write the implicated barcode pairs
min_jaccard_frag <- as.numeric(args[nn - 5])
name <- args[nn - 4] #name prefix for file naming convention
one_to_one <- args[nn - 3] #arguement for keeping bead : drop conversion 1 to 1
barcoded_tn5 <- args[nn - 2]
ti_len <- as.numeric(args[nn - 1])
barcode_prior_file <- args[nn]

# Fix to R boolean
one_to_one <- one_to_one == "True"
barcoded_tn5 <- barcoded_tn5 == "True"

# Replace the .gz convention
tblOut <- gsub(".gz$", "", tblOut)


csvFiles <- list.files(csvDir, full.names = TRUE, pattern = "_overlapCount.csv.gz$")
csvFiles <- csvFiles[sapply(lapply(csvFiles, file.info), function(x) x$size) > 0]


lapply(csvFiles, fread) %>%
  rbindlist(fill = TRUE) -> inputDF

# Only consider merging when Tn5 is the same
if (barcoded_tn5) {
  tn5_1 <- substrRight(inputDF[["barc1"]], ti_len)
  tn5_2 <- substrRight(inputDF[["barc2"]], ti_len)
  inputDF <- inputDF[tn5_1 == tn5_2]
}


# Import number of barcodes
valid_barcodes <-
  fread(hq_bc_file, col.names = c("bc"), header = FALSE)[["bc"]]
nBC <-
  fread(n_bc_file,
        col.names = c("BeadBarcode", "count"),
        sep = ",") %>%
  data.frame() %>%
   filter(BeadBarcode %in% valid_barcodes) %>%
   arrange(desc(count))
count_vec <-
  nBC$count * 2
names(count_vec) <- as.character(nBC$BeadBarcode)

sum_dt <-
  inputDF[, .(N_both = sum(n_both)), by = list(barc1, barc2)]
sum_dt$N_barc1 <- count_vec[sum_dt$barc1]
sum_dt$N_barc2 <- count_vec[sum_dt$barc2]

data.frame(sum_dt) %>% # fixed the previous divided by two in
# the upstream script (22) for overall accuracy
  mutate(jaccard_frag = round((N_both) / (N_barc1 + N_barc2 - N_both + 0.05), 5)) %>%
  filter(jaccard_frag > 0) %>%
  arrange(desc(jaccard_frag)) %>%
   data.frame() -> ovdf


# Call knee if we need to
if (min_jaccard_frag == 0) {
  message("Computing jaccard index for bead merging via a knee call--")
  two <-
    get_density_threshold_CL(head(ovdf$jaccard_frag, 1000000),
     "jaccard", logTransform = TRUE)
  min_jaccard_frag <- two[1]
  called_jaccard_frag <- two[2]
  # Write out what gets called
  write(paste0(
    "jaccard_threshold_nosafety,",
    as.character(called_jaccard_frag)
  ),
  params_file,
  append = TRUE)
}

# Append to deconvolution parameters
write(paste0("jaccard_threshold,", as.character(min_jaccard_frag)),
      params_file,
      append = TRUE)

ovdf$merged <- ovdf$jaccard_frag > min_jaccard_frag

# Merge barcodes above the knee unless there's a prior reason why we shouldn't
if (barcode_prior_file != "none") {
  # Import and assign each
  bp_df <- fread(barcode_prior_file, header = FALSE)
  feature_vec <-
    as.character(bp_df[["V2"]])
  names(feature_vec) <- as.character(bp_df[["V1"]])
  merged_df <- ovdf %>% filter(merged)
  priorf1 <-
    feature_vec[as.character(merged_df[["barc1"]])] %>% unname()
  priorf2 <-
    feature_vec[as.character(merged_df[["barc2"]])] %>% unname()
  # Find conflicts,  missing values, and
  conflicts <- priorf1 != priorf2
  sum_is_na <- sum(is.na(conflicts))
  sum_conflicts <- sum(conflicts, na.rm = TRUE)
  sum_valid <- sum(!conflicts, na.rm = TRUE)
  # Write out conflicts and filter if we observe them
  if (sum_conflicts > 0) {
    conflictFile <- gsub(
      "/final/",
      "/knee/",
      gsub(
        ".implicatedBarcodes.csv$",
        ".barcode_prior_conflicts.csv",
        tblOut
      )
    )
    write.table(
      merged_df[which(conflicts), ],
      file = conflictFile,
      row.names = FALSE,
      col.names = TRUE,
      sep = ",",
      quote = FALSE
    )
    ovdf <- ovdf[1:dim(ovdf)[1] %ni% which(conflicts), ]
  }
  # Either way, write out statistics
  out_stat_prior_df <- data.frame(
    what = c(
      "Valid merges (stil merged)",
      "Merges with 1 or more NA (still merged)",
      "Conflicted merges (not merged)"
    ),
    count = c(sum_valid, sum_is_na, sum_conflicts)
  )
  statFile <- gsub(
    "/final/",
    "/knee/",
    gsub(
      ".implicatedBarcodes.csv$",
      ".barcode_prior_stats.csv",
      tblOut
    )
  )
  write.table(
    out_stat_prior_df,
    file = statFile,
    row.names = FALSE,
    col.names = FALSE,
    sep = ",",
    quote = FALSE
  )
}

# Now filter based on the min_jaccard_frag
ovdf %>% filter(merged) %>% data.frame() -> ovdf_filt

# Guess at how wide we need to make the barcodes to handle leading zeros
guess <- ceiling(log10(dim(nBC)[1]))

nBC_keep <-
  nBC
nBC_keep$DropBarcode <- ""
nBC_keep$OverlapReads <- ""

# Loop through and eat up barcodes
idx <- 1
while (dim(nBC)[1] > 0) {
  barcode <- as.character(nBC[1, 1])
  barcode_combine <- barcode
  OverlapReads <- 0
  # Find friends that are similarly implicated and append from Barcode 1
  friendsRow1 <- which(barcode ==  ovdf_filt[, "barc1", drop = TRUE])
  if (length(friendsRow1) > 0) {
    friends1 <- as.character(ovdf_filt[friendsRow1, "barc2"])
    barcode_combine <- c(barcode_combine, friends1)
  }
  # Find friends that are similarly implicated and append from Barcode 2
  friendsRow2 <- which(barcode ==  ovdf_filt[, "barc2", drop = TRUE])
  if (length(friendsRow2) > 0) {
    friends2 <- as.character(ovdf_filt[friendsRow2, "barc1"])
    barcode_combine <- c(barcode_combine, friends2)
  }
  # If user species one to one, then only remove that one barcode
  if (one_to_one)
    barcode_combine <- barcode
  # Make a drop barcode and save our progress
  if (!barcoded_tn5) {
    dropBarcode <-
      paste0(
        name,
        "_BC",
        formatC(
          idx,
          width = guess,
          flag = "0",
          digits = 20
        ),
        "_N",
        sprintf("%02d", length(barcode_combine))
      )
  } else {
    dropBarcode <- paste0(
      name,
      "_Tn5-",
      substrRight(
        barcode_combine,
        ti_len
      ),
      "_BC",
      formatC(
        idx,
        width = guess,
        flag = "0",
        digits = 20
      ),
      "_N",
      sprintf("%02d", length(barcode_combine))
    )
  }
  # Annotate with new values
  nBC_keep[nBC_keep$BeadBarcode %in% barcode_combine, "DropBarcode"] <-
    dropBarcode
  idx <- idx + 1
  # Remove barcodes that we've dealt with
  nBC <- nBC[nBC$BeadBarcode %ni% barcode_combine, ]
}

if (barcoded_tn5) {

  params_pattern <- ".deconvolutionParams.orig.csv"
  fastqtis <- gsub(params_pattern, "", list.files(pattern = params_pattern))
  for (fastqti in fastqtis) {
    ti <- gsub("^.*-", "", fastqti)

    # Export the implicated barcodes on a ti level
    ovdfti <- ovdf[grep(paste0(ti, "$"), ovdf$"barc1"), ]
    write.table(
      ovdfti,
      file = paste0(fastqti, ".implicatedBarcodes.csv"),
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE,
      sep = ","
    )
    system(paste0("gzip ", paste0(fastqti, ".implicatedBarcodes.csv")))

    # Export the implicated barcodes
    barctrans <- nBC_keep[, c("BeadBarcode", "DropBarcode")]
    barctransti <- barctrans[grep(paste0(ti, "$"), barctrans$"BeadBarcode"), ]
    write.table(
      barctransti,
      file = paste0(fastqti, ".barcodeTranslate.tsv"),
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t"
    )

  }

} else {

  # Export the implicated barcodes
  write.table(
    ovdf,
    file = tblOut,
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep = ","
  )
  system(paste0("gzip ", tblOut))

  # Export the barcodeTranslate file
  write.table(
    nBC_keep[, c("BeadBarcode", "DropBarcode")],
    file = gsub(".implicatedBarcodes.csv$", ".barcodeTranslate.tsv", tblOut),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )

}
