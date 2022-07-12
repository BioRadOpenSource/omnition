#!/usr/bin/env Rscript
# aggregateMetrics.R
# Bio-Rad Laboratories, Inc.

# Purpose: Aggregate and parse sample metrics as a tidy df
# Usage: aggregateMetrics.R INPUTDIR

# Setting environment --------------------------------------
# -----------------------------------


# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
inputDir <- args[1]
pipelineType <- args[2]
mito_contig <- args[3]
ti <- args[4]

# Checking if inputs are set
if (!dir.exists(inputDir)) {
  stop(paste(inputDir, "does not exist."))
}



# Loading dependencies -------------------------------------
#-----------------------------------

library(tidyverse)
library(fastqcr)
library(data.table)

# Functions ------------------------------------------------
#-----------------------------------


# Function for reading in a single FastQC file
read_fastqc <- function(fastqc_zip) {
  # Retrieving data from specific modules
  modules <-
    c(
      "summary",
      "statistics",
      "Sequence Duplication Levels",
      "Sequence Length Distribution"
    )
  # Parsing FastQC zip file
  data <- qc_read(fastqc_zip, modules = modules, verbose = FALSE)
  # Getting file name
  file <- data$basic_statistics %>%
    filter(Measure == "Filename") %>%
    pull(Value) %>%
    str_replace("\\.fastq.*", "")
  if (str_detect(file, "_R[12]_")) {
    # Getting Nxf formatted sample_id
    sample <- str_replace(file, "(.*)_R[12]_.*", "\\1")
    # Getting read number
    read <-
      str_to_lower(str_replace(file, ".*_(R[12])_.*", "\\1"))
  } else {
    # Getting Nxf formatted sample_id
    sample <- str_replace(file, "(.*)_[12]\\..*", "\\1")
    # Getting read number
    read <-
      str_to_lower(str_replace(file, ".*_([12])\\..*", "\\1"))
  }
  # Extracting total number of sequences
  tot_seq <- data$basic_statistics %>%
    filter(Measure == "Total Sequences") %>%
    pull(Value) %>%
    as.numeric()
  # Averaging lengths if given as a range
  read_length_raw <- data$sequence_length_distribution$Length
  if (is.character(read_length_raw)) {
    read_length <-
      map_dbl(read_length_raw, function(x)
        mean(as.numeric(unlist(
          str_split(x, pattern = "-")
        ))))
  } else {
    read_length <- read_length_raw
  }
  # Calculating average read length
  avg_length <-
    weighted.mean(read_length, data$sequence_length_distribution$Count)
  # Calculating percent of duplication
  pct_dup <- 100 - data$total_deduplicated_percentage
  # Joining metrics into tidy df
  summary <-
    tibble(sample, read, tot_seq, avg_length, pct_dup) %>%
    pivot_longer(-c(sample, read), names_to = "metric", values_to = "value")
  return(summary)
}


# Function for reading in all FastQC files
read_fastqc_files <- function(input_dir) {
  # Creating file list
  file_pattern <- "_fastqc\\.zip"
  fastqc_files <-
    list.files(path = input_dir,
               pattern = file_pattern,
               full.names = TRUE)
  # If files are present
  if (length(fastqc_files) > 0) {
    # Importing and aggregating files into single df
    fastqc_df <-
      map_dfr(fastqc_files, read_fastqc) %>% # Using bespoke import function
      add_column(process = "fastqc", .after = 1) # Denoting process source
    return(fastqc_df)
  }
}

# Function for reading input, output, and bad reads from DEAD barcode stats
read_dead_files <- function(input_dir) {
  # Dead metadata has a column for the sample name,
  # so we don't need a per-sample parsing function
  # Creating file lists
  r1_file_pattern <- "_R1_barcode_stats\\.tsv"
  r1_dead_files <-
    list.files(path = input_dir,
               pattern = r1_file_pattern,
               full.names = TRUE)
  r2_file_pattern <- "_R2_barcode_stats\\.tsv"
  r2_dead_files <-
    list.files(path = input_dir,
               pattern = r2_file_pattern,
               full.names = TRUE)
  # If files are present
  if (length(r1_dead_files) > 0) {
    # Importing and aggregating files into single df
    dead_df <-
      map_dfr(r1_dead_files, read_tsv, col_types = cols()) %>%
      rename(metric = category) %>%
      add_column(read = "r1", .after = 1) %>%  # Adding col to specify read no.
      add_column(process = "dead", .after = 1) %>% # Denoting process source
      add_column(value = NA, .after = 4) %>% # Hide non-customer facing fields
      # filter out for just the "total" and "good" metrics
      filter(metric %in% c("total", "good")) %>%
      mutate(metric = case_when(
        metric == "total" ~ "input",
        metric == "good" ~ "output",
        TRUE ~ metric
      ))
    # If there are files for R2, i.e. an ATAC TI run
    if (length(r2_dead_files)) {
      r2_dead_df <-
        map_dfr(r2_dead_files, read_tsv, col_types = cols()) %>%
        rename(metric = category) %>%
        add_column(read = "r2", .after = 1) %>%  # Adding col to specify read no.
        add_column(process = "dead", .after = 1) %>% # Denoting process source
        add_column(value = NA, .after = 4) %>% # Hide non-customer facing fields
        # filter out for just the "total" and "good" metrics
        filter(metric %in% c("total", "good")) %>%
        mutate(metric = case_when(
          metric == "total" ~ "input",
          metric == "good" ~ "output",
          TRUE ~ metric
        ))
      dead_df <- rbind(dead_df, r2_dead_df)
    }
    return(dead_df)
  }
  return(dead_df)
}


# Function for reading in all average length files
# after Cutadapt headcrop
read_cutadapt_headcrop_files <- function(input_dir) {
  # Creating file list
  file_pattern <- "_cutadapt_headcrop_r2_length\\.csv"
  length_files <-
    list.files(path = input_dir,
               pattern = file_pattern,
               full.names = TRUE)
  # If files are present
  if (length(length_files) > 0) {
    # Importing and aggregating files into single df
    length_df <-
      map_dfr(length_files, read_csv, col_types = cols()) %>%
      pivot_longer(avg_length, names_to = "metric", values_to = "value") %>% # Tidying
      add_column(read = "r2", .after = 1) %>%  # Adding col to specify read no.
      add_column(process = "cutadapt_headcrop", .after = 1) # Denoting process source
    return(length_df)
  }
}


# Function for reading in all files containing read counts
aggregate_read_counts <- function(input_dir_chr) {
  # Creating file list
  file_pattern <- "_read_counts\\.csv"
  initial_count_files <-
    list.files(path = input_dir_chr,
               pattern = file_pattern,
               full.names = TRUE)
  # If files are present
  if (length(initial_count_files) > 0) {
    count_files <- list()
    for (file in initial_count_files) {
      if (length(readLines(file)) > 1) {
        count_files <- append(count_files, file)
      }
    }
    count_df <-
      map_df(count_files, read_csv, col_types = cols()) %>%
      arrange(sample) %>% # Sorting by sample
      relocate(count, .after = last_col()) # Putting count col last to be human friendly
    return(count_df)
  }
}

# Function for reading in alignment QC metrics from
# Picard CollectAlignmentSummaryMetrics
#nolint start
read_alignment_qc_metrics_files <- function(input_dir_chr) {
  #nolint end
  # Creating file list
  file_pattern <- "alignment_summary_qc\\.txt"
  metrics_files <-
    list.files(path = input_dir_chr,
               pattern = file_pattern,
               full.names = TRUE)
  # If files are present
  if (length(metrics_files) > 0) {
    # Importing and aggregating files into single df
    metrics_df <-
      map_dfr(
        set_names(
          metrics_files,
          getSampleId(metrics_files, ".alignment_summary_qc.txt")
        ),
        read_tsv,
        # Naming file vector with itself to create a sample col
        col_types = cols(),
        skip = 6,
        n_max = 3,
        .id = "sample"
      ) %>% # Extracting rows of interest and creating sample col with rows as filenames
      mutate(sample = str_replace(sample, paste0(".*/(.*)",
      # Converting filename to Nxf sample_id
                                                 file_pattern), "\\1")) %>%
      pivot_longer(-c(sample, CATEGORY),
                   names_to = "metric",
                   values_to = "value") %>% # Tidying
      mutate(
        read = case_when(
          .$CATEGORY == "FIRST_OF_PAIR" ~ "r1",
          .$CATEGORY == "SECOND_OF_PAIR" ~ "r2",
          .$CATEGORY == "PAIR" ~ "pair"
        )
      ) %>%
      filter(!is.na(value)) %>% # Removing unused metrics
      mutate(
        metric = str_to_lower(metric),
        # Converting decimal proportions to actual percents
        value = case_when(str_detect(metric, "pct") ~ value * 100,
                          TRUE ~ value)
      ) %>%
      select(sample, read, metric, value) %>%
      add_column(process = "picard_alignment_qc", .after = 1) # Denoting process source
    return(metrics_df)
  }
}


# Function for parsing Samtools flagstat output
parse_flagstat <- function(filename) {
  #parse text file
  flagstat <- read.table(filename, fill = T, as.is = TRUE)
  flagstat[, 1] <- as.numeric(flagstat[, 1])
  #write to something easier to comprehend
  df <- data.frame(
    "sample" = str_replace(basename(filename), "(.*)_flagstat.*", "\\1"),
    "in_total" = flagstat$V1[1],
    "secondary" = flagstat$V1[2],
    "supplementary" = flagstat$V1[3],
    "duplicates" = flagstat$V1[4],
    "mapped" = flagstat$V1[5],
    "paired in sequencing" = flagstat$V1[6],
    "read1" = flagstat$V1[7],
    "read2" = flagstat$V1[8],
    "properly_paired" = flagstat$V1[9],
    "with_itself_and_mate_mapped" = flagstat$V1[10],
    "singletons" = flagstat$V1[11],
    "with_mate_mapped_to_a_different_chr" = flagstat$V1[12],
    "with_mate_mapped_to_a_different_chr_mapQ>=5" = flagstat$V1[13],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  return(df)
}

# Function for aggregating Samtools flagstat output
read_flagstat_files <- function(input_dir_chr) {
  # Creating file list
  file_pattern <- "_flagstat\\.txt"
  flagstat_files <-
    list.files(path = input_dir_chr,
               pattern = file_pattern,
               full.names = TRUE)
  # If files are present
  if (length(flagstat_files) > 0) {
    flagstat_df <- map_df(flagstat_files, parse_flagstat) %>%
      pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
      filter(!is.na(value)) %>%
      mutate(metric = str_to_lower(metric)) %>%
      add_column(process = "samtools_flagstat", .after = 1)
    return(flagstat_df)
  }
}

parse_idxstats <- function(filename) {
  idxstats <-
    read_tsv(
      filename,
      col_names = c(
        "contig",
        "contig_length",
        "reads_aligned",
        "reads_unaligned"
      )
    )
  idxstats$sample <-
    str_replace(basename(filename), "(.*)_idxstats.txt", "\\1")
  return(idxstats)
}

# Function for aggregating Samtools idxstats output
read_idxstats_files <- function(input_dir_chr) {
  # Creating file list
  file_pattern <- "_idxstats\\.txt"
  idxstats_files <-
    list.files(path = input_dir_chr,
               pattern = file_pattern,
               full.names = TRUE)
  # If files are present
  if (length(idxstats_files) > 0) {
    idxstats_df <- map_df(idxstats_files, parse_idxstats) %>%
      pivot_longer(-c(sample, contig),
                   names_to = "metric",
                   values_to = "value") %>%
      unite(metric, c(contig, metric)) %>%
      select(sample, metric, value) %>%
      add_column(process = "samtools_idxstats", .after = 1)
    return(idxstats_df)
  }
}

parse_mark_duplicates <- function(filename) {
  #parse text file
  metrics <-
    suppressWarnings(as.data.frame(t(fread(filename,
                                           skip = "## METRICS")),
                                            stringsAsFactors = FALSE)) %>%
    rownames_to_column() %>%
    as_tibble() %>%
    mutate(sample = str_replace(basename(filename), "(.*).dedup_stats.*", "\\1")) %>%
    rename(.,  metric = rowname, value = V1)
  return(metrics)
}

read_mark_duplicates_files <- function(input_dir_chr) {
  # Creating file list
  file_pattern <- "dedup_stats\\.txt"
  dupes_files <-
    list.files(path = input_dir_chr,
               pattern = file_pattern,
               full.names = TRUE)
  # If files are present
  if (length(dupes_files) > 0) {
    dupes_df <- map_df(dupes_files, parse_mark_duplicates) %>%
      select(sample, metric, value) %>%
      add_column(process = "picard_markduplicates", .after = 1) %>%
      filter(metric != "LIBRARY") %>%
      mutate(value = as.numeric(value)) %>%
      mutate(
        metric = str_to_lower(metric),
        # Converting decimal proportions to actual percents
        value = case_when(str_detect(metric, "percent") ~ value * 100,
                          TRUE ~ value)
      )
    return(dupes_df)
  }
}

parse_deconvolution_files <- function(sample_id, input_dir_chr) {
  barcode_quants <-
    fread(paste0(input_dir_chr, "/", sample_id, ".barcodeQuantSimple.csv"))
  deconvolution_params <-
    fread(paste0(input_dir_chr, "/", sample_id, ".deconvolutionParams.csv"))
  qc_stats <-
    fread(paste0(input_dir_chr, "/", sample_id, ".QCstats.csv"))
  jaccard_frag_vec <-
    fread(
      cmd = paste0(
        "zcat < ",
        input_dir_chr,
        "/",
        sample_id,
        ".implicatedBarcodes.csv.gz",
        " | head -1000000"
      )
    )[["jaccard_frag"]]
  # Order based on V2 from highest to lowest
  setorder(barcode_quants, -V2)
  # Set fragment_thresh to the bead threshold from deconvolution_params
  fragment_thresh <-
    deconvolution_params$V2[deconvolution_params$V1 == "bead_threshold"]
  # Gets the total number of barcodes in barcode quants that are above fragment_thresh
  above_knee <- barcode_quants[V2 >= fragment_thresh, .N]
  # Set jaccard_thresh to the jaccard threshold from deconvolution_params
  jaccard_thresh <-
    deconvolution_params$V2[deconvolution_params$V1 == "jaccard_threshold"]
  # Gets the total number of barcodes in jaccard_frag_vec that are above jaccard_thresh
  above_jaccard_knee <-
    length(jaccard_frag_vec[jaccard_frag_vec >= jaccard_thresh])

  # Conditionals to find the right alignment file depending on TI configuration
  if (file.exists("fastqTIreadcounts.csv")) {
    # If TIs are used in any configuration besides superloading
    mapping_file <- read.csv("fastqTIreadcounts.csv", stringsAsFactors = FALSE)

    fastq_and_seq <- strsplit(sample_id, split = "-")

    total_read_pairs <- subset(mapping_file,
                               fastq == fastq_and_seq[[1]][1] &
                               sequence == fastq_and_seq[[1]][2])$count

  } else if (file.exists("fastqTIreadcounts_superloaded.csv")) {
    # If TIs are used in the superloading configuration
    mapping_file <- read.csv("fastqTIreadcounts_superloaded.csv",
                             stringsAsFactors = FALSE)

    fastq_and_seq <- strsplit(sample_id, split = "-")

    total_read_pairs <- subset(mapping_file,
                               fastq == fastq_and_seq[[1]][1] &
                               sequence == fastq_and_seq[[1]][2])$count

  } else {
    # If TIs are not used
    # Getting info from alignment
    metrics_file <- paste0(input_dir_chr, "/", sample_id, ".alignment_summary_qc.txt")

    # Naming file vector with itself to create a sample col
    metrics_df <- map_dfr(set_names(metrics_file, sample_id), read_tsv,
                          # Extracting rows of interest
                          # and creating sample col with rows as filenames
                          col_types = cols(), skip = 6, n_max = 3, .id = "sample") %>%
                    # Converting filename to Nxf sample_id
                    mutate(sample = str_replace(sample,
                                                paste0(".*/(.*)",
                                                       "alignment_summary_qc\\.txt"),
                                                "\\1")) %>%
                        # Tidying
                        pivot_longer(-c(sample, CATEGORY),
                                     names_to = "metric", values_to = "value")

    # Getting the total read pairs from the alignment metrics file
    total_read_pairs <- metrics_df %>%
    filter(CATEGORY == "PAIR" & metric == "TOTAL_READS") %>%
    pull(value) / 2
  }

  deconvolution_stats <- tibble(
    sample = str_replace(basename(
      paste0(sample_id, ".barcodeQuantSimple.csv")
    ), "(.*).barcodeQuantSimple.*", "\\1"),
    total_barcodes_observed = nrow(barcode_quants),
    beads_above_knee = above_knee,
    fragment_thresh_at_bead_knee = fragment_thresh,
    median_frags_per_above_knee_bead = barcode_quants[V2 >= fragment_thresh,
     median(V2)],
    bead_pairs_above_jaccard_knee = above_jaccard_knee,
    jaccard_thresh_at_knee = jaccard_thresh,
    total_cells = nrow(qc_stats),
    total_read_pairs_per_total_cells = round(total_read_pairs / nrow(qc_stats)),
    median_unique_frags_per_cell = median(qc_stats$uniqueNuclearFrags),
    median_estimated_library_size_per_cell = as.double(median(qc_stats$librarySize))
  )
  return(deconvolution_stats)
}

read_deconvolution_files <- function(input_dir_chr) {
  # Creating file list for multiple inputs
  quants_pattern <- "barcodeQuantSimple\\.csv"
  quants_files <-
    list.files(path = input_dir_chr,
               pattern = quants_pattern,
               full.names = TRUE)
  deconvolution_params_pattern <- "deconvolutionParams\\.csv"
  deconvolution_params_files <-
    list.files(path = input_dir_chr,
               pattern = deconvolution_params_pattern,
               full.names = TRUE)
  deconvolution_qc_pattern <- "QCstats\\.csv"
  deconvolution_qc_files <-
    list.files(path = input_dir_chr,
               pattern = deconvolution_qc_pattern,
               full.names = TRUE)
  jaccard_pattern <- "implicatedBarcodes\\.csv\\.gz"
  jaccard_files <-
    list.files(path = input_dir_chr,
               pattern = jaccard_pattern,
               full.names = TRUE)
  # Group files by sample ID
  quants_samples <-
    unlist(lapply(quants_files, getSampleId, ".barcodeQuantSimple.csv"))
  params_samples <-
    unlist(lapply(
      deconvolution_params_files,
      getSampleId,
      ".deconvolutionParams.csv"
    ))
  qc_samples <-
    unlist(lapply(deconvolution_qc_files, getSampleId, ".QCstats.csv"))
  jaccard_samples <-
    unlist(lapply(jaccard_files, getSampleId, ".implicatedBarcodes.csv.gz"))
  # Get sample IDs that have all files
  samples <-
    intersect(jaccard_samples, intersect(quants_samples,
     intersect(params_samples, qc_samples)))
  if (length(samples) > 0) {
    deconvolution_df <-
      map_dfr(samples, parse_deconvolution_files, input_dir_chr) %>%
      pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
      select(sample, metric, value) %>%
      add_column(process = "deconvolution", .after = 1)
    return(deconvolution_df)
  }
}

getSampleId <- function(path, extension) {
  return(basename(path) %>% str_replace(paste0(
    "(.*)", tools::file_path_sans_ext(extension), ".*"
  ), "\\1"))
}

parse_tss_file <- function(tss_file) {
  tss_enrichment <- fread(tss_file)
  tss_score <- round(max(tss_enrichment), 1)
  return(tibble(
    sample = getSampleId(tss_file, ".tss_enrichment.csv"),
    tss_score = tss_score
  ))
}

read_tss_files <- function(input_dir_chr) {
  # Creating file list
  file_pattern <- "tss_enrichment\\.csv"
  tss_files <-
    list.files(path = input_dir_chr,
               pattern = file_pattern,
               full.names = TRUE)
  # If files are present
  if (length(tss_files) > 0) {
    tss_df <-
      map_dfr(set_names(tss_files, getSampleId(tss_files,
                                               ".tss_enrichment.csv")),
                                                parse_tss_file, .id = "sample") %>%
      pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
      add_column(process = "tss_enrichment", .after = 1) %>%
      select(sample, process, metric, value)
    return(tss_df)
  }
}

read_frip_files <- function(input_dir_chr) {
  # Creating file list
  file_pattern <- "frip\\.csv"
  frip_files <-
    list.files(path = input_dir_chr,
               pattern = file_pattern,
               full.names = TRUE)
  if (length(frip_files) > 0) {
    frip_df <- map_dfr(frip_files, read_csv, col_types = cols())
  }
  return(frip_df)
}

read_frit_files <- function(input_dir_chr) {
  # Creating file list
  file_pattern <- "frit\\.csv"
  frit_files <-
    list.files(path = input_dir_chr,
               pattern = file_pattern,
               full.names = TRUE)
  if (length(frit_files) > 0) {
    frit_df <- map_dfr(frit_files, read_csv, col_types = cols())
  }
  return(frit_df)
}


read_atac_crosstalk_files <- function(input_dir_chr) {
  # Creating file list
  crosstalk_file_pattern <- "crosstalk\\.csv"
  crosstalk_files <-
    list.files(path = input_dir_chr,
               pattern = crosstalk_file_pattern,
               full.names = TRUE)
  # If files are present
  if (length(crosstalk_files) > 0) {
    crosstalk_df <-
      map_dfr(set_names(
        crosstalk_files,
        getSampleId(crosstalk_files, ".crosstalk.csv")
      ), read_csv, .id = "sample") %>%
      mutate(measurable_crosstalk =
               as.numeric(str_replace(Measurable_crosstalk, " %", ""))) %>%
      select(-Measurable_crosstalk) %>%
      mutate(estimated_total_crosstalk =
               as.numeric(str_replace(Estimated_total_crosstalk, " %", ""))) %>%
      select(-Estimated_total_crosstalk) %>%
      pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
      mutate(metric = str_to_lower(metric)) %>%
      select(sample, metric, value) %>%
      add_column(process = "crosstalk", .after = 1)
    return(crosstalk_df)
  }
}

parse_peaks <- function(peaks_file) {
  return(tibble(
    sample = getSampleId(peaks_file, ".fixedwidthpeaks.bed"),
    peaks_called = as.numeric(system(
      paste("cat", peaks_file, "| wc -l"), intern = T
    ))
  ))
}

read_peak_calls_files <- function(input_dir_chr) {
  # Creating file list
  peaks_file_pattern <- "fixedwidthpeaks\\.bed"
  peak_files <-
    list.files(path = input_dir_chr,
               pattern = peaks_file_pattern,
               full.names = TRUE)
  # If files are present
  if (length(peak_files) > 0) {
    peaks_df <-
      map_dfr(set_names(
        peak_files,
        getSampleId(peak_files, ".fixedwidthpeaks.bed")
      ), parse_peaks, .id = "sample") %>%
      pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
      select(sample, metric, value) %>%
      add_column(process = "peak_calling", .after = 1)
    return(peaks_df)
  }
}

# Function for importing and combining all dfs together
aggregate_data_atac <- function(input_dir_chr) {
  # Creating the final summary df
  agg_df <- list(
    read_fastqc_files(input_dir_chr),
    # Importing all metric files and formatting as a list of dfs
    read_cutadapt_headcrop_files(input_dir_chr),
    read_dead_files(input_dir_chr),
    read_alignment_qc_metrics_files(input_dir_chr),
    read_flagstat_files(input_dir_chr),
    read_idxstats_files(input_dir_chr),
    read_mark_duplicates_files(input_dir_chr),
    read_deconvolution_files(input_dir_chr),
    read_tss_files(input_dir_chr),
    read_frip_files(input_dir_chr),
    read_frit_files(input_dir_chr),
    read_atac_crosstalk_files(input_dir_chr),
    read_peak_calls_files(input_dir_chr),
    aggregate_read_counts(input_dir_chr)
  ) %>%
    reduce(bind_rows) %>% # Joining all dfs together by appending rowwise
    arrange(sample) %>% # Sorting by sample
    relocate(value, .after = last_col()) # Putting value col last to be human friendly

  # The following chunk of code is for calculating summary metrics for use in the report
  if (ti) {
    for (sampfile in list.files(pattern = ".sampleBeadFiltSummary.csv")) {
      samp <- gsub(".sampleBeadFiltSummary.csv", "", sampfile)
      sample_bead_filt_table <- read.csv(sampfile)
      # add total cell summary metric to df
      agg_df <- agg_df %>% add_row(
        sample = samp,
        process = "summary",
        metric = "total_cells",
        value = sample_bead_filt_table$TotalCells
      )
      read_counts <- read.csv(
        list.files(pattern = "fastqTIreadcounts"), stringsAsFactors = FALSE)
      # Total read pairs is the sum of counts for this sample name
      reads_per_cell <-  read_counts %>%
        filter(sample == samp) %>%
        summarize(count = sum(count)) %>%
        pull(count) /
        sample_bead_filt_table$TotalCells
      # add reads per cell summary metric to df
      agg_df <- agg_df %>% add_row(
        sample = samp,
        process = "summary",
        metric = "mean_valid_read_pairs_per_cell",
        value = reads_per_cell
      )
      # add uniq nuc fragments per cell summary metric to df
      agg_df <- agg_df %>% add_row(
        sample = samp,
        process = "summary",
        metric = "median_unique_nuclear_fragments_per_cell",
        value = sample_bead_filt_table$MedianCellFrags
      )
      # filter the fastqtireadcounts file to the current sample
      # and create a list of the relevant QC stats files
      samp_to_fastqti <-  read_counts %>%
        filter(sample == samp)
      fastqtis <- paste(
        samp_to_fastqti$fastq, "-", samp_to_fastqti$sequence, ".QCstats.csv", sep = "")
      # read in all the QC stats files and bind them together
      sample_level_inf <- fastqtis %>%
        lapply(read_csv) %>%
        bind_rows
      # add estimated lib size summary metric to df
      agg_df <- agg_df %>% add_row(
        sample = samp,
        process = "summary",
        metric = "median_estimated_library_size_per_cell",
        value = as.double(median(sample_level_inf$librarySize))
      )
    }
  # When not using TIs do the following
  } else {
    for (samp in unique(
      agg_df[agg_df$process == "deconvolution", ]$sample)) {
      # calculates the frip score for the summary process for the reports
      tot_cells <- agg_df[agg_df$sample == samp &
        agg_df$metric == "total_cells", "value"]
      # add total cells summary metric to df
      agg_df <- agg_df %>% add_row(
        sample = samp,
        process = "summary",
        metric = "total_cells",
        value = tot_cells[[1]]
      )
      # calculates the frip score for the summary process for the reports
      read_per_cell <- agg_df[agg_df$sample == samp &
        agg_df$metric == "total_read_pairs_per_total_cells", "value"]
      # add reads per cell summary metric to df
      agg_df <- agg_df %>% add_row(
        sample = samp,
        process = "summary",
        metric = "mean_valid_read_pairs_per_cell",
        value = read_per_cell[[1]]
      )
      # calculates the frip score for the summary process for the reports
      nuc_frag_per_cell <- agg_df[agg_df$sample == samp &
        agg_df$metric == "median_unique_frags_per_cell", "value"]
      # add uniq nuc fragments per cell summary metric to df
      # (this code is out of place to preserve formatting)
      agg_df <- agg_df %>% add_row(
        sample = samp,
        process = "summary",
        metric = "median_unique_nuclear_fragments_per_cell",
        value = nuc_frag_per_cell[[1]]
      )
      # calculates the frip score for the summary process for the reports
      nuc_frag_per_cell <- agg_df[agg_df$sample == samp &
        agg_df$metric == "median_estimated_library_size_per_cell", "value"]
      # add estimated lib size summary metric to df
      agg_df <- agg_df %>% add_row(
        sample = samp,
        process = "summary",
        metric = "median_estimated_library_size_per_cell",
        value = nuc_frag_per_cell[[1]]
      )
    }
  }
  # Other summary metrics can be done same way for samples w/ and w/o TIs
  for (samp in unique(agg_df[agg_df$process == "picard_alignment_qc", ]$sample)) {
    # get total reads and /100 so we get percents directly
    tot_reads <- agg_df[
        agg_df$sample == samp &
        agg_df$metric == "total_reads" &
        agg_df$read == "pair",
        "value"
      ] /
      100
    # calculates reads aligned for the summary process for the reports
    reads_aligned <- agg_df[
        agg_df$sample == samp &
        agg_df$metric == "reads_aligned_in_pairs" &
        agg_df$read == "pair",
        "value"
      ] /
      tot_reads
    # add aligned read percent summary metric to df
    agg_df <- agg_df %>% add_row(
      sample = samp,
      process = "summary",
      metric = "aligned_read_pair_pct",
      value = reads_aligned[[1]]
    )
    # calculates duplicate reads aligned for the summary process for the reports
    reads_dup <- agg_df[
        agg_df$sample == samp &
        agg_df$metric == "read_pair_duplicates",
        "value"
      ] *
      # need *2 because orig values is in pairs
      2 /
      tot_reads
    # add duplicate read summary metric to df
    agg_df <- agg_df %>% add_row(
      sample = samp,
      process = "summary",
      metric = "aligned_duplicate_read_pair_pct",
      value = reads_dup[[1]]
    )
    # calculates mito reads aligned for the summary process for the reports
    reads_mito <- agg_df[
        agg_df$sample == samp &
        agg_df$process == "samtools_idxstats",
      ] %>%
      filter(str_detect(metric,
                        paste0("^(.*\\.)?", mito_contig, "_reads_(un)?aligned")))  %>%
      summarize(value = sum(value)) %>% # adding up all mitochondrial reads
      pull(value) /
      tot_reads[[1]]
    # add mito reads summary metric to df
    agg_df <- agg_df %>% add_row(
      sample = samp,
      process = "summary",
      metric = "aligned_mitochondrial_read_pair_pct",
      value = reads_mito
    )
  }
  for (samp in unique(agg_df[agg_df$process == "frip", ]$sample)) {
    # calculates the frip score for the summary process for the reports
    frip_score <- agg_df[agg_df$sample == samp &
                    agg_df$metric == "readsInPeaks", "value"] /
                    agg_df[agg_df$sample == samp &
                    agg_df$metric == "totalDedupedReads", "value"]
    # add frip summary metric to df
    agg_df <- agg_df %>% add_row(
      sample = samp,
      process = "summary",
      metric = "fraction_of_reads_in_peaks",
      value = frip_score[[1]]
    )
  }
  agg_df <- arrange(agg_df, sample)
  return(agg_df)
}


# Analysis -----------------------------------------------
#-------------------------------------

# Importing data and aggregating together

if (pipelineType == "ATAC") {
  cat("[PROGRESS]: Importing and aggregating ATAC sample metrics.\n")
  data <- aggregate_data_atac(inputDir)
}


# Checking to make sure data exists before writing
if (!is_empty(data)) {
  cat("[PROGRESS]: Writing summary dataframe to file.\n")
  # Writing tidied df to file
  write_csv(data, "metric_summary.csv")
} else {
  cat("[PROGRESS]: No sample metrics files found.\n")
}
