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

# Checking if inputs are set
if (!dir.exists(inputDir)) {
  stop(paste(inputDir, "does not exist."))
}



# Loading dependencies -------------------------------------
#-----------------------------------

library(tidyverse)
library(fastqcr)
library(data.table)
library(R.utils)

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
      map_dbl(read_length_raw, function(x) {
        mean(as.numeric(unlist(
          str_split(x, pattern = "-")
        )))
      })
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
    tibble(sample, tot_seq, avg_length, pct_dup) %>%
    pivot_longer(-c(sample), names_to = "metric", values_to = "value") %>%
    mutate(metric = paste(metric, read, sep = "_"))
  return(summary)
}

# Function for reading mixed species counts RDS file
read_mixed_species_files <- function(input_dir) {
  # Creating mixed species file list
  file_pattern <- "\\.species_mix_counts\\.rds"
  mixed_species_files <-
    list.files(
      path = input_dir,
      pattern = file_pattern,
      full.names = TRUE
    )
  if (length(mixed_species_files) > 0) {
    mixed_species_df <- tibble()
    for (sample_file in mixed_species_files) {
      sample_file_name <- sub("\\.species_mix_counts\\.rds$", "", sample_file) %>%
        substr(4, nchar(.))
      mixed_species_file <- readRDS(sample_file)
      species1 <- mixed_species_file$species1
      species2 <- mixed_species_file$species2
      process_name <- "crosstalk"
      sample <- c(
        sample_file_name, sample_file_name,
        sample_file_name, sample_file_name
      )
      process <- c(process_name, process_name, process_name, process_name)
      metric <- c(
        "mixed_cells", "observed_mixed_species_multiplets",
        "estimated_total_multiplets", "cell_purity"
      )
      # don't include shadow cells in calc if crosstalk threshold = shadow threshold
      cols <- colnames(mixed_species_file$summaryTable)
      shadow_species1 <- paste("Shadow_", species1, "_cells", sep = "")
      shadow_species2 <- paste("Shadow_", species2, "_cells", sep = "")
      if (shadow_species1 %in% cols && shadow_species2 %in% cols) {
        mixed_cells_value <- mixed_species_file$summaryTable["Mixed_cells"] +
          mixed_species_file$summaryTable[shadow_species1] +
          mixed_species_file$summaryTable[shadow_species2]
      } else {
        mixed_cells_value <- mixed_species_file$summaryTable["Mixed_cells"]
      }
      observed_mixed_species_value <-
        mixed_species_file$summaryTable["Measurable_crosstalk"][[1]]
      value <- c(
        mixed_cells_value[[1]], observed_mixed_species_value,
        mixed_species_file$summaryTable["Estimated_total_crosstalk"][[1]],
        mixed_species_file$summaryTable["Cell_purity"][[1]]
      )
      mixed_species_df <- rbind(
        mixed_species_df,
        data.frame(sample, process, metric, value,
          stringsAsFactors = TRUE
        )
      )
    }
    return(mixed_species_df)
  }
}

# Function for reading in all FastQC files
read_fastqc_files <- function(input_dir) {
  # Creating file list
  file_pattern <- "_fastqc\\.zip"
  fastqc_files <-
    list.files(
      path = input_dir,
      pattern = file_pattern,
      full.names = TRUE
    )
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
read_debarcoder_files <- function(input_dir) {
  # Dead metadata has a column for the sample name,
  # so we don't need a per-sample parsing function
  # Creating file lists
  r1_file_pattern <- "_R1_barcode_stats\\.tsv"
  r1_debarcoder_files <-
    list.files(
      path = input_dir,
      pattern = r1_file_pattern,
      full.names = TRUE
    )
  r2_file_pattern <- "_R2_barcode_stats\\.tsv"
  r2_debarcoder_files <-
    list.files(
      path = input_dir,
      pattern = r2_file_pattern,
      full.names = TRUE
    )
  # If files are present
  if (length(r1_debarcoder_files) > 0) {
    # Importing and aggregating files into single df
    debarcoder_df <-
      map_dfr(r1_debarcoder_files, read_tsv, col_types = cols()) %>%
      rename(metric = category) %>%
      rename(value = count) %>%
      add_column(process = "debarcoder", .after = 1) %>% # Denoting process source
      # filter out for just the "total" and "good" metrics
      filter(metric %in% c("total", "good")) %>%
      mutate(metric = case_when(
        metric == "total" ~ "input_r1",
        metric == "good" ~ "output_r1",
        TRUE ~ metric
      ))
    # If there are files for R2
    if (length(r2_debarcoder_files)) {
      r2_debarcoder_df <-
        map_dfr(r2_debarcoder_files, read_tsv, col_types = cols()) %>%
        rename(metric = category) %>%
        rename(value = count) %>%
        add_column(process = "debarcoder", .after = 1) %>% # Denoting process source
        # filter out for just the "total" and "good" metrics
        filter(metric %in% c("total", "good")) %>%
        mutate(metric = case_when(
          metric == "total" ~ "input_r2",
          metric == "good" ~ "output_r2",
          TRUE ~ metric
        ))
      debarcoder_df <- rbind(debarcoder_df, r2_debarcoder_df)
    }
    return(debarcoder_df)
  }
  return(debarcoder_df)
}

# Function for reading in all average length files after Cutadapt
read_cutadapt_trim_files <- function(input_dir) {
  # Creating file list
  file_pattern <- "_cutadapt_trim_r2_length\\.csv"
  length_files <-
    list.files(
      path = input_dir,
      pattern = file_pattern,
      full.names = TRUE
    )
  # If files are present
  if (length(length_files) > 0) {
    # Importing and aggregating files into single df
    length_df <-
      map_dfr(length_files, read_csv, col_types = cols()) %>%
      pivot_longer(avg_length, names_to = "metric", values_to = "value") %>% # Tidying
      mutate(metric = paste(metric, "r2", sep = "_")) %>%
      add_column(process = "cutadapt_trim", .after = 1) # Denoting process source
    return(length_df)
  }
}


# Function for reading in all STAR log files generated during alignment
read_star_align_files <- function(input_dir) {
  # Creating file list
  file_pattern <- "_Log\\.final\\.out"
  star_files <-
    list.files(
      path = input_dir,
      pattern = file_pattern,
      full.names = TRUE
    )
  # If files are present
  if (length(star_files) > 0) {
    # Importing and aggregating files into single df
    star_df <-
      suppressWarnings(
        map_dfr(
          set_names(star_files, star_files),
          read_tsv,
          # Suppressing import warnings because STAR log files have non-standard
          # formatting and naming file vector with itself to create a sample col
          col_names = c("metric", "value"),
          col_types = cols(),
          skip = 5,
          .id = "sample"
        )
      ) %>% # Skipping log header and creating sample col with rows as filenames
      # Converting filename to Nxf sample_id
      mutate(sample = str_replace(sample, paste0(".*/(.*)", file_pattern), "\\1")) %>%
      filter(!is.na(value)) %>% # Removing STAR log section headers imported as df rows
      mutate(
        metric = str_to_lower(metric),
        # Parsing and tidying metric text
        metric = str_replace_all(metric, "number ?|of ", ""),
        metric = str_replace_all(metric, "% | %| \\||,|\\(|\\)", ""),
        metric = str_replace_all(metric, " |-", "_"),
        metric = str_replace_all(metric, ":_", "-"),
        metric = str_replace_all(metric, "/", "."),
        metric = case_when(
          str_detect(value, "%") ~ paste0("pct_", metric),
          # Appending prefix to percent vars
          TRUE ~ metric
        ),
        metric = case_when(
          metric == "input_reads" ~ "filtered_reads",
          TRUE ~ metric
        ),
        value = as.numeric(str_extract(value, "[\\d\\.]+"))
        # Removing any non-numeric characters from value col and reformatting as numeric
      ) %>%
      select(sample, metric, value) %>% # Ordering cols of interest
      add_column(process = "star_align", .after = 1) # Denoting process source
    return(star_df)
  }
}


# Function for reading in all post-alignment metric files produced by Picard
read_picard_files <- function(input_dir) {
  # Creating file list
  file_pattern <- "\\.rna_seq_metrics\\.txt"
  metrics_files <-
    list.files(
      path = input_dir,
      pattern = file_pattern,
      full.names = TRUE
    )
  # If files are present
  if (length(metrics_files) > 0) {
    # Importing and aggregating files into single df
    metrics_df <-
      map_dfr(
        set_names(metrics_files, metrics_files),
        read_tsv,
        # Naming file vector with itself to create a sample col
        col_types = cols(),
        skip = 6,
        n_max = 1,
        .id = "sample"
      ) %>% # Extracting rows of interest and creating sample col with rows as filenames
      # Converting filename to Nxf sample_id
      mutate(sample = str_replace(sample, paste0(".*/(.*)", file_pattern), "\\1")) %>%
      pivot_longer(-sample, names_to = "metric", values_to = "value") %>% # Tidying
      filter(!is.na(value)) %>% # Removing unused metrics
      mutate(
        metric = str_to_lower(metric),
        # Converting decimal proportions to actual percents
        value = case_when(
          str_detect(metric, "pct") ~ value * 100,
          TRUE ~ value
        )
      ) %>%
      select(sample, metric, value) %>%
      add_column(process = "picard", .after = 1) # Denoting process source
    return(metrics_df)
  }
}


# Function for reading in all expression stats files
# nolint start
read_summarize_expression_files <- function(input_dir) {
  # nolint end
  # Creating file list
  # nolint start
  file_pattern_normal <- "\\.final.scrnaseq_counts\\.csv"
  file_pattern_merging <- "\\.final.scrnaseq_counts.csv"
  allowlist_pattern_normal <- "_barcode_allowlist\\.csv"
  allowlist_pattern_merging <- "_updated_barcode_allowlist\\.csv"

  file_regex <- paste(file_pattern_normal,
    file_pattern_merging,
    allowlist_pattern_normal,
    allowlist_pattern_merging,
    sep = "|"
  )
  # nolint end
  count_files_regex <- paste(file_pattern_normal, file_pattern_merging, sep = "|")
  allowlist_files_regex <- paste(allowlist_pattern_normal,
    allowlist_pattern_merging,
    sep = "|"
  )
  count_files <-
    list.files(
      path = input_dir,
      pattern = count_files_regex,
      full.names = TRUE
    )
  allowlist_files <-
    list.files(
      path = input_dir,
      pattern = allowlist_files_regex,
      full.names = TRUE
    )
  # If allowlist files exist
  if (length(allowlist_files) > 0) {
    allow_df <-
      basename(allowlist_files) %>%
      str_remove_all(., file_regex) %>%
      set_names(allowlist_files, .) %>%
      map_dfr(
        .,
        read_csv,
        col_names = "barcode",
        # Naming file vector with itself to create a sample col
        col_types = cols(),
        .id = "sample"
      ) %>%
      mutate(above_knee = TRUE)
  }
  # If files are present
  if (length(count_files) > 0) {
    # Importing and aggregating files into single df
    count_df <- basename(count_files) %>%
      str_remove_all(., file_regex) %>%
      set_names(count_files, .) %>%
      map_dfr(
        .,
        read_csv,
        # Naming file vector with itself to create a sample col
        col_types = cols(),
        .id = "sample"
      ) %>%
      # This empty mutate is to make identical object to old function
      # mutate returns a different type of tibble than map_dfr/read_csv
      mutate(.) %>%
      left_join(allow_df, by = c("sample", "barcode"))

    # Break pipe to filter out below-knee barcodes and count total reads in cells
    reads_in_cells <- count_df %>%
      filter(above_knee) %>%
      group_by(sample) %>%
      summarize(reads_in_cells = sum(input_reads), .groups = "keep")

    # If mixed species, create species-specific reads in cells
    species_cols <- str_subset(names(count_df), "input_.*_reads")
    if (length(species_cols) == 2) {
      species_names <- str_replace_all(species_cols, "input_|_reads", "")
      sp_reads_in_cells <- count_df %>%
        filter(above_knee) %>%
        select(sample, all_of(species_cols)) %>%
        group_by(sample) %>%
        summarize(
          !!paste0(
            species_names[1], "_reads_in_cells"
          ) := sum(!!rlang::sym(species_cols[1])),
          !!paste0(
            species_names[2], "_reads_in_cells"
          ) := sum(!!rlang::sym(species_cols[2])),
          .groups = "keep"
        )
      reads_in_cells <- left_join(reads_in_cells, sp_reads_in_cells)
    }

    # Resume pipe
    count_df <- count_df %>%
      select(-c(barcode, above_knee)) %>% # Dropping unused barcode col
      group_by(sample) %>%
      # Calculating per sample counts across all barcodes
      summarise(across(-matches("sample"), sum), .groups = "keep") %>%
      # Join reads_in_cells
      left_join(reads_in_cells, by = "sample") %>%
      mutate(reads_mapped_transcriptome = genic_reads) %>%
      ungroup() %>%
      pivot_longer(
        -sample,
        names_to = "metric",
        values_to = "value"
      ) %>% # Tidying
      mutate(metric = case_when(
        metric == "input_reads" ~ "deduplicated_reads",
        TRUE ~ metric
      )) %>%
      # Removing gene count because
      # it doesn't deal with redundancy across barcodes
      filter(metric != "genes") %>%
      add_column(
        process = "summarize_expression",
        .after = 1
      ) # Denoting process source

    # Calculating the percent of input reads for each read category
    pct_df <- count_df %>%
      filter(
        !(metric %in%
          c(
            "genes",
            "umi",
            "transcripts",
            str_match(count_df$metric, "input_.*_reads")
          )
        )
      ) %>%
      group_by(sample) %>%
      mutate(
        value = value / value[metric == "deduplicated_reads"] *
          100,
        # Calculating percent based on the number of input reads
        # Changing metric names from raw counts to percent
        metric = paste0("pct_", metric)
      ) %>%
      ungroup() %>%
      filter(metric != "pct_deduplicated_reads") # Removing denominator row

    # Do species specific calculations
    if (length(species_cols) == 2) {
      sp1_df <- count_df %>%
        filter(str_detect(metric, species_names[1]) &
          !str_detect(metric, "umi")) %>%
        group_by(sample) %>%
        mutate(value = value / value[str_detect(metric, "input")] * 100) %>%
        filter(!str_detect(metric, "input")) %>%
        mutate(
          sample = paste0(sample, "_", species_names[1]),
          metric = "fraction_reads_in_cells"
        )
      sp2_df <- count_df %>%
        filter(str_detect(metric, species_names[2]) &
          !str_detect(metric, "umi")) %>%
        group_by(sample) %>%
        mutate(value = value / value[str_detect(metric, "input")] * 100) %>%
        filter(!str_detect(metric, "input")) %>%
        mutate(
          sample = paste0(sample, "_", species_names[2]),
          metric = "fraction_reads_in_cells"
        )
    }

    # Combining dfs together and rename
    # reads_mapped_transcriptome and reads_in_cells
    final_df <- count_df %>%
      bind_rows(pct_df) %>%
      arrange(sample) %>%
      filter(!metric %in% c("reads_in_cells", "reads_mapped_transcriptome")) %>%
      mutate(metric = case_when(
        metric == "pct_reads_mapped_transcriptome" ~
          "reads_mapped_transcriptome",
        metric == "pct_reads_in_cells" ~
          "fraction_reads_in_cells",
        TRUE ~ metric
      ))

    # Bind species-specific reads in cells if they exist
    if (exists("sp1_df")) {
      final_df <- final_df %>%
        bind_rows(sp1_df)
    }
    if (exists("sp2_df")) {
      final_df <- final_df %>%
        bind_rows(sp2_df)
    }
    return(final_df)
  }
}

# Function for reading in all files containing numcell analysis results
read_cell_calling_files <- function(input_dir) {
  # Creating file list
  file_pattern <- "_numcell_analysis\\.csv"
  numcell_files <-
    list.files(
      path = input_dir,
      pattern = file_pattern,
      full.names = TRUE
    )
  # If files are present
  if (length(numcell_files) > 0) {
    # Importing and aggregating files into single df
    numcell_df <-
      map_df(numcell_files, read_csv, col_types = cols()) %>%
      mutate(variable = str_to_lower(variable)) %>% # Making text more parsable
      # Adding metric col for joining with others
      mutate(variable = paste0("called_cells_", variable)) %>% # Tidying
      rename(metric = variable) %>%
      add_column(process = "cell_calling", .after = 1) # Denoting process source
    return(numcell_df)
  }
}


# Function for reading in all files containing expression matrix feature counts
# nolint start
read_count_matrix_features_files <- function(input_dir) {
  # nolint end
  # Creating file list
  file_pattern <- "_matrix_features\\.csv"
  feature_files <-
    list.files(
      path = input_dir,
      pattern = file_pattern,
      full.names = TRUE
    )
  # If files are present
  if (length(feature_files) > 0) {
    # Importing and aggregating files into single df
    feature_df <-
      map_df(feature_files, read_csv, col_types = cols()) %>%
      # Denoting process source
      add_column(process = "count_matrix_features", .after = 1)
    return(feature_df)
  }
}


# Function for reading in all files containing read counts
aggregate_read_counts <- function(input_dir_chr) {
  # Creating file list
  file_pattern <- "_read_counts\\.csv"
  initial_count_files <-
    list.files(
      path = input_dir_chr,
      pattern = file_pattern,
      full.names = TRUE
    )
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
      mutate(metric = ifelse(metric == "read_count",
        paste(metric, read, sep = "_"), metric
      )) %>%
      select(-c("read"))
    return(count_df)
  }
}

# Function for aggregating bead merging data
read_bead_merging <- function(input_dir_chr) {
  # Creating file list
  file_pattern <- "\\.bead_merging_metadata\\.csv"
  merging_files <-
    list.files(
      path = input_dir_chr,
      pattern = file_pattern,
      full.names = TRUE
    )
  # If files are present
  if (length(merging_files) > 0) {
    merging_df <-
      map_dfr(set_names(merging_files, merging_files),
        read_csv,
        col_types = cols(), .id = "sample"
      ) %>%
      mutate(sample = str_replace(sample, paste0(".*/(.*)", file_pattern), "\\1")) %>%
      mutate(percent_raw_unique_deconvolution_reads = round(
        raw_unique_edges / total_edges * 100,
        2
      )) %>%
      mutate(percent_corrected_unique_deconvolution_reads = round(
        corrected_unique_edges / total_edges * 100,
        2
      )) %>%
      mutate(percent_above_knee_unique_deconvolution_reads = round(
        above_knee_unique_edges / total_edges * 100,
        2
      )) %>%
      pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
      add_column(process = "edge_metadata", .after = 1)
    return(merging_df)
  }
}

# Function for reading in all files with cell edge counts
# nolint start
read_cell_edge_count_files <- function(input_dir) {
  # nolint end
  # Creating file list
  file_pattern <- ".cell_edge_counts\\.csv"
  cell_edge_count_files <-
    list.files(
      path = input_dir,
      pattern = file_pattern,
      full.names = TRUE
    )
  # If files are present
  if (length(cell_edge_count_files) > 0) {
    # Importing and aggregating files into single df
    cell_edge_count_df <-
      map_df(cell_edge_count_files, read_csv, col_types = cols())
    return(cell_edge_count_df)
  }
}

# Function for aggregating and calculating per-cell bead merging metrics
aggregate_bead_merging <- function(input_dir_chr) {
  # Creating file list
  # nolint start
  file_pattern_merging <- "\\.final.scrnaseq_counts.csv"
  # nolint end
  count_files <-
    list.files(
      path = input_dir_chr,
      pattern = file_pattern_merging,
      full.names = TRUE
    )
  allowlist_pattern_normal <- "_barcode_allowlist\\.csv"
  allowlist_pattern_merging <- "_updated_barcode_allowlist\\.csv"
  allowlist_files_regex <- paste(allowlist_pattern_normal,
    allowlist_pattern_merging,
    sep = "|"
  )
  allowlist_files <-
    list.files(
      path = input_dir_chr,
      pattern = allowlist_files_regex,
      full.names = TRUE
    )
  if (length(count_files) > 0 && length(allowlist_files) > 0) {
    # Importing and aggregating files into single df
    count_df <-
      map_dfr(
        set_names(count_files, count_files),
        read_csv,
        # Naming file vector with itself to create a sample col
        col_types = cols(),
        .id = "sample"
      ) %>%
      mutate(
        sample = str_replace(sample, paste0(
          ".*/(.*)",
          file_pattern_merging
        ), "\\1")
      )

    allow_df <-
      basename(allowlist_files) %>%
      str_remove_all(., allowlist_files_regex) %>%
      set_names(allowlist_files, .) %>%
      map_dfr(
        .,
        countLines,
        .id = "sample"
      ) %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("sample") %>%
      setNames(., c("sample", "cells"))
    count_df <- count_df %>% filter(count_df$sample %in% allow_df$sample)

    allow_cells_with_merges_files <-
      basename(allowlist_files) %>%
      str_remove_all(., allowlist_files_regex) %>%
      set_names(allowlist_files, .)

    summary_df <- tibble()
    for (sampleId in unique(count_df$sample)) {
      filtered_allow_files <- allowlist_files[grep(sampleId, allowlist_files)]
      filtered_allow_df <- read.table(filtered_allow_files, header = FALSE)
      sample_tr <- gsub("-", "", gsub("_", "", sampleId))
      colnames(filtered_allow_df) <- c("barcode")
      sample_allow_list <- filtered_allow_df[
        grep(sample_tr, filtered_allow_df$barcode),
      ]
      allow_cells_with_merges_calc <-
        sum(
          lapply(
            strsplit(
              str_replace(sample_allow_list, sample_tr, ""), "N"
            ),
            `[[`, 2
          ) > 1
        )
      above_knee <- allow_df$cells[which(allow_df$sample == sampleId)]
      adf <- count_df %>%
        filter(sample == sampleId) %>%
        slice_head(n = above_knee) %>%
        filter(str_detect(barcode, sample_tr)) %>%
        summarize(
          n_cells = above_knee,
          cells_with_merges = allow_cells_with_merges_calc
        ) %>%
        add_column("sample" = sampleId, .before = 1) %>%
        mutate(
          percent_cells_with_multiple_beads =
            round(cells_with_merges / n_cells * 100, 2)
        )
      summary_df <- rbind(summary_df, adf)
    }
    summary_df <- summary_df %>%
      pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
      add_column("process" = "bead_merging", .after = 1)
    return(summary_df)
  }
}

# Function for aggregating median umis per oligo
read_umi_per_oligo <- function(input_dir_chr) {
  # Creating file list
  file_pattern <- "\\.median_umi_per_dimer\\.csv"
  metric_files <-
    list.files(
      path = input_dir_chr,
      pattern = file_pattern,
      full.names = TRUE
    )
  # If metric files are present
  if (length(metric_files) > 0) {
    merging_df <-
      map_dfr(set_names(metric_files, metric_files),
        read_csv,
        col_types = cols(), .id = "sample"
      ) %>%
      mutate(sample = str_replace(sample, paste0(".*/(.*)", file_pattern), "\\1")) %>%
      pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
      add_column(process = "bead_merging", .after = 1)
    return(merging_df)
  }
}

# Function for aggregating bead merging data
read_umi_per_barcode_stats <- function(input_dir_chr) {
  # Creating file list
  file_pattern <- "_umi_per_barcode_stats\\.csv"
  metric_files <-
    list.files(
      path = input_dir_chr,
      pattern = file_pattern,
      full.names = TRUE
    )

  # If files are present
  if (length(metric_files) > 0) {
    merging_df <-
      map_dfr(set_names(metric_files, metric_files),
        read_csv,
        col_types = cols(), .id = "sample"
      ) %>%
      mutate(sample = str_replace(sample, paste0(".*/(.*)", file_pattern), "\\1")) %>%
      pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
      add_column(process = "umi_counting", .after = 1)
    return(merging_df)
  }
}

# Function for aggregating false connection stats
read_false_connections <- function(input_dir_chr) {
  # Creating file list
  file_pattern <- "\\.false_connections\\.csv"
  metric_files <-
    list.files(
      path = input_dir_chr,
      pattern = file_pattern,
      full.names = TRUE
    )
  # If metric files are present
  if (length(metric_files) > 0) {
    merging_df <-
      map_dfr(set_names(metric_files, metric_files),
        read_csv,
        col_types = cols(), .id = "sample"
      ) %>%
      mutate(sample = str_replace(sample, paste0(".*/(.*)", file_pattern), "\\1")) %>%
      pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
      add_column(process = "false_connections", .after = 1)
    return(merging_df)
  }
}

# Function for importing and combining all dfs together
aggregate_data_rna <- function(input_dir_chr) {
  # Creating the final summary df
  list(
    read_fastqc_files(input_dir_chr),
    # Importing all metric files and formatting as a list of dfs
    read_mixed_species_files(input_dir_chr),
    read_cutadapt_trim_files(input_dir_chr),
    read_debarcoder_files(input_dir_chr),
    read_star_align_files(input_dir_chr),
    read_picard_files(input_dir_chr),
    read_summarize_expression_files(input_dir_chr),
    read_cell_calling_files(input_dir_chr),
    read_count_matrix_features_files(input_dir_chr),
    aggregate_read_counts(input_dir_chr),
    read_bead_merging(input_dir_chr),
    read_cell_edge_count_files(input_dir_chr),
    aggregate_bead_merging(input_dir_chr),
    read_umi_per_barcode_stats(input_dir_chr),
    read_false_connections(input_dir_chr),
    read_umi_per_oligo(input_dir_chr)
  ) %>%
    reduce(bind_rows) %>% # Joining all dfs together by appending rowwise
    arrange(sample) %>% # Sorting by sample
    relocate(value, .after = last_col()) # Putting value col last to be human friendly
}

# create demoniator
generate_denominator <- function(data, sample_id, total_reads_metric) {
  denominator <- ((data %>%
    filter(sample == sample_id, metric == total_reads_metric) %>%
    select(value))[[1]])
  return(denominator)
}

# Function precalculate the % no_DO_dimer as needed from the debarcoder stats
precalculate_no_DO_dimer_pct <- function(data, inputDir) {
  # Dead metadata has a column for the sample name,
  # so we don't need a per-sample parsing function
  # Creating file lists
  r1_file_pattern <- "_R1_barcode_stats\\.tsv"
  r1_debarcoder_files <-
    list.files(
      path = inputDir,
      pattern = r1_file_pattern,
      full.names = TRUE
    )
  # If files are present
  if (length(r1_debarcoder_files) > 0) {
    # Importing and aggregating files into single df
    debarcoder_df <-
      map_dfr(r1_debarcoder_files, read_tsv, col_types = cols()) %>%
      rename(metric = category) %>%
      rename(value = count) %>%
      add_column(process = "debarcoder", .after = 1) %>% # Denoting process source
      # filter out for just the "good" and "edge" metrics
      filter(metric %in% c("good", "edge"))
    for (sample_id in unique(debarcoder_df$sample)) {
      # number of total R1 CBC w/ and w/o edges
      denominator <- (debarcoder_df %>%
        filter(sample == sample_id, metric == "good") %>%
        select(value))[[1]]
      # number of total R1 CBC w/ and w/o edges - number w/o edges
      numerator <- denominator - ((debarcoder_df %>%
        filter(sample == sample_id, metric == "edge") %>%
        select(value))[[1]])
      # precalculate the % of no corresponding edges and append to summary df
      final_val <- round((numerator / denominator) * 100, 2)
      pct_df <- tibble(
        sample = sample_id, process = "bead_merging",
        metric = "percent_no_DO_dimer", value = final_val
      )
      data <- rbind(data, pct_df)
    }
  }
  return(data)
}

# extract the value for the specified process/metric and calculate
append_summary_df <- function(
    data, filtered_process, sample_id,
    filtered_metric, denominator, summary_metric,
    summary_df, is_pct, process = "summary") {
  final_val <- calculate_summary(
    data, filtered_process, sample_id,
    filtered_metric, denominator, summary_metric, is_pct
  )
  pct_df <- tibble(
    sample = sample_id, process = process,
    metric = summary_metric, value = final_val
  )
  summary_df <- rbind(summary_df, pct_df)
  return(summary_df)
}

calculate_summary <- function(
    data, filtered_process, sample_id,
    filtered_metric, denominator, summary_metric, is_pct) {
  df <- data %>% filter(process == filtered_process)
  val <- (df %>%
    filter(sample == sample_id, metric == filtered_metric) %>%
    select(value))[[1]]
  if (is_pct == TRUE) {
    final_val <- val / denominator * 100
  } else {
    final_val <- val / denominator
  }
  return(final_val)
}

# extract the value for the specified process/metric
append_pct_summary_df <- function(
    data, filtered_process, sample_id,
    filtered_metric, summary_metric, summary_df) {
  pct <- calculate_pct_summary(
    data, filtered_process, sample_id,
    filtered_metric, summary_metric
  )
  pct_df <- tibble(
    sample = sample_id, process = "summary",
    metric = summary_metric, value = pct
  )
  summary_df <- rbind(summary_df, pct_df)
  return(summary_df)
}

calculate_pct_summary <- function(
    data, filtered_process,
    sample_id, filtered_metric, summary_metric) {
  df <- data %>% filter(process == filtered_process)
  pct <- (df %>%
    filter(sample == sample_id, metric == filtered_metric) %>%
    select(value))[[1]]
  return(pct)
}

estimate_downsampled_cDNA <- function(
    data, sample_id, debarcoder_pct,
    cutadapt_trim_pct, summary_df) {
  # Check if we have a downsampled cDNA
  downsamplecDNA_df <- data %>% filter(process == "downsample_cDNA")
  if (nrow(downsamplecDNA_df) > 0) {
    estimated_raw_downsample_reads <- (downsamplecDNA_df %>%
      filter(sample == sample_id, metric == "input") %>%
      select(value))[[1]] / (cutadapt_trim_pct / 100) / (debarcoder_pct / 100)
    estimated_raw_dsample_reads_df <- tibble(
      sample = sample_id, process = "summary",
      metric = "estimated_raw_cDNA_counts", value = estimated_raw_downsample_reads
    )
    summary_df <- rbind(summary_df, estimated_raw_dsample_reads_df)
    return(summary_df)
  } else {
    return(summary_df)
  }
}

estimate_downsampled_DOS <- function(data, sample_id, debarcoder_pct, summary_df) {
  # Estimate raw DO read input if we have downsampled DOs
  downsampleDOs_df <- data %>% filter(process == "downsample_deconvolution_oligos")
  if (nrow(downsampleDOs_df) > 0) {
    estimated_raw_edges_counts <- (downsampleDOs_df %>%
      filter(sample == sample_id) %>%
      filter(metric == "input") %>%
      select(value))[[1]] / (debarcoder_pct / 100)
    estimated_raw_edges_counts <- round(estimated_raw_edges_counts, digits = 3)
    estimated_raw_edges_counts_df <- tibble(
      sample = sample_id,
      process = "summary", metric = "estimated_raw_DO_counts",
      value = estimated_raw_edges_counts
    )
    summary_df <- rbind(summary_df, estimated_raw_edges_counts_df)
    return(summary_df)
  } else {
    return(summary_df)
  }
}

calculate_estimates <- function(
    data, sample_id, denominator, cutadapt_trim_df,
    total_reads_metric, summary_df) {
  # calculate debarcode pct
  debarcoder_pct <- calculate_summary(
    data, "debarcoder", sample_id,
    "output_r1", denominator, "debarcoder_pct", TRUE
  )
  cutadapt_trim_pct <- 100
  # if cutadapt trim was run calculate pct
  if (nrow(cutadapt_trim_df) > 0) {
    cutadapt_trim_pct <- calculate_summary(
      data, "cutadapt_trim", sample_id,
      "output", denominator, "cutadapt_trim_pct", TRUE
    )
  }
  summary_df <- estimate_downsampled_cDNA(
    data, sample_id, debarcoder_pct,
    cutadapt_trim_pct, summary_df
  )
  summary_df <- estimate_downsampled_DOS(data, sample_id, debarcoder_pct, summary_df)
  estimated_total_reads_count <- (data %>%
    filter(sample == sample_id) %>%
    filter(process == "fastqc") %>%
    filter(metric == total_reads_metric) %>%
    select(value))[[1]] / (cutadapt_trim_pct / 100) / (debarcoder_pct / 100)
  estimated_total_reads_count_df <- tibble(
    sample = sample_id,
    process = "summary", metric = "estimated_total_raw_input_reads",
    value = estimated_total_reads_count
  )
  summary_df <- rbind(summary_df, estimated_total_reads_count_df)
  return(summary_df)
}

# precalculate percents for pipeline summary table
precalculate_pcts <- function(data) {
  # create an empty summary tibble
  summary_df <- tibble()
  # create optional cutadapt_trim_df
  cutadapt_trim_df <- data %>% filter(process == "cutadapt_trim")
  # generate denom and reads metric for fastqc
  fastqc_df <- data %>% filter(process == "fastqc")
  # generate reads per cell denom
  matrix_df <- data %>% filter(process == "count_matrix_features")
  for (sample_id in unique(fastqc_df$sample)) {
    total_reads_metric <- "read_count_r1"
    denominator <- generate_denominator(
      fastqc_df, sample_id,
      total_reads_metric
    )
    summary_df <- calculate_estimates(
      data, sample_id, denominator,
      cutadapt_trim_df, total_reads_metric, summary_df
    )
    # precalculate input read %, input_reads_pct will always equal 100
    summary_df <- append_summary_df(
      data, "fastqc", sample_id,
      total_reads_metric, denominator, "input_pct", summary_df, TRUE
    )
    # precalculate cutadapt_trim_pctm optional output
    if (nrow(cutadapt_trim_df) > 0) {
      summary_df <- append_summary_df(
        data, "cutadapt_trim", sample_id,
        "output", denominator, "cutadapt_trim_pct", summary_df, TRUE
      )
    }
    # precalculate debarcoder_pct
    summary_df <- append_summary_df(
      data, "debarcoder", sample_id,
      "output_r1", denominator, "debarcoder_pct", summary_df, TRUE
    )
    # precalculate star_align_pct
    summary_df <- append_summary_df(
      data, "star_align", sample_id,
      "output", denominator, "star_align_pct", summary_df, TRUE
    )
    # precalculate umi_tools_dedup_pct
    summary_df <- append_summary_df(
      data, "umi_tools_dedup", sample_id,
      "output", denominator, "umi_tools_dedup_pct", summary_df, TRUE
    )
    # grab the already existing summarize_expression %s and add them to the summary df
    summary_df <- append_pct_summary_df(
      data, "summarize_expression", sample_id,
      "pct_intergenic_reads", "reads_mapped_confidently_intergenic", summary_df
    )
    summary_df <- append_pct_summary_df(
      data, "summarize_expression", sample_id,
      "pct_intronic_reads", "reads_mapped_confidently_intronic", summary_df
    )
    summary_df <- append_pct_summary_df(
      data, "summarize_expression", sample_id,
      "pct_coding_reads", "reads_mapped_confidently_exonic", summary_df
    )
    # precalculate valid_barcodes
    # intentional duplicate of debarcoder_pct
    summary_df <- append_summary_df(
      data, "debarcoder", sample_id,
      "output_r1", denominator, "valid_barcodes", summary_df, TRUE
    )
    # grab estimated cells from matrix_df
    estimated_cells <- generate_denominator(matrix_df, sample_id, "estimated_cells")
    # precalculate mean_reads_per_cell
    summary_df <- append_summary_df(
      data, "fastqc", sample_id,
      total_reads_metric, estimated_cells, "mean_reads_per_cell", summary_df, FALSE
    )
  }
  return(summary_df)
}

# Analysis -----------------------------------------------
#-------------------------------------

# Set data.table threading
setDTthreads()

# Importing data and aggregating together

if (pipelineType == "3' RNA Droplet") {
  cat("[PROGRESS]: Importing and aggregating 3' RNA sample metrics.\n")
  data <- aggregate_data_rna(inputDir)
  precalculated_df <- precalculate_pcts(data)
  # check for no DO-dimer output
  file_pattern_no_DO_dimer <- "\\_R1_CBC_no_DO_dimer.tsv"
  no_DO_dimer_files <- list.files(
    path = inputDir,
    pattern = file_pattern_no_DO_dimer,
    full.names = TRUE
  )
  # if no corresponding DO-dimer
  if (length(no_DO_dimer_files) > 0) {
    data <- precalculate_no_DO_dimer_pct(data, inputDir)
  }
  data <- data %>% bind_rows(precalculated_df)
}

# Checking to make sure data exists before writing
if (!is_empty(data)) {
  cat("[PROGRESS]: Writing summary dataframe to file.\n")
  # Writing tidied df to file
  write_csv(data, "metric_summary.csv")
} else {
  cat("[PROGRESS]: No sample metrics files found.\n")
}
