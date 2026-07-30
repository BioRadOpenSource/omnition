#!/usr/bin/env Rscript

library(tidyverse)

comArgs <- commandArgs(TRUE)

input_dir <- comArgs[1]
output_dir <- comArgs[2]
combinatorial <- comArgs[3]

atacPipelineSummaryTable <-
  function(metrics_file,
           sample_id,
           multiplexed,
           ti_config) {
    tryCatch(
      expr = {
        # read metrics file
        metrics <- read_csv(metrics_file,
          col_types = list(
            sample = col_character(),
            process = col_character(),
            metric = col_character(),
            value = col_double()
          )
        )
        # If multiplexed, the sample_id can represent multiple
        # fastq-ti and must be grouped
        if (multiplexed && sample_id == "SuperloadedSample") {
          samples <- unique(metrics$sample)
        } else if (multiplexed &&
                   !sample_id == "SuperloadedSample") {
          samples <-
            ti_config %>%
             filter(sample == sample_id) %>%
             mutate(fastqti = paste(fastq, sequence, sep = "-")) %>%
              pull(fastqti)
        } else {
          samples <- sample_id
        }
        # Input to BWA should be equal to output from DEAD
        bwa_input <- metrics %>%
          filter(process == "bwa", metric == "input", sample %in% samples) %>%
          summarize(bwa_input = sum(value)) %>%
           pull(bwa_input)
        # Output from BWA are properly paired reads
        bwa_output <- metrics %>%
          filter(process == "bwa", metric == "output", sample %in% samples) %>%
          summarize(bwa_output = sum(value)) %>%
           pull(bwa_output)
        # Assemble frags input should be equal to BWA proper pairs output
        assemble_frags_input <- metrics %>%
          filter(str_detect(process, "assemble_fragments"),
                 metric == "input",
                 sample %in% samples) %>%
          summarize(assemble_frags_input = sum(value)) %>%
           pull(assemble_frags_input)
        # BWA output reads are removed in ASSEMBLE_FRAGMENTS if below
        # mapq threshold OR insert size is too large
        assemble_frags_output <- metrics %>%
          filter(
            str_detect(process, "assemble_fragments"),
            metric == "output",
            sample %in% samples
          ) %>%
          summarize(assemble_frags_output = sum(value)) %>%
           pull(assemble_frags_output)
        # ANNOTATE_FRAGMENTS input should be equal to ASSEMBLE_FRAGMENTS output
        annotate_frags_input <- metrics %>%
          filter(
            str_detect(process, "^assemble_fragments"),
            metric == "input",
            sample %in% samples
          ) %>%
          summarize(annotate_frags_input = sum(value)) %>%
           pull(annotate_frags_input)
        # ASSEMBLE_FRAGMENTS output is filtered to remove the blocklist
        # during ANNOTATE_FRAGMENTS
        annotate_frags_output <- metrics %>%
          filter(
            str_detect(process, "^assemble_fragments"),
            metric == "output",
            sample %in% samples
          ) %>%
          summarize(annotate_frags_output = sum(value)) %>%
           pull(annotate_frags_output)
        # REANNOTATE_FRAGMENTS input should be equal to ANNOTATE_FRAGMENTS output
        reanno_frags_input <- metrics %>%
          filter(
            str_detect(process, "reannotate_fragments"),
            metric == "input",
            sample %in% samples
          ) %>%
          summarize(reannotate_frags_input = sum(value)) %>%
           pull(reannotate_frags_input)
        # ANNOTATE_FRAGMENTS output has below knee barcodes and within-droplet
        # duplicates removed during REANNOTATE_FRAGMENTS
        reanno_frags_output <- metrics %>%
          filter(
            str_detect(process, "reannotate_fragments"),
            metric == "output",
            sample %in% samples
          ) %>%
          summarize(reannotate_frags_output = sum(value)) %>%
           pull(reannotate_frags_output)
        # Bring everything together into a data frame
        summary_df <- enframe(
          name = "Metric",
          value = "Value",
          c(
            "Reads with Valid Barcodes" = bwa_input,
            "Pairs Aligned with High Quality" = assemble_frags_output,
            "High Quality Pairs Remaining after Blocklist Removal" =
             annotate_frags_output,
            "Deduplicated Above-Knee Read Pairs" = reanno_frags_output
          )
        ) %>%
          mutate(
            Value = case_when(
              !str_detect(Metric, "Reads with Valid Barcodes") ~
               (Value / Value[Metric == "Reads with Valid Barcodes"]) * 100,
              TRUE ~ Value
            ),
            Value = case_when(
              !str_detect(Metric, "Reads with Valid Barcodes") ~ paste0(format(
                round(Value, 1), nsmall = 1, trim = TRUE, scientific = FALSE
              ), "%"),
              TRUE ~ as.character(Value)
            )
          )
        return(summary_df)
      },
      error = function(e) {
        write(
          paste(
            "[WARNING]: Unable to produce pipeline summary table for",
            sample_id
          ),
          stderr()
        )
        write(e, stderr())
        return(NULL)
      }
    )
  }

batchPipelineSummaryTable <-
  function(multiplex_samples,
           metrics_file,
           multiplexed,
           ti_config) {
    ## create pipeline summary for all samples
    values <- map(
        multiplex_samples,
        atacPipelineSummaryTable,
        metrics_file = metrics_file,
        multiplexed = multiplexed,
        ti_config = ti_config
    )
    ## add sample names to list
    names(values) <- multiplex_samples
    ## remove NULLs from list
    values <- values %>% discard(is.null)
    ## set value column to sample name
    values <-
      map2(values, names(values), function(x, y) {
        names(x)[2] <- y
        return(x)
      })
    ## reduce and join, then transpose
    values <- reduce(values, inner_join, by = "Metric")
    ## convert back to tibble and handle naming
    values <-
      as_tibble(cbind(nms = names(values), t(values))) %>%
       `colnames<-`(c("Sample", .[1, 2:ncol(.)])) %>%
        .[-1, ]
    return(values)
  }

# Calculate reads per index percentages and add new column
readsPerIndex <- function(inputfile) {
  index_df <- read.csv(inputfile)

  index_df <- group_by(index_df, sample) %>%
    mutate(percent = round(count / sum(count) * 100, 1))

  return(index_df)
}

# name of metric summary file
metrics_file <- "metric_summary.csv"
metrics_df <-
  read.csv(metrics_file, stringsAsFactors = FALSE)

# import nextflow parameters
params <- read.csv("tmp/params.csv", header = F, stringsAsFactors = FALSE)
params$V1 <- gsub(" ", "", params$V1)
params <- as.data.frame(t(params), stringsAsFactors = FALSE)
colnames(params) <- as.character(unlist(params[1, ]))
params <- params[-1, ]

####################################
###### PIPELINE SUMMARY TABLE ######
####################################

# Grab full sample names
multiplex_samples <- list.files(input_dir, pattern = "insert_size_metrics.txt")
multiplex_samples <- gsub(".insert_size_metrics.txt", "", multiplex_samples)

## Combinatorial multiplexing
ti_file <- list.files(input_dir, full.names = TRUE, pattern = "*.validated.csv")
if (combinatorial == "true" && (length(file.exists(ti_file)) > 0)) {
    # import the TI config
    ti_config <- read.csv(ti_file)
    # Set multiplexing param
    multiplexed <- TRUE
    ti_df <- metrics_df %>%
              filter(process == "bwa", metric == "input")
    ti_df
    ti_samples <- unique(ti_df$sample)
    ti_samples

## Combinatorial superloading
} else if (combinatorial == "true" && !(length(file.exists(ti_file)) > 0) &&
 "SuperloadedSample" %in% multiplex_samples) {
    multiplex_samples <- "SuperloadedSample"
    ti_config <- NULL
    # Set multiplexing param
    multiplexed <- TRUE
    ti_df <- metrics_df %>%
              filter(process == "bwa", metric == "input")
    ti_df
    ti_samples <- unique(ti_df$sample)
    ti_samples

## Non-combinatorial
} else {
    ti_config <- NULL
    multiplexed <- FALSE
}

# Build aggregated pipeline summary table
tryCatch(
    expr = {
        ## build and write table
        pipeline_summary_table <- batchPipelineSummaryTable(multiplex_samples,
         metrics_file, multiplexed, ti_config)
        pipeline_summary_table
        if (combinatorial == "true") {
            index_pipeline_summary_table <- batchPipelineSummaryTable(ti_samples,
             metrics_file, FALSE, NULL)
            index_pipeline_summary_table
            pipeline_summary_table <-
             rbind(pipeline_summary_table, index_pipeline_summary_table)
        }
        write_csv(pipeline_summary_table, file = "pipeline_summary_table.csv")
    },
    error = function(e) {
        write("[WARNING]: Unable to produce pipeline summary table.", stderr())
        write(e, stderr())
    }
)

####################################
###### PIPELINE SUMMARY TABLE ######
####################################

if (combinatorial == "true") {
    tryCatch(
        expr = {
            index_count_file <- Sys.glob(file.path(getwd(), "fastqTIreadcounts*.csv"))
            reads_per_index_df <- readsPerIndex(index_count_file)
            write_csv(reads_per_index_df,
             file = "fastqTIreadcountsfinal.csv"
            )
        },
        error = function(e) {
            write("[WARNING]: Unable to index count table.", stderr())
            write(e, stderr())
        }
    )
}
