#!/usr/bin/env Rscript
# Caleb Lareau
# Bio-Rad Laboratories, Inc.

options(warn = -1)
options(datatable.fread.input.cmd.message = FALSE)
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tools)))
suppressMessages(suppressWarnings(library(GenomicRanges)))

"%ni%" <- Negate("%in%")

options(scipen = 999)

# This script processes fragments and annotates them with the bead barcode ID
# based on a separate dictionary file

substrRight <- function(x, n = 6) {
  substr(x, nchar(x) - n + 1, nchar(x))
}

# Input Output
args <- commandArgs(trailingOnly = FALSE)
nn <- length(args)

# Import parameters using logic from the end
frag_bedpe_file <- args[nn - 3]
barcode_translate_file <- args[nn - 2]
frag_anno_file_out <- args[nn - 1]
sumstats_file_out <- args[nn]


# Import frags and annotate with cell barcode ID
frags <-
  fread(
    cmd = paste0("zcat < ", frag_bedpe_file),
    col.names = c(
      "read_name",
      "chr",
      "start",
      "end",
      "orientation",
      "support",
      "bead_barcode"
    ),
    header = FALSE
  )
#remove the orientation column to maintain compatibility with downstream code
frags[, x := NULL]
barcode_translate <-
  fread(
    barcode_translate_file,
    col.names = c("bead_barcode", "cell_barcode"),
    header = FALSE
  )
mdf <- merge(frags, barcode_translate, by = "bead_barcode")


# Group by the cell barcode and 1) get representative read 2) get summary stats


# Part 1 -- deduplicate reads over the cell barcode
pcr_dup_df <- mdf[, list(.I[1], sum(support)), by = list(chr, start, end, cell_barcode)]
pcr_dup_df2 <-
  data.table(pcr_dup_df[, c("chr", "start", "end", "cell_barcode")],
             read_name = mdf$read_name[pcr_dup_df[["V1"]]],
             support = pcr_dup_df[, "V2"])
rm(pcr_dup_df)


# Part 2 -- count fragments for summary statistics
dedup_count <- pcr_dup_df2[, .(n_unique = .N), by = cell_barcode]
raw_count <- mdf[, .(n_total = .N), by = cell_barcode]
merge_ss <- merge(raw_count, dedup_count, by = "cell_barcode")
merge_ss$chr <- as.character(mdf[1, "chr"])

# Export deduplicated fragments
frag_anno_file_out_nocompress <- frag_anno_file_out
pcr_dup_df2 <- pcr_dup_df2[order(start)]
write.table(
  pcr_dup_df2,
  file = frag_anno_file_out_nocompress,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

# Export summary statistics
write.table(
  merge_ss,
  file = sumstats_file_out,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

# Export counts
input_frags <- nrow(frags)
output_frags <- nrow(pcr_dup_df2)
write.table(
  input_frags,
  file = "input_frags.tsv",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
write.table(
  output_frags,
  file = "output_frags.tsv",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
