#!/usr/bin/env Rscript
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
blocklist_file <- args[nn - 4]
frag_bedpe_file <- args[nn - 3]
read_bead_file <- args[nn - 2]
annotated_out_file <- args[nn - 1]
unique_count_out_file <- args[nn]

# Set data.table threading
setDTthreads()

# Import frags and annotate with bead
frags <-
  fread(
    cmd = paste0("zcat < ", frag_bedpe_file),
    col.names = c("chr", "start", "end", "read_name", "orientation", "support")
  )
bead_read <-
  fread(
    cmd = paste0("zcat < ", read_bead_file),
    col.names = c("read_name", "bead_id")
  ) %>%
   na.omit() %>%
    unique()
mdf <- merge(frags, bead_read, by = "read_name") %>% na.omit()

# Filter for fragments overlapping the blocklist
bl <- fread(blocklist_file, col.names = c("chr", "start", "end"))
if (nrow(bl) > 0) {
  bl <- bl %>% data.frame() %>% makeGRangesFromDataFrame()
  ov_bl <- findOverlaps(bl, makeGRangesFromDataFrame(mdf))
  blocklist_reads <- mdf$read_name[subjectHits(ov_bl)]
  mdf <- mdf[read_name %ni% blocklist_reads]
}

# NOTE: the majority of the missing read names are instances where the
# MAPQ filter was not surpassed.

# Quantify the number of unique fragments per barcode
pcr_dup_df <-
  mdf[, .(count = .N), by = list(chr, start, end, bead_id)]
out_bead_quant <- pcr_dup_df[, .(nUnique = .N), by = bead_id]

# Export tables
write.table(
  mdf,
  file = annotated_out_file,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
write.table(
  out_bead_quant,
  file = unique_count_out_file,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

# Export read counts
input_frags <- nrow(frags)
output_frags <- nrow(mdf)
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
