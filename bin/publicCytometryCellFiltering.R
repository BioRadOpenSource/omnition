#!/usr/bin/env Rscript

# setup -------------------------------------------------------------------

# Libraries
library(BUSpaRse)
library(tidyverse)
library(Seurat)
library(argparse)
library(data.table)
library(Matrix)

# read in command line arguments
parser <- ArgumentParser(description = "Process a run report")
parser$add_argument("--inputdir",
                   help = "input directory of results to process")
parser$add_argument("--outputdir",
                   help = "directory to receive output")
parser$add_argument("--sample",
                   help = "Name of sample currently being analyzed")
parser$add_argument("--matrix_prefix",
                   help = "prefix for sparse matrix files")
opts <- parser$parse_args()
message(opts)

# Setup Directories
data_dir <- opts$inputdir
out_dir <- opts$outputdir

# pipeline ----------------------------------------------------------------

##### Read in allowlist ####
allowlist_path <- file.path(data_dir,
                           paste0(opts$sample, "_barcode_allowlist.csv"))
allowlist <- fread(allowlist_path, header = FALSE)

##### Read counter data ####
raw_cells <- BUSpaRse::read_count_output(opts$inputdir,
                                        opts$matrix_prefix,
                                        tcc = FALSE)

#### Raw data exploration ####
# Retain row names
raw_tib <- as.matrix(raw_cells) %>% t() %>% as_tibble(rownames = NA)
write.csv(raw_tib, file.path(out_dir, "raw_cells.csv"), row.names = TRUE)

num_deduped_reads <- sum(rowSums(raw_tib))
write.csv(data.frame(num_deduped_reads = num_deduped_reads),
          file.path(out_dir, "deduped_reads.csv"))

umi_percell <- rowSums(raw_tib)
sorted_umi_percell <- sort(umi_percell, decreasing = TRUE)
umi_percell_indexed <- tibble(index = seq_along(sorted_umi_percell),
                             nUMI = sorted_umi_percell)

#### Create the Seurat object ####
# Create with minimal features to avoid RNA naming, then add ADT assay properly
cells <- CreateSeuratObject(raw_cells, project = opts$sample, assay = "ADT")

# Print the number of cells we're starting with
message(paste("Number of cells in initial Seurat object: ",
              ncol(cells), sep = ""))
n_initial_cells <- ncol(cells)

#### Filter for cells on allowlist ####
message("Filtering for cells on allowlist")
cell_filter <- !is.na(match(cells@assays$ADT@data@Dimnames[[2]],
                           allowlist$V1))
cells <- cells[, cell_filter]
message(paste("Number of cells remaining after allowlist subsetting: ",
              ncol(cells), sep = ""))

# Save filtered cells in same format as raw_cells.csv
filtered_tib <- as.matrix(cells@assays$ADT@data) %>%
  t() %>%
  as_tibble(rownames = NA)
write.csv(filtered_tib,
          file.path(out_dir, "filtered_cells.csv"),
          row.names = TRUE)

# Calculate median UMI per cell from filtered count matrix
umi_per_cell_filtered <- rowSums(filtered_tib)
median_umi_per_cell <- median(umi_per_cell_filtered)

# Calculate median ADTs per cell from filtered count matrix
# For each cell (row), count number of nonzero ADTs (columns)
adts_per_cell <- apply(filtered_tib, 1, function(x) sum(x != 0))
median_adts_per_cell <- median(adts_per_cell)

# Update standardized metrics to include median UMI per cell and median ADTs per cell
standardized_metrics <- data.frame(
  sample = c(opts$sample, opts$sample, opts$sample),
  process = c("cell_filtering", "cell_filtering", "cell_filtering"),
  metric = c("num_deduped_reads", "median_umi_per_cell", "median_adts_per_cell"),
  value = c(num_deduped_reads, median_umi_per_cell, median_adts_per_cell),
  stringsAsFactors = FALSE
)
write.csv(standardized_metrics,
          file.path(out_dir, paste0(opts$sample, "_cell_filtering.csv")),
          row.names = FALSE,
          quote = FALSE)

# Save filtered count matrix in sparse matrix format (same as coreBuildCountMatrix.R)
message("Saving filtered count matrix in sparse matrix format")
filtered_matrix <- cells@assays$ADT@data
barcodes <- colnames(filtered_matrix)
features <- rownames(filtered_matrix)

# Write sparse matrix files with same prefix format
matrix_prefix <- paste0(opts$sample, ".filtered")
writeMM(filtered_matrix, file = file.path(out_dir, paste0(matrix_prefix, ".mtx")))
write.table(barcodes,
            file = file.path(out_dir, paste0(matrix_prefix, ".barcodes.tsv")),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
write.table(features,
            file = file.path(out_dir, paste0(matrix_prefix, ".genes.tsv")),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write_rds(cells, file.path(out_dir, "cell_calling.RDS"))
write_rds(umi_percell, file.path(out_dir, "umi_percell.RDS"))
write_rds(n_initial_cells, file.path(out_dir, "n_initial_cells.RDS"))
