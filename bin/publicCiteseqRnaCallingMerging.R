#!/usr/bin/env Rscript

# Step 1: Load required libraries
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(Seurat))

# Step 2: Command line argument parsing
parser <- ArgumentParser(
  description = "Merge ADT count data using allowlist and barcode translation"
)
parser$add_argument(
  "--sample_id",
  required = TRUE,
  help = "Sample ID for file naming"
)
parser$add_argument(
  "--allowlist",
  required = TRUE,
  help = "Path to barcode allowlist CSV file"
)
parser$add_argument(
  "--translate",
  required = TRUE,
  help = "Path to barcode translate TSV file"
)
parser$add_argument(
  "--matrix_prefix",
  required = TRUE,
  help = "Prefix for matrix files (without extension)"
)
args <- parser$parse_args()

message("Starting ADT barcode merging and translation...")

# Step 3: Load allowlist file
message("Loading allowlist file...")
allowlist <- fread(args$allowlist, header = FALSE)
colnames(allowlist) <- "DropBarcode"
message(paste("Loaded", nrow(allowlist), "barcodes from allowlist"))

# Step 4: Load translate file
message("Loading barcode translate file...")
translate <- fread(args$translate, header = TRUE)
message(paste("Loaded", nrow(translate), "barcode translations"))

# Step 5: Filter translate file to only include DropBarcodes in allowlist
message("Filtering translate file using allowlist...")
filtered_translate <- translate[translate$DropBarcode %in% allowlist$DropBarcode]
message(paste("Filtered translate table has", nrow(filtered_translate), "entries"))

# Step 6: Load sparse matrix from three files
message("Loading sparse matrix files...")

# Load barcodes
barcodes_file <- paste0(args$matrix_prefix, ".barcodes.txt")
if (!file.exists(barcodes_file)) {
  barcodes_file <- paste0(args$matrix_prefix, ".barcodes.tsv")
}
barcodes <- fread(barcodes_file, header = FALSE)$V1
message(paste(
  "Loaded", length(barcodes), "barcodes"
))

# Load genes
genes_file <- paste0(args$matrix_prefix, ".genes.txt")
if (!file.exists(genes_file)) {
  genes_file <- paste0(args$matrix_prefix, ".genes.tsv")
}
genes <- fread(genes_file, header = FALSE)$V1
message(paste(
  "Loaded", length(genes), "genes"
))

# Load matrix
matrix_file <- paste0(args$matrix_prefix, ".mtx")
sparse_matrix <- readMM(matrix_file)
message(paste(
  "Loaded sparse matrix with dimensions:",
  nrow(sparse_matrix), "x", ncol(sparse_matrix)
))

# Matrix is barcodes x genes - set row and column names
rownames(sparse_matrix) <- barcodes
colnames(sparse_matrix) <- genes

# Step 7: Keep sparse matrix format for memory efficiency
message("Keeping matrix in sparse format for memory efficiency...")
message(paste(
  "Matrix remains sparse - dimensions:",
  nrow(sparse_matrix), "x", ncol(sparse_matrix)
))

# Step 8: Create raw_cells.csv using sparse matrix
message("Creating raw_cells.csv from sparse matrix...")
raw_tib <- as_tibble(as.matrix(sparse_matrix), rownames = NA)
write.csv(raw_tib, "raw_cells.csv", row.names = TRUE)
message("Generated raw_cells.csv")

# Step 9: Generate CellFiltering-style outputs using sparse operations
num_deduped_reads <- sum(Matrix::rowSums(sparse_matrix))
write.csv(data.frame(num_deduped_reads = num_deduped_reads), "deduped_reads.csv")
message("Generated deduped_reads.csv")

umi_percell <- Matrix::rowSums(sparse_matrix)
sorted_umi_percell <- sort(umi_percell, decreasing = TRUE)

# Step 10: Aggregate by DropBarcode BEFORE filtering (unfiltered intermediate matrix)
message("Aggregating by DropBarcode before filtering...")
# Map ALL bead barcodes to DropBarcodes using the full translate table
bead_barcodes <- rownames(sparse_matrix)
dropbarcode_matches <- match(bead_barcodes, translate$BeadBarcode)
dropbarcode_column <- translate$DropBarcode[dropbarcode_matches]

# For unfiltered matrix, we want to include ALL barcodes
# Those without translations get synthetic DropBarcodes
message("Assigning synthetic DropBarcodes to barcodes without translations...")

# Get the starting barcode number: count unique strings in field 2 of translate file
current_barcode <- length(unique(translate[[2]]))
message(paste("Starting barcode number:", current_barcode))

# Vectorized approach: identify which barcodes need synthetic IDs
needs_synthetic <- is.na(dropbarcode_column)
n_synthetic <- sum(needs_synthetic)

if (n_synthetic > 0) {
  # Generate all synthetic barcodes at once
  synthetic_numbers <- seq(current_barcode + 1, current_barcode + n_synthetic)

  # Check for overflow on the last synthetic barcode
  if (
    nchar(as.character(synthetic_numbers[n_synthetic])) >
    nchar(as.character(current_barcode))
  ) {
    message("Barcode value overflow")
  }

  synthetic_barcodes <- paste0(args$sample_id, "BC", synthetic_numbers, "N1")

  # Assign DropBarcodes: use existing ones or synthetic ones
  dropbarcodes_unfiltered <- dropbarcode_column
  dropbarcodes_unfiltered[needs_synthetic] <- synthetic_barcodes
} else {
  # No synthetic barcodes needed
  dropbarcodes_unfiltered <- dropbarcode_column
}

# Keep ALL rows for unfiltered matrix
unfiltered_matrix <- sparse_matrix

# Set DropBarcode as rownames
rownames(unfiltered_matrix) <- dropbarcodes_unfiltered

# Step 11: Sum rows with the same DropBarcode using sparse operations (unfiltered)
message("Summing counts by DropBarcode using sparse operations (unfiltered)...")
unique_dropbarcodes_unfiltered <- unique(dropbarcodes_unfiltered)
n_unique_unfiltered <- length(unique_dropbarcodes_unfiltered)
n_genes <- ncol(unfiltered_matrix)

agg_matrix_unfiltered <- Matrix::sparseMatrix(
  i = match(dropbarcodes_unfiltered, unique_dropbarcodes_unfiltered),
  j = seq_along(dropbarcodes_unfiltered),
  x = 1,
  dims = c(n_unique_unfiltered, length(dropbarcodes_unfiltered))
)

unfiltered_aggregated_sparse <- agg_matrix_unfiltered %*% unfiltered_matrix
rownames(unfiltered_aggregated_sparse) <- unique_dropbarcodes_unfiltered
colnames(unfiltered_aggregated_sparse) <- colnames(unfiltered_matrix)

message(paste(
  "Unfiltered aggregated matrix dimensions:",
  nrow(unfiltered_aggregated_sparse),
  "DropBarcodes x", ncol(unfiltered_aggregated_sparse),
  "genes"
))

# Step 11b: Output unfiltered aggregated matrix in Matrix Market format
message("Saving unfiltered aggregated matrix...")
unfiltered_prefix <- paste0(args$sample_id, ".unfiltered")


# Transpose the unfiltered matrix for output compatibility
unfiltered_agg_transposed <- t(unfiltered_aggregated_sparse)
message(paste(
  "Transposed unfiltered aggregated matrix to:",
  nrow(unfiltered_agg_transposed), "x",
  ncol(unfiltered_agg_transposed)
))

writeMM(unfiltered_agg_transposed, paste0(unfiltered_prefix, ".mtx"))
write.table(
  colnames(unfiltered_agg_transposed),
  paste0(unfiltered_prefix, ".barcodes.tsv"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)
write.table(
  rownames(unfiltered_agg_transposed),
  paste0(unfiltered_prefix, ".genes.tsv"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)
message("Saved unfiltered.mtx, unfiltered.barcodes.tsv, and unfiltered.genes.tsv")

# Step 11c: Now filter the aggregated matrix to only include DropBarcodes in allowlist
message("Filtering aggregated matrix to only include DropBarcodes in allowlist...")
# Convert allowlist to hash set for O(1) lookup performance
allowlist_set <- as.character(allowlist$DropBarcode)
dropbarcodes_in_allowlist <- rownames(unfiltered_aggregated_sparse) %in% allowlist_set
final_matrix_sparse <- unfiltered_aggregated_sparse[
  dropbarcodes_in_allowlist, , drop = FALSE
]

message(paste(
  "Filtered aggregated matrix dimensions:",
  nrow(final_matrix_sparse), "DropBarcodes x", ncol(final_matrix_sparse), "genes"
))

# Convert to regular matrix only for final output
final_matrix <- as.matrix(final_matrix_sparse)

# Step 12: Create Seurat object for filtering simulation
seurat_matrix <- t(final_matrix)
cells <- CreateSeuratObject(seurat_matrix, project = args$sample_id, assay = "ADT")

# Ensure orig.ident is set to the sample ID (though it should already be correct)
cells@meta.data$orig.ident <- args$sample_id

n_initial_cells <- ncol(cells)
message(paste(
  "Created Seurat object with",
  n_initial_cells,
  "cells using processed final_matrix"
))

message(paste(
  "Final matrix dimensions:",
  nrow(final_matrix), "DropBarcodes x", ncol(final_matrix), "genes"
))

# Step 13: Save the final matrix
message("Saving final matrix to: filtered_cells.csv")

# Write as CSV with row names (DropBarcodes) and column names (genes)
write.csv(final_matrix, "filtered_cells.csv", row.names = TRUE)

# Calculate median UMI per cell from final filtered matrix
umi_per_cell_filtered <- rowSums(final_matrix)
median_umi_per_cell <- median(umi_per_cell_filtered)

# Calculate median ADTs per cell from final filtered matrix
adts_per_cell <- apply(final_matrix, 1, function(x) sum(x != 0))
median_adts_per_cell <- median(adts_per_cell)

# Update standardized metrics to include median UMI per cell and median ADTs per cell
standardized_metrics <- data.frame(
  sample = c(args$sample_id, args$sample_id, args$sample_id),
  process = c("cell_filtering", "cell_filtering", "cell_filtering"),
  metric = c("num_deduped_reads", "median_umi_per_cell", "median_adts_per_cell"),
  value = c(num_deduped_reads, median_umi_per_cell, median_adts_per_cell),
  stringsAsFactors = FALSE
)
write.csv(standardized_metrics,
          paste0(args$sample_id, "_cell_filtering.csv"),
          row.names = FALSE,
          quote = FALSE)
message("Updated metrics CSV with median UMI per cell and median ADTs per cell")

# Output filtered count matrix in Matrix Market format (sparse)
filtered_prefix <- paste0(args$sample_id, ".filtered")

# Transpose the filtered matrix for RNA output compatibility
final_matrix_sparse_transposed <- t(final_matrix_sparse)
message(paste(
  "Transposed filtered matrix to:",
  nrow(final_matrix_sparse_transposed), "x",
  ncol(final_matrix_sparse_transposed)
))

writeMM(final_matrix_sparse_transposed, paste0(filtered_prefix, ".mtx"))
write.table(
  colnames(final_matrix_sparse_transposed),
  paste0(filtered_prefix, ".barcodes.tsv"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)
write.table(
  rownames(final_matrix_sparse_transposed),
  paste0(filtered_prefix, ".genes.tsv"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

# Also save summary statistics
summary_file <- paste0(args$sample_id, "_merge_summary.txt")
summary_text <- paste0(
  "ADT Merging Summary for Sample: ", args$sample_id, "\n",
  "================================\n",
  "Original allowlist entries: ", nrow(allowlist), "\n",
  "Original translate entries: ", nrow(translate), "\n",
  "Filtered translate entries: ", nrow(filtered_translate), "\n",
  "Original matrix barcodes: ", length(bead_barcodes), "\n",
  "Barcodes with existing translations: ", sum(!is.na(dropbarcode_column)), "\n",
  "Barcodes with synthetic translations: ", sum(is.na(dropbarcode_column)), "\n",
  "Unfiltered aggregated DropBarcodes: ", nrow(unfiltered_aggregated_sparse), "\n",
  "Final filtered DropBarcodes: ", nrow(final_matrix), "\n",
  "Genes: ", ncol(final_matrix), "\n",
  "Total counts in unfiltered matrix: ", sum(unfiltered_aggregated_sparse), "\n",
  "Total counts in final filtered matrix: ", sum(final_matrix), "\n"
)

writeLines(summary_text, summary_file)
message(paste("Summary saved to:", summary_file))

# Step 14: Save CellFiltering-style outputs
message("Saving CellFiltering-style outputs...")
saveRDS(cells, "cell_calling.RDS")
saveRDS(umi_percell, "umi_percell.RDS")
saveRDS(n_initial_cells, "n_initial_cells.RDS")
message("Saved cell_calling.RDS, umi_percell.RDS, and n_initial_cells.RDS")

message("ADT merging complete!")

message(paste(
  "Final output saved with",
  nrow(final_matrix), "DropBarcodes and", ncol(final_matrix), "genes"
))
