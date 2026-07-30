#!/usr/bin/env Rscript

#Use of this source code is governed by an MIT
#license that can be found in the LICENSE file or at
#https://opensource.org/licenses/MIT.

# Code: basic pipeline for expression matrix analysis

# setup -------------------------------------------------------------------

# Libraries
library(BUSpaRse)
library(tidyverse)
library(Seurat)
library(GGally)
library(svglite)
library(gridExtra)
library(pheatmap)
library(argparse)
library(rjson)
library(grid)
library(ape)
library(purrr)
library(data.table)

# read in command line arguments
parser <- ArgumentParser(description = "Process a run report")
parser$add_argument("--inputdir",
                    help = "input directory of results to process",
                    required = TRUE)
parser$add_argument("--outputdir",
                    help = "directory to receive output",
                    required = TRUE)
parser$add_argument("--sample",
                    help = "Name of sample currently being analyzed")
opts <- parser$parse_args()

print(opts)

# Setup Directories
out_dir <- opts$outputdir

# Set defaults
dpi <- 18
max_dim <- 10 # number PCs for clustering and dimension reduction
pairwise_plot_markers <- c("CD3", "CD4", "CD8", "CD45RA", "CD45RO", "CD45",
                          "CD19", "CD11b", "CD14", "CD16")
tsne_aspect_ratio <- 1.5 # Aspect ratio for large tsne plots

# functions ---------------------------------------------------------------

# Common ggplot themes
histogram_theme <- theme(axis.title = element_text(size = 16),
                        axis.ticks.y = element_blank(),
                        plot.title = element_text(size = 18),
                        panel.background = element_rect(fill = "white"),
                        panel.grid.major.x = element_blank(),
                        panel.grid.major.y = element_line(colour = "#EBEBEB"),
                        plot.margin = margin(10, 10, 10, 10))

violin_theme <- theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
                     axis.text.y = element_text(size = 16, hjust = 1),
                     axis.title = element_text(size = 16),
                     axis.ticks.x = element_blank(),
                     panel.background = element_rect(fill = "white"),
                     panel.grid.major.x = element_blank(),
                     panel.grid.major.y = element_line(colour = "#EBEBEB"))

# Helper function to convert matrix to tibble
matrix_to_tibble <- function(matrix_data) {
  matrix_data %>% t() %>% as_tibble(rownames = NA)
}

# Helper function to check if matrix has data
has_matrix_data <- function(matrix_data) {
  !is.null(matrix_data) && nrow(matrix_data) > 0
}

# Helper function to extract ADT assay data
extract_adt_data <- function(cells_obj) {
  adt_assay <- cells_obj@assays$ADT

  if (has_matrix_data(adt_assay@scale.data)) {
    return(matrix_to_tibble(adt_assay@scale.data))
  }

  if (has_matrix_data(adt_assay@data)) {
    return(matrix_to_tibble(adt_assay@data))
  }

  return(matrix_to_tibble(adt_assay@counts))
}

# Helper function to extract default assay data
extract_default_assay_data <- function(cells_obj) {
  default_assay <- DefaultAssay(cells_obj)
  assay_data <- cells_obj@assays[[default_assay]]

  if (has_matrix_data(assay_data@scale.data)) {
    return(matrix_to_tibble(assay_data@scale.data))
  }

  if (has_matrix_data(assay_data@data)) {
    return(matrix_to_tibble(assay_data@data))
  }

  return(matrix_to_tibble(assay_data@counts))
}

# Function to extract scaled data as tibble (used multiple times)
extract_scaled_data <- function(cells_obj) {
  # Check if ADT assay exists
  adt_exists <- "ADT" %in% names(cells_obj@assays)

  if (adt_exists) {
    return(extract_adt_data(cells_obj))
  }

  # If no ADT assay, use the default assay
  return(extract_default_assay_data(cells_obj))
}

# Function to create feature plots for a single antibody
create_feature_plots <- function(cells, feature_name, tsne_success, umap_success, dpi) {
  tryCatch({
    if (tsne_success) {
      FeaturePlot(cells, features = feature_name, min.cutoff = 0, max.cutoff = "q95",
                 pt.size = 0.5, reduction = "tsne", cols = c("blue", "red"))
      ggsave(paste0("tsne.", feature_name, ".png"),
            width = 4, height = 4, dpi = dpi * 4)
    }
    if (umap_success) {
      FeaturePlot(cells, features = feature_name, min.cutoff = 0, max.cutoff = "q95",
                 ncol = 3, pt.size = 0.5, reduction = "umap", cols = c("blue", "red"))
      ggsave(paste0("umap.", feature_name, ".png"),
            width = 4, height = 4, dpi = dpi * 4)
    }
    # Downsample to 1000 cells if we have more, otherwise use all
    downsample_n <- min(1000, ncol(cells))
    VlnPlot(subset(cells, downsample = downsample_n), features = feature_name, ncol = 3)
    ggsave(paste0("vln.", feature_name, ".png"), dpi = dpi * 3)
    RidgePlot(cells, features = feature_name, ncol = 2)
    ggsave(paste0("ridge.", feature_name, ".png"), dpi = dpi * 3)
  }, error = function(e) {
    print(paste("Warning: Failed to create feature plots for", feature_name))
    print(e)
  })
}

# Helper function to write CSV files with consistent parameters
write_csv_file <- function(data, filepath, with_rownames = FALSE) {
  write.table(data, file = filepath, row.names = with_rownames,
             col.names = TRUE, sep = ",")
}

# Function to check if any string in list1 is contained in any string in list2
contains_marker <- function(column_name, markers) {
  any(sapply(markers, function(marker) grepl(marker, column_name)))
}

# This extracts dimensionality reduction matrices from Seurat object
get_reductions <- function(reduction, obj) {
  obj@reductions[[reduction]]@cell.embeddings[, 1:2] %>%
    data.table::as.data.table(., keep.rownames = "barcode", key = "barcode")
}

MakePairwisePlot <- function(dat, to_save, outfile_name = "plot.png",
                            start_col = 1, end_col = ncol(dat),
                            width = 15, height = 15, dpi = 300) {
    # Makes pairwise biaxial plots specific to the format of the cells tibble
    # Args:
    #   dat: seurat object of interest. Assumed particular configuration
    #   to_save: boolean indicating whether to save the plot of interest
    #   outfile_name: the name of the file to be saved. Otherwise ignore.
    #   start_col: the column to start the pairwise plots on
    #   end_col: the column to end the pairwise plots on
    # Returns:
    #   the plot of interest
    # Note:
    #   The format as of now has the first column and last 4 columns omitted
    plt <- ggpairs(dat[, c(start_col:end_col)])

    if (to_save) {
        ggsave(file = outfile_name,
               plot = plt,
               width = width,
               height = height,
               dpi = dpi)
    }

    return(plt)
}

FeatureSummarize <- function(dat) {
    # Takes a tibble of cells by features as input as returns a data frame
    #   summarizing expression information
    # Args: dat: tibble of cells by features
    # Returns: data frame of features by min, 1st quartile, median, mean,
    #   3rd quartile, and max
    # Note: function assumes that the barcodes are in the first column of the
    #   data frame

    # Make the list of summary information
    if (ncol(dat) < 2) {
      warning("FeatureSummarize: Data has fewer than 2 columns")
      return(data.frame())
    }
    feature_cols <- dat[, 2:ncol(dat)]
    ab_sums <- lapply(feature_cols, summary)

    # Convert to data frame more efficiently
    ab_sums <- do.call(rbind, ab_sums)
    return(as.data.frame(ab_sums))
}

# Take the scaled count matrix, a list of all the antibodies to make
# histograms from, and a list of unscaled antibodies to make all-zero
# histograms for
AntibodyHistograms <- function(scaled_matrix, all_antibodies,
                              zero_antibodies) {
  ab_histograms <- lapply(all_antibodies, function(antibody_name) {
    if (!(antibody_name %in% zero_antibodies)) {
      # If the antibody is one that has scaled data
      antibody <- reshape2::melt(scaled_matrix[antibody_name],
                                id.vars = NULL)
      ggplot(antibody, aes(x = value)) +
        geom_histogram(fill = "#1E5FAA", bins = 30) +
        ggtitle(antibody_name) +
        xlab("Expression") +
        ylab("Count") +
        histogram_theme
    } else {
      # If the antibody is not scaled and thus needs to be approximated with
      # all zeroes
      artificial_column <- tibble(rep(0, nrow(scaled_matrix)))
      colnames(artificial_column) <- c(antibody_name)
      antibody <- reshape2::melt(artificial_column, id.vars = NULL)
      ggplot(antibody, aes(x = value)) +
        scale_x_continuous(limits = c(-0.1, 1)) +
        geom_histogram(fill = "#1E5FAA") +
        ggtitle(antibody_name) +
        xlab("Expression") +
        ylab("Count") +
        histogram_theme
    }
  })
  return(ab_histograms)
}

# pipeline ----------------------------------------------------------------

# Read in post cell calling filtered results
cells <- readRDS("cell_calling.RDS")
umi_percell <- readRDS("umi_percell.RDS")
n_initial_cells <- readRDS("n_initial_cells.RDS")

# Record the number post calling but pre selection
n_post_calling <- ncol(cells)

# Always further subset by UMI count, in most cases this probably won't remove
# any more cells
print(paste0("max nUMI per cell is ", max(umi_percell),
            ", proceeding with min umi cutoff"))

# margin = 2 means normalize by columns (cells)
# FindVariableFeatures can fail with very few features (e.g., only 3 ADTs)
# because loess smoothing requires sufficient data points
n_features <- nrow(cells)
print(paste0("Number of features: ", n_features))

variable_features_success <- FALSE
if (n_features >= 10) {
  # Only attempt FindVariableFeatures if we have enough features
  tryCatch({
    cells <- FindVariableFeatures(cells)
    top10 <- head(VariableFeatures(cells), 10)
    LabelPoints(plot = VariableFeaturePlot(cells), points = top10,
               repel = TRUE, xnudge = 0, ynudge = 0)
    ggsave("variable_features_plot.png", dpi = dpi)
    variable_features_success <- TRUE
  }, error = function(e) {
    print("Warning: FindVariableFeatures failed, likely due to too few features")
    print(e)
    print("Continuing with all features as variable features")
  })
} else {
  print("Too few features for FindVariableFeatures, using all features")
}

cells <- ScaleData(cells, do.center = FALSE, do.scale = FALSE)

# Determine the number of PCs to compute based on available features
n_features <- nrow(cells)
n_cells <- ncol(cells)
# PCA requires at least 2 features
# and the number of PCs must be less than both features and cells
max_pcs <- min(50, n_features - 1, n_cells - 1)
max_pcs <- max(1, max_pcs)  # Ensure at least 1 PC

print(paste0("Computing PCA with ", max_pcs, " components (",
             n_features, " features, ", n_cells, " cells available)"))

pca_success <- FALSE
if (n_features >= 2 && n_cells >= 2) {
  tryCatch({
    cells <- RunPCA(cells, features = rownames(cells), npcs = max_pcs)
    pca_success <- TRUE
  }, error = function(e) {
    print("Warning: PCA failed")
    print(e)
    # Set max_pcs to 0 to indicate no PCA available
    max_pcs <<- 0
  })
} else {
  print(paste0("Insufficient data for PCA (need at least 2 features and 2 cells)"))
  max_pcs <- 0
}

# Since the cells Seurat object is umi cut as soon as it's created, we don't
# get a Seurat object with all cells
# Write object to file
saveRDS(cells, file = "umi_cut.cells.rds")

# Create tibble from processed Seurat object
out <- extract_scaled_data(cells)
# This output isn't actually cd45 filtered
write.csv(out, "umi_cut.cells.csv")

names(out) <- sub("/", "-", names(out)) # ggplot doesn't seem to like the "/"

# Antibody strength
ab_sums <- FeatureSummarize(out)
ab_sums <- ab_sums[order(ab_sums$`3rd Qu.`, decreasing = TRUE), ]
tryCatch({
  pheatmap(ab_sums[, c(-1, -6)],
           filename = "ab.staining.summary.heatmap.png",
           fontsize_row = 5,
           fontsize_col = 10,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = opts$sample)
}, error = function(e) {
  print("Warning: Heatmap generation failed")
  print(e)
})

# Antibody strength violin plot, not currently displayed
max_violin_per_plot <- 13
antibodies <- colnames(out)
for (x in seq(1, ncol(out), by = max_violin_per_plot)) {
  end_col <- min(x + max_violin_per_plot - 1, ncol(out))
  ab_violin_data <- reshape2::melt(out[antibodies[x:end_col]],
                                  id.vars = NULL)
  ab_violin_plot <- ggplot(ab_violin_data, aes(x = variable, y = value)) +
                    ylim(0, 2) +
                    xlab("Antibody") +
                    ylab("Expression") +
                    geom_violin(trim = TRUE, fill = "#1E5FAA") +
                    violin_theme
  filename <- paste0("ab.staining.violin.", str_pad(x, 3, pad = "0"), ".png")
  ggsave(filename = filename, plot = ab_violin_plot, width = 20, height = 5,
         dpi = dpi)
}

# Get the list of antibodies that won't get a histogram in the above loop
# because they don't have data in the out object due to low counts
zero_antibodies <- setdiff(rownames(cells), colnames(out))
# Antibody strength histograms, same data as violin plot but offers more
# precise information. Currently displayed
all_antibodies <- rownames(cells)
ab_histograms <- AntibodyHistograms(out, all_antibodies, zero_antibodies)
# Arrange the histograms into a grid
ab_hist_plot <- do.call(grid.arrange, c(ab_histograms, c(ncol = 5)))
hist_height <- 4 * ceiling(length(antibodies) / 5)
ggsave(filename = "ab.staining.histograms.png", plot = ab_hist_plot,
       width = 20, height = hist_height, limitsize = FALSE, dpi = dpi * 3)

# PCA/information checking
if (pca_success) {
  tryCatch({
    ElbowPlot(cells) # Sharp drop at 10
    ggsave("pca_plot.png", dpi = dpi)
  }, error = function(e) {
    print("Warning: ElbowPlot failed")
    print(e)
  })
} else {
  print("Skipping PCA plots - PCA was not successful")
}

tryCatch({
  results_table <- read.csv(paste0(opts$sample, "_numcell_analysis.csv"),
                           header = TRUE)
  sample_filter <- results_table$sample == opts$sample &
                  results_table$variable == "kneedle_cumfrac_knee"
  kneedle <- results_table$value[sample_filter]
  cumfrac_file <- paste0(opts$sample, "_allalgos_cumfrac.csv")
  kneedle_cumfrac_frame <- read.csv(cumfrac_file, header = TRUE,
                                   stringsAsFactors = FALSE)

  # Get cumulative fraction plot
  # Show the calculated knee points on a plot
  cumfrac_plot <- ggplot(data = kneedle_cumfrac_frame,
                        aes(x = `Cell.Index`, y = `Cumulative.Fraction`)) +
    geom_line() +
    geom_vline(xintercept = kneedle) +
    # geom_text(aes(x = kneedle, y = 0, label = kneedle),
    #          angle = 90, size = 3, vjust = -0.5, hjust = 0) +
    labs(title = paste0(opts$sample, " Cumulative Fraction"))
  print(cumfrac_plot)
  cumfrac_filename <- paste0(opts$sample, "_allalgos_cumfrac.png")
  ggsave(cumfrac_filename, path = opts$outputdir)
  cumfrac_csv <- file.path(opts$outputdir,
                          paste0(opts$sample, "_allalgos_cumfrac.csv"))
  write_csv_file(kneedle_cumfrac_frame, cumfrac_csv)

  # Get cumulative log log fraction plot
  # Show the calculated knee points on a plot
  log10numi_frame <- read.csv(paste0(opts$sample, "_allalgos_loglog.csv"),
                             header = TRUE)
  log10_plot <- ggplot(data = log10numi_frame,
                      aes(x = `Log10.Cell.Index`, y = `Log10.nUMI`)) +
    geom_line() +
    geom_vline(data = results_table,
              mapping = aes(xintercept = log10(kneedle))) +
    geom_text(data = results_table,
             mapping = aes(x = log10(results_table$value), y = 0,
                          label = results_table$variable),
             angle = 90, size = 3, vjust = -0.5, hjust = 0) +
    labs(title = paste0(opts$sample, " Log10 nUMI"))
  loglog_filename <- paste0(opts$sample, "_allalgos_loglog.png")
  ggsave(loglog_filename, path = opts$outputdir)

  loglog_csv <- file.path(opts$outputdir,
                         paste0(opts$sample, "_allalgos_loglog.csv"))
  write_csv_file(log10numi_frame, loglog_csv)
}, error = function(e) {
  print("Warning: Could not generate cumulative fraction or loglog plots")
  print("Missing input files (numcell_analysis allalgos_cumfrac or allalgos_loglog)")
  print(e)
})

print(paste0("Number of cells remaining after filtering: ", ncol(cells)))
n_post_filtering <- ncol(cells)
# Write a csv file with the post selection and post filtering counts
select_filtering_table <- data.frame(
  sample = rep(opts$sample, 3),
  variable = c("n_initial_cells", "n_post_calling", "n_post_filtering"),
  value = c(n_initial_cells, n_post_calling, n_post_filtering)
)
select_filter_csv <- file.path(opts$outputdir, "select_filtering.csv")
write_csv_file(select_filtering_table, select_filter_csv)

# Generate the scaled matrix after processing
filtered_out <- extract_scaled_data(cells)

# Pairwise plotting
v <- sapply(colnames(filtered_out), contains_marker,
           markers = pairwise_plot_markers)
# Select pairwise plot - only if we have matching markers
if (sum(v) >= 2) {
  tryCatch(
    MakePairwisePlot(dat = filtered_out[, names(v)[v]],
                     to_save = TRUE,
                     outfile_name = "pairwise.plot.select.png",
                     width = 15,
                     height = 15,
                     dpi = 72),
    error = function(e) {
      print("Error in creation of pairwise plot")
      print(e)
    }
  )
} else {
  print("Skipping pairwise plot - insufficient matching markers found")
}

# Clustering
# Seurat RunTSNE will fail if the max dims given is greater than the number of
# dimensions computed via PCA
# So, if max_dims is too big, we reduce it to the number of dimensions
# calculated by PCA
if (pca_success) {
  computed_dims <- ncol(Embeddings(object = cells[["pca"]]))
  if (max_dim > computed_dims) {
    max_dim <- computed_dims
    print(paste0("WARNING: max_dim too high, reduced to: ", computed_dims))
  }

  # Also ensure max_dim doesn't exceed the number of features
  n_features <- nrow(cells)
  if (max_dim >= n_features) {
    max_dim <- max(1, n_features - 1)
    print(paste0("WARNING: max_dim adjusted for feature count, set to: ", max_dim))
  }
} else {
  print("WARNING: PCA not available, setting max_dim to 0")
  max_dim <- 0
}

# Get features from the appropriate assay
adt_available <- ("ADT" %in% names(cells@assays) &&
                  !is.null(cells@assays$ADT@scale.data) &&
                  nrow(cells@assays$ADT@scale.data) > 0)
if (adt_available) {
  total_features <- cells@assays$ADT@scale.data %>% rownames()
} else {
  # Use the default assay if ADT is not available or doesn't have scaled data
  default_assay <- DefaultAssay(cells)
  scale_data_available <- (!is.null(cells@assays[[default_assay]]@scale.data) &&
                           nrow(cells@assays[[default_assay]]@scale.data) > 0)
  if (scale_data_available) {
    total_features <- cells@assays[[default_assay]]@scale.data %>% rownames()
  } else {
    total_features <- rownames(cells)
  }
}

if (pca_success && max_dim > 0) {
  tryCatch({
    cells <- Seurat::FindNeighbors(cells, reduction = "pca",
                                  dims = 1:max_dim) %>%
             FindClusters(resolution = 0.1, n.start = 10)
    cells <- BuildClusterTree(cells, reorder = TRUE, verbose = TRUE,
                             dims = 1:max_dim)
    png("cluster_tree.png")
    PlotClusterTree(object = cells)
    dev.off()  # Close the PNG device

  # Cluster tabulation
  cluster_table <- cells@meta.data$seurat_clusters %>% table()
  cluster_table <- tibble(cluster = names(cluster_table),
                         percentage = (cluster_table / sum(cluster_table)) *
                           100,
                         count = cluster_table)
  write.csv(cluster_table, "cluster_freq_table.csv")

  # find markers for every cluster compared to all remaining cells, report
  # only the positive ones
  markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25,
                           logfc.threshold = 0.25)
  markers %>% group_by(cluster) %>% top_n(n = 2)

  # Note: >900 genes was used in the training set followed by a test set to
  # deal with dropouts. Not sure if this is the best way to go about it.
  top <- markers %>% group_by(cluster) %>% top_n(n = 3)

  DoHeatmap(subset(cells, downsample = 300), features = total_features) +
    NoLegend() + theme(text = element_text(size = 6))
  ggsave("cluster_heatmap_total.png", dpi = dpi * 8)
  DoHeatmap(subset(cells, downsample = 300), features = top$gene) +
    NoLegend()
  ggsave("cluster_heatmap_top_genes.png", dpi = dpi * 8)

  # Per-cluster heatmaps
  # Step in to the per_cluster_heatmaps when making the heatmaps, then return
  # to the output directory when done
  dir.create("per_cluster_heatmaps", showWarnings = FALSE)
  setwd("per_cluster_heatmaps")
  clusters <- sort(as.numeric(levels(unique(markers$cluster))))
  for (i in clusters) {
    curr_markers <- markers[markers$cluster == i, ] %>%
      top_n(n = 100) %>%
      .[["gene"]] %>%
      .[complete.cases(.)]
    DoHeatmap(subset(cells, downsample = 300), features = curr_markers,
             size = 2) +
      NoLegend() +
      theme(text = element_text(size = 6))
    cluster_filename <- paste0("cluster.", i, ".png")
    ggsave(cluster_filename, width = 4, height = 2, dpi = dpi)
  }
  # Return to output directory before leaving try block
  setwd(out_dir)
  },
  error = function(e) {
    print(e)
    print("Clustering failed.")
    # Ensure we return to output directory even on error
    setwd(out_dir)
  })
} else {
  print("Skipping clustering - insufficient dimensions or PCA not available")
}

# Change back to out_dir even if there's a failure in clustering
setwd(out_dir)


# Dimr
tsne_successful <- FALSE
umap_successful <- FALSE

if (pca_success && max_dim > 0) {
  tsne_out <- tryCatch({
    cells <- Seurat::RunTSNE(object = cells, dims = 1:max_dim,
                            check_duplicates = FALSE)
    "success"  # Return success indicator instead of cells object
  },
  error = function(e) {
    print(e)
    print("tsne generation failed")
    "error"
  })
  tsne_successful <- (tsne_out == "success")

  umap_out <- tryCatch({
    cells <- Seurat::RunUMAP(object = cells, dims = 1:max_dim)
    "success"  # Return success indicator instead of cells object
  },
  error = function(e) {
    print(e)
    print("umap generation failed")
    "error"
  })
  umap_successful <- (umap_out == "success")
} else {
  print("Skipping TSNE and UMAP - PCA not available or insufficient dimensions")
}

# Dimr and heatmap viz
if (umap_successful) {
  DimPlot(cells, reduction = "umap", label = TRUE)
  ggsave("labeled_umap.png", height = 6, width = 6 * tsne_aspect_ratio,
         dpi = dpi * 6)
}
if (tsne_successful) {
  DimPlot(cells, reduction = "tsne", label = TRUE)
  ggsave("labeled_tsne.png", height = 6, width = 6 * tsne_aspect_ratio,
         dpi = dpi * 6)
}

# Dimr plotting for all markers
tryCatch({
  dir.create("per_gene_viz", showWarnings = FALSE)
  setwd("per_gene_viz")

  for (i in all_antibodies) {
    create_feature_plots(cells, i, tsne_successful, umap_successful, dpi)
  }

  # Return to output directory
  setwd(out_dir)
}, error = function(e) {
  print("Error in per-gene visualization")
  print(e)
  # Ensure we return to output directory even on error
  setwd(out_dir)
})

# Output (already in out_dir)
saveRDS(cells, "fully_processed.cells.rds")
out <- extract_scaled_data(cells)
write.csv(out, "fully_processed.cells.csv", row.names = TRUE)

# Write Seurat metadata table
meta <- cells@meta.data %>%
  data.table::as.data.table(., keep.rownames = "barcode")

# Go through each dimensionality reduction and merge them into single data.table
if (length(names(cells@reductions)) > 0) {
  tryCatch({
    reduction_data <- names(cells@reductions) %>%
      # Extract dim reduction matrices from seurat object
      lapply(., FUN = get_reductions, obj = cells) %>%
      # This goes through each table and merges them one by one with barcode as key
      purrr::reduce(., .f = function(x, y) x[y, on = "barcode"])

    # Now we have all dim reductions in single table, we can merge it to meta table
    meta <- meta[reduction_data, on = "barcode"]
  }, error = function(e) {
    print("Warning: Failed to merge dimensionality reduction data with metadata")
    print(e)
  })
}

seurat_meta <- paste0(opts$sample, "_seurat_metadata.csv.gz")
data.table::fwrite(meta, seurat_meta, quote = TRUE)

# Written at the end in case any of the settings get changed during the
# analysis (max_dim, for instance)
print("preserve settings")
curr_settings <- tibble(setting_names = "number_pc", setting_values = max_dim)
write.csv(curr_settings, "thisrun.settings.csv")
