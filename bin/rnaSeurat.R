#!/usr/bin/env Rscript
# Bio-Rad Laboratories, Inc.

# Purpose: Cluster cells using a feature count matrix

# Setting environment ------------------------------------------
# -------------------------------

# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
sampleId <- args[1]
dataDir <- args[2]
mixed <- args[3]


# Checking if inputs are set
if (is.na(sampleId)) {
  stop("Sample ID must be provided.")
} else if (!dir.exists(dataDir)) {
  stop(paste(dataDir, "does not exist."))
}

# Loading dependencies ----------------------------------------
# --------------------------------

packages <- c(
  "Seurat", "Matrix", "magrittr", "data.table", "SeuratWrappers", "stringr"
)
invisible(lapply(packages, library, character.only = TRUE))

data.table::setDTthreads()
future::plan("multicore", workers = parallel::detectCores())

set.seed(-1)

# Functions ----------------------------------------------------
# -------------------------------

# Function for finding the elbow of a scree plot
# using the second derivative of the curve
find_elbow <- function(stdev_vec) {
  # Creating coordinates from input vector
  x_vec <-
    seq_along(stdev_vec) # Creating vector of corresponding PCs
  # nolint start
  y_vec <- stdev_vec
  # nolint end
  # Fitting polynomial equation to the data to smooth the curve
  fit_eq <- lm(y_vec ~ poly(x_vec, 7))
  # Simulating data using fit model equation
  y_sim_vec <- predict.lm(fit_eq, list(x = x_vec))
  # Normalizing data so both data sets are constrained to between 0 and 1
  x_norm_vec <- (x_vec - min(x_vec)) / (max(x_vec) - min(x_vec))
  y_norm_vec <-
    (y_sim_vec - min(y_sim_vec)) / (max(y_sim_vec) - min(y_sim_vec))
  # Calculating the difference curve coordinates
  y_diff_vec <- max(y_norm_vec) - y_norm_vec - x_norm_vec
  # Finding local maximas using second derivative
  y_maxima_idx <- which(diff(sign(diff(y_diff_vec))) == -2) + 1
  # Retreiving values from local maximas
  y_maxima_vec <- y_diff_vec[y_maxima_idx]
  # If maxima were identified
  if (length(y_maxima_vec) != 0) {
    # Identifying the local maxima(s) with the highest magnitude
    elbow_int <-
      as.integer(names(which(y_maxima_vec == max(y_maxima_vec))))
    # If more than elbow has the largest inflection point, use the first one
    if (length(elbow_int) > 1) {
      elbow_int <- elbow_int[1]
    }
    # If no maxima were identified
  } else {
    elbow_int <- length(stdev_vec)
    warning("No elbow found. Using all available components.")
  }
  # Outputting result
  return(elbow_int)
}

# Function for reading in Bio-Rad Laboratories, Inc. formatted sparse matrices
# Inputs are expected as a vector of paths
# Assumes one sample per data.dir
# Assumes files are named as ${sampleId}.(filtered).barcodes.tsv(.gz),
# ${sampleId}.(filtered).features.tsv(.gz), ${sampleId}.(filtered).mtx.gz
readBioRad <-
  function(data.dirs = NULL,
           feature.column = 1,
           unique.features = TRUE,
           append.sampleId = TRUE,
           mixed) {
    cat("Data directory is:", data.dirs, "\n")
    data.dirs <- as.list(data.dirs)
    obj <-
      lapply(
        data.dirs,
        readBioRadSample,
        mixed,
        feature.column,
        append.sampleId
      )
    obj <- do.call(cbind, obj)
    return(obj)
  }

readBioRadSample <-
  function(data.dir,
           mixed,
           feature.column = 1,
           append.sampleId = TRUE) {
    mtxFile <- list.files(data.dir, pattern = "*.mtx(.gz)?")
    barcodeFile <- list.files(data.dir, pattern = "*.barcodes.tsv(.gz)?")
    featureFile <- list.files(data.dir, pattern = "*.genes.tsv(.gz)?")

    if (any(c(length(barcodeFile), length(featureFile), length(mtxFile)) > 1)) {
      stop("Input directory contains more than one barcode, feature, or matrix file.")
    }

    sampleName <- gsub("(.filtered)?.mtx.gz", "", mtxFile)
    cat("matrix file name", mtxFile, "\n")
    mtx <- Matrix::readMM(mtxFile)
    # Read in barcodes ands features
    barcodes <- data.table::fread(barcodeFile, header = FALSE)
    if (append.sampleId) {
      barcodes[[1]] <- paste(barcodes[[1]], sampleName, sep = "-")
    }

    features <- data.table::fread(featureFile, header = FALSE)
    if (!all(dim(mtx) == c(nrow(features), nrow(barcodes)))) {
      stop("Count matrix dimensions not equal to length of features or barcodes.")
    }
    if (feature.column > ncol(features)) {
      write(
        "[WARNING]: feature.column requested is greater than
         columns in features.tsv, defaulting to first column.",
        stderr()
      )
      feature.column <- 1
    }

    if (mixed) {
      # appending species name to gene name in mixed species experiments
      # Need to assume column 1 = ensembl id, column 2 is the gene symbol
      if (feature.column == 2) {
        features[, species := as.character(NA)] # initialize column
        features[str_detect(features[[1]],
          pattern = "^ENSG[0-9]"), species := "Homo-sapiens"]
        features[str_detect(features[[1]],
          pattern = "^ENSMUSG[0-9]"), species := "Mus-musculus"]
        features[, gene_species := paste0(features[[feature.column]], "-", species)]
      } else {
        cat("No Ensembl ID associated with the genes, proceeding without labeling.")
      }

      actualSpecies <- unique(features$species)
      expectedSpecies <- c("Homo-sapiens", "Mus-musculus")
      if (length(actualSpecies) > 2) {
        cat("More than 2 species detected. This is not supported.")
      } else if (!setequal(actualSpecies, expectedSpecies)) {
        cat("Species detected are not in the supported list of species.")
      }

      rownames(mtx) <- features$gene_species
    } else {
      rownames(mtx) <- features[[feature.column]]
    }

    colnames(mtx) <- barcodes[[1]]
    if (class(mtx) == "ngTMatrix") {
      mtx <- mtx * 1
    }
    mtx <- as(mtx, Class = "dgCMatrix")
    return(mtx)
  }

# This extracts dimensionality reduction matrices from Seurat object
get_reductions <- function(reduction, obj) {
  obj@reductions[[reduction]]@cell.embeddings[, 1:2] %>%
    data.table::as.data.table(., keep.rownames = "barcode", key = "barcode")
}

# Analysis ---------------------------------------------------
# ---------------------------------

# Importing and transforming data
# Reading in sparse matrix using matrix.mtx, features.tsv,
#  and barcodes.tsv files from input dir
matrix <- readBioRad(data.dir = dataDir, feature.column = 2, mixed = mixed)

# Converting matrix to Seurat object
obj <- Seurat::CreateSeuratObject(
  counts = matrix,
  project = sampleId,
  min.cells = 1,
  min.features = 1
)


tryCatch({
  # Normalizing data, scaling data, and selecting features
  # This function will take the count slot and log1p transform
  # the log1p, (log(1+count), transformation writes into data slot
  obj <- Seurat::SCTransform(obj, verbose = FALSE, seed.use = 1448145)
  # Setting number of components to use for PCA
  if (ncol(matrix) >= 50) {
    npcs <- 50
  } else {
    npcs <- ncol(matrix)
  }

  # Running PCA and outputting ", npcs, " components
  # Dimensionality reduction with PCA and outputting the specified number of PCs
  obj <-
    Seurat::VariableFeatures(object = obj) %>%
    Seurat::RunPCA(obj, features = ., npcs = npcs, seed.use = 42)

  # Identifying inflection point of PCA scree plot
  pca <- obj[["pca"]]
  elbow <- find_elbow(pca@stdev)

  # Setting the optimum number of PCs to use downstream
  if (elbow < npcs) {
    opt_npcs <- elbow + 1
    cat(
      paste0(
        "[PROGRESS]: PCA scree plot inflection at PC",
        elbow,
        ". Using ",
        opt_npcs,
        " components for downstream analyses.\n"
      )
    )
  } else {
    opt_npcs <- elbow
    cat(
      paste0(
        "[PROGRESS]: No inflection point identified in PCA scree plot. Using ",
        opt_npcs,
        " components for downstream analyses.\n"
      )
    )
  }

  # Running dimensionality reduction with UMAP
  # Dimensionality reduction with UMAP using optimum number of PCs
  obj <- Seurat::RunUMAP(obj, dims = 1:opt_npcs, verbose = FALSE, seed.use = 42)

  # Constructing K-nearest neighbors graph
  # Using K-nearest neighbors (KNN) to calculate inter-cell distances
  obj <- Seurat::FindNeighbors(obj, dims = 1:opt_npcs, verbose = FALSE)

  # Identifying clusters
  # Identifying clusters using KNN graph
  obj <- Seurat::FindClusters(obj, verbose = FALSE, random.seed = 0)

  # Saving clustered object
  # Saving analysis
  saveRDS(obj, file = paste0(sampleId, "_seurat.rds"))

  # Saving metadata
  # Write Seurat metadata table
  meta <- obj@meta.data %>%
    data.table::as.data.table(., keep.rownames = "barcode")

  # Go through each dimensionality reduction and merge them into single data.table
  meta <- names(obj@reductions) %>%
    # Extract dim reduction matricies from seurat object
    lapply(., FUN = get_reductions, obj = obj) %>%
    # This goes through each table and merges them one by one with barcode as key
    purrr::reduce(., .f = function(x, y) x[y, on = "barcode"]) %>%
    # Now we have all dim reductions in single table, we can merge it to meta table
    meta[., on = "barcode"]

  seurat_meta <- paste0(sampleId, "_seurat_metadata.csv.gz")
  data.table::fwrite(meta, seurat_meta, quote = TRUE)

  seurat_simple <- meta[, c("UMAP_1", "UMAP_2", "seurat_clusters", "nCount_RNA")]
  data.table::fwrite(seurat_simple, paste0(sampleId, "_umap.csv"))

  # Differences are calculated from the data slot
  top_features <- SeuratWrappers::RunPrestoAll(obj,
    min.pct = 0.1, logfc.threshold = 0, min.cells.group = 0,
    return.thresh = 1, base = 2, random.seed = 1
  )

  # Converts to data.table
  data.table::setDT(top_features)

  # Fixed cutoffs for top features
  top_feature_count <- 250
  p_value_cutoff <- 0.01

  # Sort by avg_log2FC descending
  data.table::setorder(top_features, -avg_log2FC)
  # Filter table by p value
  feature_filter <- top_features[p_val_adj < p_value_cutoff, ]
  # Select top features for each cluster
  feature_filter <- feature_filter[, head(.SD, top_feature_count), by = cluster]
  # Keep unique features
  feature_filter <- unique(feature_filter[, c("gene")])
  # Drop features not in the filtered list
  top_features <- top_features[feature_filter, on = "gene", nomatch = 0L]
  # put the p_val_adj in scientific notation and to use 3 sig figs
  top_features$p_val_adj <- formatC(signif(top_features$p_val_adj, digits = 3),
    digits = 3, format = "e"
  )
  # make the avg_log2FC only use 3 sig figs
  top_features$avg_log2FC <- signif(top_features$avg_log2FC, digits = 3)
  # Reshape data from to wide format
  top_features_wide <- data.table::dcast(top_features,
    gene ~ cluster,
    value.var = c("avg_log2FC", "p_val_adj"),
    drop = FALSE
  )

  # This will group and order padj & log2FC by cluster number
  colnames(top_features_wide) %>%
    stringr::str_extract(., "^gene|[0-9]*$") %>%
    as.numeric(.) %>%
    # The NA is for the gene column since its not numeric,
    #  this will keep it as 1st col
    order(., na.last = FALSE) %>%
    # Updates order by reference
    data.table::setcolorder(top_features_wide, .)
  # Sort alphabetically by features
  data.table::setorder(top_features_wide, gene)

  # Saving the top features
  data.table::fwrite(top_features_wide,
    paste0(sampleId, "_top_features.csv"),
    na = "NA"
  )

}, error = function(e) {
  # Handle the error
  message("An error occurred: ", e$message)

  message_file <- paste(sampleId, "SEURAT_messages.txt", sep = "_")

  message_prefix <- paste0("Warning: [3' RNA Droplet] Sample ", sampleId)

  writeLines(paste0(
    message_prefix,
    " failed during Seurat analysis. ",
    "Placeholder data used for plots & table under results tab. ",
    "Please review user guide for more information."
  ), message_file)

  # Export empty objects for report
  saveRDS(obj, file = paste0(sampleId, "_seurat.rds"))

  meta <- data.table(
    barcode = c("ERROR"),
    orig.ident = c("ERROR"),
    nCount_RNA = c(0),
    nFeature_RNA = c(0),
    nCount_SCT = c(0),
    nFeature_SCT = c(0),
    SCT_snn_res.0.8 = c(0),
    seurat_clusters = c(0),
    PC_1 = c(0),
    PC_2 = c(0),
    UMAP_1 = c(0),
    UMAP_2 = c(0)
  )
  seurat_meta <- paste0(sampleId, "_seurat_metadata.csv.gz")
  data.table::fwrite(meta, seurat_meta, quote = TRUE)

  seurat_simple <- data.table(
    UMAP_1 = c(0),
    UMAP_2 = c(0),
    seurat_clusters = c(0),
    nCount_RNA = c(0)
  )
  data.table::fwrite(seurat_simple, paste0(sampleId, "_umap.csv"))


  top_features_wide <- data.table::data.table(
    gene = "ERROR", avg_log2FC_0 = 0, p_val_adj_0 = 0
  )
  data.table::fwrite(
    top_features_wide, paste0(sampleId, "_top_features.csv"),
    na = "NA"
  )

  # Exit the script
  quit("no", status = 0, runLast = FALSE)
})
