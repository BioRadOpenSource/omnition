#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(ArchR))
set.seed(1)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))


## Load the BSgenome reference package that was made
suppressPackageStartupMessages(library(BSgenome.ref.na.1.0))

parser <- ArgumentParser()
parser$add_argument("-g", "--gtf", help = "Path to gtf file.")
parser$add_argument("-b", "--blocklist", help = "Path to blocklist bed file.")
parser$add_argument("-p", "--peaks", help = "Path to peaks bed file.")
parser$add_argument("-c", "--cpus", default = 1, type = "integer",
 help = "Number of CPUs to use in multiprocessing.")
parser$add_argument("-t", "--tile_size", default = 500,
 help = "Bin size for genome tiling.")
parser$add_argument("-d", "--demo", help = "Boolean to indicate if data is demo data.")
parser <- parser$parse_args()

gtf <- parser$gtf
blocklist <- parser$blocklist
cpus <- parser$cpus
peaks <- parser$peaks
tile_size <- parser$tile_size
demo <- parser$demo

demo <- demo == "True"

## Setting default number of Parallel threads to 16.
addArchRThreads(threads = as.numeric(cpus))
## Setting Parallel threads to 1 can fix some errors

## removing the requirement for a specific prefix for chromosome names
addArchRChrPrefix(chrPrefix = FALSE)

## Importing the blocklist bed file as a granges object
gblocklist <- import(blocklist)

## Using the reference BSgenome package and blocklist to set the reference genome
genomeAnnotation <-
  createGenomeAnnotation(
    genome = BSgenome.ref.na.1.0,
    filter = FALSE,
    filterChr = NULL,
    blacklist = gblocklist
  )

## Importing the gtf file into a granges object
granges <- import(gtf)

## Adding a column named symbol to the granges object
# which is what ArchR uses to relate between TSS, exons, and genes
granges$symbol <- mcols(granges)[, c("gene_id")]

## Setting geneAnnotation for the reference using the gtf files
geneAnnotation <- createGeneAnnotation(TSS =
 resize(granges[mcols(granges)[, c("type")] ==
 "transcript"], width = 1),
  exons = granges[mcols(granges)[, c("type")] == "exon"],
  genes = granges[mcols(granges)[, c("type")] == "transcript"])

## Setting path to samples and getting sample names
sample <- list.files(pattern = "\\.final.bam$")
samplename <-
  gsub(".final.bam", "", list.files(pattern = "\\.final.bam$"))

## Demo data specific settings
if (demo) {
  addGeneScoreMat <- TRUE

  geneScoreMatList <- list(useGeneBoundaries = FALSE,
    extendUpstream = c(1, 5e+05),
    extendDownstream = c(1, 5e+05))
} else {
  addGeneScoreMat <- FALSE

  geneScoreMatList <- list(useGeneBoundaries = TRUE)
}

## Creating arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = sample,
  sampleNames = samplename,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  minTSS = 0,
  minFrags = 1,
  addTileMat = TRUE,
  TileMatParams = list(tileSize = tile_size),
  addGeneScoreMat = addGeneScoreMat,
  GeneScoreMatParams = geneScoreMatList,
  excludeChr = c("None"),
  bcTag = "XC",
  # threads = 1,
  subThreading = FALSE,
  bamFlag = list(
    isMinusStrand = FALSE,
    isProperPair = TRUE,
    isDuplicate = FALSE
  ),
  force = TRUE
)

## Creating archr project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "ArchR",
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  copyArrows = TRUE #This is recommened so that you maintain an unaltered
  # copy for later usage. Can be removed to lower memory usage
)

## Set parameters for umap
LSImatrix <- "TileMatrix"
LSIselection <- "top"
LSIbinarize <- TRUE
LSIvarFeatures <- 25000
LSIdimsToUse <- c(1:30)
LSIiterations <- 2
if (length(proj$cellNames) < 40) {
  UMAPneighbors <- length(proj$cellNames)
} else {
  UMAPneighbors <- 40
}
UMAPverbosity <- TRUE

## Dimensionality reduction
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = LSImatrix,
  name = "IterativeLSI",
  firstSelection = LSIselection,
  binarize = LSIbinarize,
  varFeatures = LSIvarFeatures,
  dimsToUse = LSIdimsToUse,
  iterations = LSIiterations
)

## Clustering
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")

## UMAP embedding
proj <-
  addUMAP(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    nNeighbors = UMAPneighbors,
    verbose = UMAPverbosity
  )

## Color the UMAP based on a specific variable for each cluster
p1 <-
  plotEmbedding(
    ArchRProj = proj,
    colorBy = "cellColData",
    name = "Clusters",
    embedding = "UMAP"
  )

plotPDF(
  p1,
  name = paste0("Plot-", samplename, "-UMAP-Clusters.pdf"),
  ArchRProj = proj,
  addDOC = FALSE,
  width = 5,
  height = 5
)

## Creating and saving an R data frame that has umap
# coordinates and cluster ID for each barcode
umap_and_cluster <- getEmbedding(proj)
umap_and_cluster$Clusters <-
  proj$Clusters[which(!is.na(proj$Clusters))]
save(umap_and_cluster,
     file = paste0("ArchR/", samplename, "_umap_and_clusterID.rda"))

## Differential Chromatin Accessibility wrapped in tryCatch
tryCatch(
  expr = {
    ### Import peaks as GenomicRanges
    peaks_gr <-
      makeGRangesFromDataFrame(read_tsv(peaks, col_names =
       c("seqname", "start", "end")))
    ### Add peakset to ArchR Project
    proj <- addPeakSet(ArchRProj = proj, peakSet = peaks_gr)
    ### Add a Peak Matrix
    proj <- addPeakMatrix(ArchRProj = proj)
    ### Create marker heatmap
    #### Get marker features
    markers_peaks <- getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "PeakMatrix",
      groupBy = "Clusters",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon"
    )
    #### As pdf
    marker_heatmap <- plotMarkerHeatmap(
      seMarker = markers_peaks,
      plotLog2FC = ifelse(ncol(markers_peaks) <= 2, TRUE, FALSE),
      cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
      transpose = FALSE
    )
    #### As data frame
    marker_heatmap_df <- plotMarkerHeatmap(
      seMarker = markers_peaks,
      plotLog2FC = ifelse(ncol(markers_peaks) <= 2, TRUE, FALSE),
      cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
      transpose = FALSE,
      returnMatrix = TRUE
    )
    ### Export PDF
    plotPDF(
      marker_heatmap,
      name = paste0("Plot-", samplename, "-peak-heatmap.pdf"),
      ArchRProj = proj,
      addDOC = FALSE,
      width = 8,
      height = 6
    )
    ### Export data frame
    write_tsv(
      rownames_to_column(as.data.frame(marker_heatmap_df), var = "peak"),
      file = paste0("ArchR/", samplename, "-peak-heatmap.tsv")
    )
  },
  error = function(e) {
    cat("[WARNING]: Unable to detect differential chromatin accessibility.\n")
    print(e)
  }
)


## Function to generate sparse matrices
getSparseMatrix <- function(arrow_file = ArrowFile, tile_size = 500) {
    mat <- ArchR:::.getMatFromArrow(ArrowFile = arrow_file)
    # Get rownames
    rows <- h5read(arrow_file, "/TileMatrix")$Info$FeatureDF
    rows <- paste0(rows$seqnames, ":", rows$start, "-", rows$start + tile_size - 1)
    # Add dimnames
    rownames(mat) <- rows
    dimnames(mat) <- list(rows, colnames(mat))
    return(mat)
}

## Get sparse matrices from arrow files and write to disk
tryCatch(
  expr = {
    # get arrow file name
    ArrowFile <- list.files("ArchR/ArrowFiles/", pattern = "*.arrow")
    # create mtx
    mtx <- getSparseMatrix(arrow_file = ArrowFile, tile_size = tile_size)
    # write sparse matrix to disk
    writeMM(mtx, file = paste0("ArchR/", samplename, ".mtx"))
    write.table(rownames(mtx), file = paste0("ArchR/", samplename, ".rownames.tsv"),
     sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(colnames(mtx), file = paste0("ArchR/", samplename,
     ".colnames.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  },
  error = function(e) {
    cat("[WARNING]: Unable to write tile matrix to disk.\n")
    print(e)
  }
)

## Saving archr project
proj <- saveArchRProject(ArchRProj = proj)
