#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument("countsFile", nargs = 1, help = "Barcode counts file.")
parser$add_argument(
  "-o",
  "--output_directory",
  type = "character",
  default = getwd(),
  help = "Directory to write output files. [Default: Current working directory]"
)
parser$add_argument("-p",
                    "--prefix",
                    type = "character",
                    default = "alignments",
                    help = "Prefix for filenames. [Default: alignments]")

parser <- parser$parse_args()

# Set data.table threading
setDTthreads()

countsFile <- parser$countsFile
outDir <- parser$output_directory
if (!dir.exists(outDir)) {
  dir.create(outDir)
}
prefix <- parser$prefix

#read counts
bcCount <- fread(countsFile)
bcCount <- bcCount[order(uniqueNuclearFrags, decreasing = TRUE)]

#write beads per droplet into a column
getDropSize <- function(x) {
  spl <- unlist(strsplit(x, split = "_"))
  if (length(spl) > 1) {
    dropSize <- as.numeric(gsub("N", "", spl[length(spl)]))
  } else{
    dropSize <- 1
  }
  return(dropSize)
}

bcCount$beadsInDrop <- sapply(bcCount$DropBarcode, getDropSize)

fwrite(bcCount, file = paste0(outDir, "/", prefix, ".cell_data.csv"))
