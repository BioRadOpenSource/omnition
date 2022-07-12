#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))

#argparse
parser <- ArgumentParser()
parser$add_argument("countsFile", nargs = 1, help = "Barcode counts file.")
parser$add_argument("species1", nargs = 1, help = "Species ID 1")
parser$add_argument("species2", nargs = 1, help = "Species ID 2")
parser$add_argument("-y",
                    "--assay",
                    type = "character",
                    default = "atac",
                    help = "Assay used.")
parser$add_argument(
  "-a",
  "--allowlist",
  type = "character",
  default = "NULL",
  help = "Barcode allowlist file used if filtering barcodes."
)
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
parser$add_argument(
  "-ct",
  "--crosstalk_threshold",
  type = "double",
  default = 0.9,
  #nolint start
  help = "Minimum fraction of UMIs attributed to a single species for a cell to be called as that species. [Default: 0.9]"
  #nolint end
)
parser$add_argument(
  "-cb",
  "--combinatorial_barcode",
  type = "character",
  default = "MULTIPLEXED",
  help = "Combinatorial barcode."
)
parser <- parser$parse_args()

countsFile <- parser$countsFile
allowlist <- parser$allowlist
species1 <- parser$species1
species2 <- parser$species2
prefix <- parser$prefix
outDir <- parser$output_directory
crosstalkThreshold <- parser$crosstalk_threshold
assay <- parser$assay
barcode <- parser$combinatorial_barcode

#input above knee barcodes
counts <- fread(countsFile)
featS1 <- "Percent_%s_fragments"
featS2 <- "unique"
featS3 <- "Frags"

#count things
if (!allowlist == "NULL") {
  validBarcodes <- fread(allowlist, header = F)
  counts <- counts %>% filter(barcode %in% validBarcodes$V1)
}

if (!barcode == "MULTIPLEXED") {
  counts <- counts %>% filter(grepl(barcode, DropBarcode))
}

counts <- counts %>%
  mutate(!!sprintf(featS1, species1) := !!as.name(paste0(featS2, species1, featS3)) /
           uniqueNuclearFrags) %>%
  mutate(!!sprintf(featS1, species2) := !!as.name(paste0(featS2, species2, featS3)) /
           uniqueNuclearFrags)

counts$cellType <-
  cut(
    counts[[sprintf(featS1, species1)]],
    breaks = c(
      -Inf,
      as.numeric(1 - crosstalkThreshold),
      as.numeric(crosstalkThreshold),
      Inf
    ),
    labels = c(
      species2,
      "mixed",
      species1
    )
  )

s1cells <- sum(counts$cellType == species1)
s2cells <- sum(counts$cellType == species2)
mixedCells <- sum(counts$cellType == "mixed")
if (any(s1cells == 0, s2cells == 0)) {
  cellPurity <- 100
} else{
  s1ReadsInS1Cells <-
    sum(counts %>%
          filter(cellType == species1) %>%
          select(!!as.name(paste0(
            featS2, species1, featS3
          ))))
  s2ReadsInS2Cells <-
    sum(counts %>%
          filter(cellType == species2) %>%
          select(!!as.name(paste0(
            featS2, species2, featS3
          ))))
  s1ReadsInS2Cells <-
    sum(counts %>%
          filter(cellType == species2) %>%
          select(!!as.name(paste0(
            featS2, species1, featS3
          ))))
  s2ReadsInS1Cells <-
    sum(counts %>%
          filter(cellType == species1) %>%
          select(!!as.name(paste0(
            featS2, species2, featS3
          ))))
  cellPurity <-
    signif((s1ReadsInS1Cells + s2ReadsInS2Cells) / (
      s1ReadsInS1Cells + s2ReadsInS2Cells + s1ReadsInS2Cells + s2ReadsInS1Cells
    ) * 100,
    3
    )
}

crosstalk <-
  signif(mixedCells / (s1cells + s2cells + mixedCells) * 100,
         3)
crosstalk_total <-
  signif(2 * mixedCells / (s1cells + s2cells + mixedCells) * 100,
         3)

df <- data.frame(
  s1cells,
  s2cells,
  Mixed_cells = mixedCells,
  Measurable_crosstalk = paste(crosstalk, "%"),
  Estimated_total_crosstalk = paste(crosstalk_total, "%"),
  Cell_purity = cellPurity
)

colnames(df)[c(1:2)] <-
  c(
    paste0(species1, "_cells"),
    paste0(species2, "_cells")
  )

# write rds
saveRDS(
  list(
    counts = counts,
    summaryTable = df,
    species1 = species1,
    species2 = species2
  ),
  file = paste0(outDir, "/", prefix, ".species_mix_counts.rds")
)

write.table(
  df,
  file = paste0(outDir, "/", prefix, ".crosstalk.csv"),
  sep = ",",
  row.names = F,
  quote = F
)

fwrite(counts, file = paste0(outDir, "/", prefix, ".species_mix_counts.csv"))
