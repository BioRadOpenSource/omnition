#!/usr/bin/env Rscript
# Caleb Lareau
# Bio-Rad Laboratories, Inc.

options(warn = -1)

suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(tools)))
suppressMessages(suppressWarnings(library(tidyr)))

"%ni%" <- Negate("%in%")

options(scipen = 999)

#-----------------
# Parse arguments
#-----------------

args <- commandArgs(trailingOnly = TRUE)
if (file_path_sans_ext(basename(args[1])) == "R") {
  i <- 2
} else {
  # Rscript
  i <- 0
}

tempQcIn <- args[i + 1]
mitoIn <- args[i + 2]
firstQCout <- args[i + 3]
species_mix <- args[i + 4]

# For devel only
if (FALSE) {
  #nolint start
  tempQcIn <- "/Users/clareau/dat/Research/BuenrostroResearch/lareau_dev/deconvolution/tests/deconvolution2/temp/jaccardPairsForIGV.basicQC-temp.tsv"
  mitoIn <- "/Users/clareau/dat/Research/BuenrostroResearch/lareau_dev/deconvolution/tests/deconvolution2/temp/filt_split/jaccardPairsForIGV.chrM_frag.sumstats.tsv"
  firstQCout <- "/Users/clareau/dat/Research/BuenrostroResearch/lareau_dev/deconvolution/tests/deconvolution2/final/jaccardPairsForIGV.basicQC.tsv"
  species_mix <- "none"
  #nolint end
}

nuc_in <-
  fread(tempQcIn,
        col.names = c("cell_barcode", "n_total", "n_unique", "chromosome"))
mito_in <-
  fread(mitoIn,
        col.names = c("cell_barcode", "n_total", "n_unique", "chromosome"))

# Summarize counts
# No species mixing
nuc_sum <- nuc_in %>%
 group_by(cell_barcode) %>%
  summarize(totalNuclearFrags = sum(n_total),
            uniqueNuclearFrags = sum(n_unique))

mito_sum <- mito_in %>%
 group_by(cell_barcode) %>%
  summarize(totalMitoFrags = sum(n_total),
            uniqueMitoFrags = sum(n_unique))
basic_qc <- left_join(nuc_sum, mito_sum, by = "cell_barcode")

if (species_mix == "true") {
  species_names <-
    nuc_in %>%
     select(chromosome) %>%
     extract(chromosome, into = c("species", "chr"),
     "^([^.]+)\\.(.*)") %>%
      summarise(names = sort(unique(species)))
  nuc_in$s1 <-
    as.numeric(grepl(species_names$names[1], nuc_in$chromosome))
  nuc_in$s2 <-
    as.numeric(grepl(species_names$names[2], nuc_in$chromosome))
  s1_df <- nuc_in %>%
   filter(s1 == 1) %>%
    group_by(cell_barcode) %>%
    summarize(totalS1Frags = sum(n_total),
              uniqueS1Frags = sum(n_unique))
  names(s1_df) <-
    c(
      "cell_barcode",
      paste0("total", species_names$name[1], "Frags"),
      paste0("unique", species_names$name[1], "Frags")
    )
  s2_df <- nuc_in %>%
  filter(s2 == 1) %>%
    group_by(cell_barcode) %>%
    summarize(totalS2Frags = sum(n_total),
              uniqueS2Frags = sum(n_unique))
  names(s2_df) <-
    c(
      "cell_barcode",
      paste0("total", species_names$name[2], "Frags"),
      paste0("unique", species_names$name[2], "Frags")
    )
  QCstats <- left_join(basic_qc, s1_df, by = "cell_barcode") %>%
    left_join(s2_df, by = "cell_barcode")
} else {
  QCstats <- basic_qc
}

QCstats[is.na(QCstats)] <- 0
write.table(
  QCstats,
  file = firstQCout,
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)
