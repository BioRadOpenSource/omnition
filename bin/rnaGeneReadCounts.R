#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(magrittr))

parser <- ArgumentParser()

parser$add_argument("samplesDir", nargs = 1,
                    help = "Path to directory containing mtx formatted count matrices.")
parser$add_argument(
  "-o",
  "--output_directory",
  type = "character",
  default = getwd(),
  help = "Directory to write output files. [Default: Current working directory]"
)

parser <- parser$parse_args()

samplesDir <- parser$samplesDir

#every .mtx.gz is assumed to be a counts matrix with
# an associated genes.tsv and barcodes.tsv
samples <-
  gsub(".mtx.gz", "", list.files(samplesDir, pattern = ".mtx.gz"))

for (sampleId in samples) {
  genes <-
    read.table(
      paste0(samplesDir, "/", sampleId, ".genes.tsv"),
      col.names = c("gene_id", "gene_name")
    )

  # get the ensembl stable ID prefix for the gene
  # remove the regex that was searched for to isolate the species prefix
  genes %<>% mutate(species = str_extract(gene_id, regex(".*(?<=G[0-9])")))
  genes$species <- gsub("G[0-9]", genes$species, replacement = "")
  counts <- readMM(paste0(samplesDir, "/", sampleId, ".mtx.gz"))
  genes$count <- rowSums(counts)
  names(genes)[names(genes) == "count"] <- sampleId
  if (!exists("counts_table")) {
    counts_table <- genes
  } else {
    counts_table <- full_join(counts_table, genes)
  }
}

counts_table[is.na(counts_table)] <- 0
species <- unique(counts_table$species)
# check if mixed species, remove the species column from the output .tsv
if (length(species) > 1) {
  for (i in seq_len(length(species))) {
    write.table(
      split(counts_table, counts_table$species)[[i]] %>% select(-species),
      file = paste("species",
       as.character(i),
        "_gene_umi_counts_per_sample.tsv", sep = ""),
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
      sep = "\t"
    )
  }
} else {
  write.table(
  counts_table %>% select(-species),
  file = "gene_umi_counts_per_sample.tsv",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)
}
