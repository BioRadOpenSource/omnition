#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(devtools))

parser <- ArgumentParser()
parser$add_argument("-g", "--gtf", help = "Path to gtf file.")
parser <- parser$parse_args()
gtf <- parser$gtf

# Importing the gtf into a granges object to get a list of chromosomes in the gtf file
granges <- import(gtf)

## Getting a list of all the chromosomes present in the GTF file
chr_gtf <- levels(seqnames(granges))

## Adding fa so they can be compared with the chromosome
# files made by splitting up the reference fasta
chr_gtf_fa <- paste0(chr_gtf, ".fa")

## Getting all the fasta files created by breaking up the reference fasta
chrs_fa <- list.files(path = "./refchrs", pattern = "\\.fa")

## Creating a list of only the chromosomes present in both the fasta and
# gtf to avoid future conflict between the two
chrs <- intersect(chr_gtf_fa, chrs_fa)

# Start of the line for setting the fasta sequences to use in the seed file
str_chrs <- "seqnames: c("

for (c in chrs) {
  # pasting each file without .fa into one long string to be used in the seed file
  str_chrs <- paste0(str_chrs, "\"", substr(c, 0, nchar(c) - 3), "\",")
}

# Removing the last , and adding a closing ) to the seqnames line
str_chrs <- paste0(substr(str_chrs, 0, nchar(str_chrs) - 1), ")")

## Writing the contents of the seed file that will be used to create the
# BSgenome package
seedfile <- file("BSgenome.ref-seed")
writeLines(
  c(
    "Package: BSgenome.ref.na.1.0",
    "Title: Reference genome created from ref.fa",
    "Description: Reference genome created from ref.fa.",
    "Version: 1.0",
    "Author: na",
    "Maintainer: na <Jane.Developer@some.domain.net>",
    "License: GPL (>= 2)",
    "organism: na",
    "common_name: na",
    "provider: na",
    "release_date: na",
    "source_url: na",
    "organism_biocview: na",
    "genome: ref",
    "BSgenomeObjname: ref",
    "circ_seqs: character(0)",
    str_chrs,
    "seqs_srcdir: ./refchrs",
    "seqfiles_suffix: .fa",
    "seqfiles_prefix:"
  ),
  seedfile
)
close(seedfile)

## Making the BSgenome pacakge for use as a reference for ArchR with the seed file
forgeBSgenomeDataPkg("BSgenome.ref-seed")

## Build BSgenome the package
build("BSgenome.ref.na.1.0")
