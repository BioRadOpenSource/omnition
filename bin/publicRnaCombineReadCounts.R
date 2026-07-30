#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
options(warn = -1)

parser <- ArgumentParser()

parser$add_argument("genome", nargs = 1, help = "Genome read counts file.")
parser$add_argument("unaligned", nargs = 1, help = "Unaligned read counts file.")
parser$add_argument("-p", "--prefix",
                    default = "mySample",
                    help = "Sample name to prepend onto outputs. [Default: mySample]")

parser <- parser$parse_args()

genomeFile <- parser$genome
unalFile <- parser$unaligned
prefix <- parser$prefix

genome <- fread(genomeFile)
unal <- fread(unalFile)

# Set data.table threading
setDTthreads()

#append dup vs non onto alignment columns
colnames(genome)[2:ncol(genome)] <- paste0(colnames(genome)[2:ncol(genome)], ".raw")

#create data.tables
genomeCount <- data.table(barcode = genome$barcode,
                          alignedReads = genome[, rowSums(.SD, na.rm = TRUE),
                          .SDcols = grep("readCount", names(genome))]
                          )
unalCount <- data.table(barcode = unal$barcode,
                        unalignedReads = unal[, rowSums(.SD, na.rm = TRUE),
                        .SDcols = grep("readCount", names(unal))]
                        )

rm(list = c("genome", "unal"))

#set keys
setkey(genomeCount, barcode)
setkey(unalCount, barcode)

#merge then sort
cmb <- Reduce(function(...) merge(..., all = TRUE), list(genomeCount, unalCount))
cmb <- setorder(cmb, -"alignedReads", na.last = TRUE)

#and write
fwrite(cmb, file = paste0(prefix, ".counts_per_barcode.csv"), sep = ",")
