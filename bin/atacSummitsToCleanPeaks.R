#!/usr/bin/env Rscript

# input: macs2 summit calls, blocklist regions, integer
# for bp padding, chromosome sizes and top n
# output: top n refined peak calls with fixed width, no blocklist

suppressPackageStartupMessages(library(argparse))

"%ni%" <- Negate("%in%")

parser <- ArgumentParser()

parser$add_argument("summit_file", nargs = 1, help = "Summits from MACS2.")
parser$add_argument("chrom_sizes_file", nargs = 1, help = "Chromosome sizes file.")
parser$add_argument("-b", "--blocklist", help = "Blocklist .bed file for genome used.")
parser$add_argument(
  "-w",
  "--peak_width",
  type = "double",
  default = 250,
  help = "Fixed width in base pairs for output peaks. [Default: 250]"
)
parser$add_argument("-n",
                    "--n",
                    type = "double",
                    default = 999999999,
                    help = "Unknown n argument. [Default: 999999999]")
parser$add_argument(
  "-f",
  "--fdr_threshold",
  type = "double",
  default = 0.01,
  help = "False discovery rate threshold. [Default: 0.01]"
)
parser$add_argument(
  "-o",
  "--output_directory",
  type = "character",
  default = getwd(),
  help = "Directory to write output files. [Default: Current working directory]"
)
parser$add_argument("-m",
                    "--name",
                    type = "character",
                    default = "peaks",
                    help = "Prefix for output file. [Default: peaks]")

args <- parser$parse_args()

print(args)

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tools))

makePeaksDF <-
  function(summit_file,
           peak_width,
           blocklist_bed_file,
           chrom_sizes_file,
           n,
           fdr_threshold) {
    # Make GRanges of peaks and blocklist
    pad <- round(as.numeric(peak_width) / 2)
    summit_file <- read.table(summit_file)
    peakdf <-
      summit_file[summit_file$V5 > -1 * log10(fdr_threshold), ]
    # Import sizes and remove chromosomes
    sizedf <-
      read.table(chrom_sizes_file, stringsAsFactors = FALSE)
    names(sizedf) <- c("seqnames", "end")
    sizedf$start <- 0
    chrranges <- makeGRangesFromDataFrame(sizedf)
    chrs <- as.character(sizedf[, 1])
    chrs <- chrs[chrs %ni% c("chrY", "chrM", "MT")]
    # Set peaks
    peakdf <- peakdf[peakdf[, 1] %in% chrs, ]
    peaks <-
      makeGRangesFromDataFrame(setNames(
        data.frame(peakdf[, 1], peakdf[, 2] - pad, peakdf[, 2] + pad, peakdf[, 5]),
        c("seqnames", "start", "end", "score")
      ), keep.extra.columns = TRUE)
    # Check for and fix out of bounds peaks; GenomicRanges
    # will Warn if out of bounds ranges exist
    chr_widths <- width(ranges(chrranges))
    names(chr_widths) <- seqnames(chrranges)
    seqlengths(peaks) <-
      chr_widths[match(names(seqlengths(peaks)), names(chr_widths))]
    peaks <-
      trim(peaks) # does the actual trimming, all out of bounds
    # ranges will now start or end at the end of the contig
    bdf <-
      data.frame(fread(input = blocklist_bed_file, header = FALSE))
    if (nrow(bdf > 0)) {
      bg <-
        makeGRangesFromDataFrame(setNames(data.frame(bdf[, 1], bdf[, 2],
                                                     bdf[, 3]),
                                                      c("seqnames", "start", "end")))
      # Remove blocklist and peaks out of bounds
      peaks <-
        #nolint start
        peaks[!(1:length(peaks) %in% data.frame(findOverlaps(peaks, bg))$queryHits)]
      #nolint end
      peaks <- subsetByOverlaps(peaks, chrranges, type = "within")
      peaks <- sortSeqlevels(peaks)
    }
    peaks <- sort(peaks)
    # Filter peaks based on summit score
    #nolint start
    keep_peaks <- 1:length(peaks)
    #nolint end
    while (!(isDisjoint(peaks[keep_peaks]))) {
      # Fast variable access
      chr_names <- as.character(seqnames(peaks[keep_peaks]))
      starts <- start(peaks[keep_peaks])
      ends <- end(peaks[keep_peaks])
      scores <- mcols(peaks)$score
      # See if consecutive peaks are overlapping
      overlap_next <-
        intersect(which(chr_names[1:(length(keep_peaks) - 1)] ==
                          chr_names[2:(length(keep_peaks))]),
                  which(ends[1:(length(keep_peaks) - 1)] >=
                   starts[2:(length(keep_peaks))]))
      # Compare consectuive peaks
      overlap_previous <- overlap_next + 1
      overlap_comparison <-
        scores[keep_peaks[overlap_previous]] > scores[keep_peaks[overlap_next]]
      discard <-
        keep_peaks[c(overlap_previous[!overlap_comparison],
                     overlap_next[overlap_comparison])]
      keep_peaks <- keep_peaks[keep_peaks %ni% discard]
    }
    peaks <- sortSeqlevels(peaks)
    peaks <- sort(peaks)
    # Export the final result by making a data frame; getting
    # the top (or as many) n peaks
    # based on the score and then resort based on genomic position.
    fP <- data.frame(peaks[keep_peaks], rank = seq_along(keep_peaks))
    nout <- min(as.numeric(n), dim(fP)[1])
    odf <- head(fP[order(fP$score, decreasing = TRUE), ], nout)
    return(odf)
  }

odf <-
  makePeaksDF(
    args$summit_file,
    args$peak_width,
    args$blocklist,
    args$chrom_sizes,
    n = args$n,
    fdr_threshold = args$fdr_threshold
  )

# Write to file
write.table(
  odf[sort(odf$rank, decreasing = FALSE, index.return = TRUE)$ix, c(1, 2, 3)],
  file = paste0(args$output_directory, "/", args$name, ".fixedwidthpeaks.bed"),
  col.names = FALSE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)
