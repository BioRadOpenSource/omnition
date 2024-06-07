#!/usr/bin/env python3

import argparse
import statistics as stats
import pysam


def getargs():
    parser = argparse.ArgumentParser(
        description="Script to count CBCs with biased UMIs."
    )

    options = parser.add_argument_group("Options")
    options.add_argument(
        "-i",
        "--input-bam-file",
        required=True,
        help="Deduplicated bam file with XM and XC tags",
    )
    options.add_argument("-s", "--sample-id", required=True, help="Sample ID")
    options.add_argument(
        "-t", "--barcode-tag", required=True, help="Tag used for barcodes."
    )
    options.add_argument("-u", "--umi-tag", required=True, help="Tag used for UMIs.")

    return parser.parse_args()


def inputbam(filename, barcode_tag, umi_tag):
    # Import the bam file
    inbam = pysam.AlignmentFile(filename, "rb")
    umidict = dict()
    nts = ["A", "C", "G", "T", "N"]
    startreads = 0

    for read in inbam.fetch():
        # Pull CBC and UMI tags
        if read.has_tag(barcode_tag):
            mycbc = read.get_tag(barcode_tag)
            myumi = read.get_tag(umi_tag)
            # Create a dictionary to record nucleotide content
            # at each position for each umi
            if umidict.get(mycbc) is None:
                umidict[mycbc] = [[0, 0, 0, 0, 0] for x in range(len(myumi))]
            for j in range(len(myumi)):
                umidict[mycbc][j][nts.index(myumi[j])] += 1
            startreads += 1
    startcbcs = len(umidict.keys())

    return (umidict, startcbcs, startreads)


def countbiasedumis(sample, umidict, inputfile, startcbcs, startreads):
    outfile = inputfile.replace(".bam", "_biased_umis.csv")

    halves = []
    wholes = []
    allcbcs = umidict.keys()
    for i in allcbcs:
        # Count the number of biased NTs (>80% biased at that position)
        failcount = 0
        for j in range(len(umidict[i])):
            if max(umidict[i][j]) >= sum(umidict[i][j]) * 0.8:
                failcount += 1
        umicount = sum(umidict[i][0])
        if failcount >= len(umidict[i]) / 2:
            halves.append(umicount)
        if failcount == len(umidict[i]):
            wholes.append(umicount)
    # Write the output
    F = open(outfile, "w")
    F.write(
        "Sample,Input_CBCs,Input_UMIs,50_Filtered_CBCs,50_Total_Filtered_UMIs,"
        "50_Max_UMIs_per_CBC,"
        "50_Mean_UMIs_per_CBC,50_Median_UMIs_per_CBC,100_Filtered_CBCs,"
        "100_Total_Filtered_UMIs,"
        "100_Max_UMIs_per_CBC,100_Mean_UMIs_per_CBC,100_Median_UMIs_per_CBC\n"
    )
    F.write(
        "{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(
            sample,
            startcbcs,
            startreads,
            len(halves),
            sum(halves),
            max(halves),
            stats.mean(halves),
            stats.median(halves),
            len(wholes),
            sum(wholes),
            max(wholes),
            stats.mean(wholes),
            stats.median(wholes),
        )
    )
    F.close()


if __name__ == "__main__":
    args = getargs()
    cbcs, startcbcs, startreads = inputbam(
        args.input_bam_file, args.barcode_tag, args.umi_tag
    )
    countbiasedumis(args.sample_id, cbcs, args.input_bam_file, startcbcs, startreads)
