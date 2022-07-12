#!/usr/bin/env python3

import pyfastx
import argparse


# This function defines the input arguments.
def getargs():
    parser = argparse.ArgumentParser(
        description="Script to count the reads for each TI and check against the config."
    )

    options = parser.add_argument_group("Options")
    options.add_argument("-i", "--input-fastq", required=True, help="Input fastq file")
    options.add_argument(
        "-l", "--index-list", required=True, type=str, help="Reference list of indexes"
    )
    options.add_argument(
        "-s", "--sample-id", required=True, help="The Sample ID for the fastq file."
    )
    options.add_argument(
        "-m", "--messages", required=True, help="The messages file to add possible warnings to."
    )

    return parser.parse_args()


# Imports and parses the reference index list
def importTIindex(index):

    TIlist = list()
    TInames = dict()
    TIlen = 0

    with open(index, "r") as f:
        data = f.readline()
        data = f.readline()
        while data:
            if len(data) > 0:
                data = data.strip().split(",")
                TIlist.append(data[1])
                TInames[data[1]] = data[0]
                if len(data[1]) != TIlen:
                    if TIlen == 0:
                        TIlen = len(data[1])
                    else:
                        print(
                            "[ERROR] Reference TI sequences are of unequal lengths.  Please check pipeline parameter file."
                        )
                        exit(1)
            data = f.readline()

    return (TIlist, TIlen, TInames)


# Reads fastq file line by line and counts reads for each TI
def inputFastqs(filename, TIlen, refTIlist):
    TIsfound = []
    TIreadcounts = []

    for name, seq, qual in pyfastx.Fastq(filename, build_index=False):
        myTI = name.split("_")[0][(TIlen * -1) :]
        # Check that TIs found are in the reference list
        if myTI not in refTIlist:
            print(
                "[ERROR] TI "
                + myTI
                + " found in data not present in reference list. Please check run parameters."
            )
            exit(1)
        if myTI not in TIsfound:
            TIsfound.append(myTI)
            TIreadcounts.append(0)
        TIreadcounts[TIsfound.index(myTI)] += 1

    return (TIsfound, TIreadcounts)

# Add warnings to messages file if warnings exist
def writewarnings(sample, refTInames, realTIs, TIcounts):
    with open(args.messages, 'a') as f:
        for i in range(len(realTIs)):
            if TIcounts[i] < 10000:
                f.write("WARN: [ATAC] "
                    + "Sample,TI "
                    + sample
                    + ","
                    + refTInames[realTIs[i]]
                    + ","
                    + realTIs[i]
                    + " was detected to have low-counts. Count = "
                    + str(TIcounts[i])
                    + " (Threshold = 10000)"
                    + "\n"
                )

args = getargs()

# Load in input files
refTIlist, refTIlen, refTInames = importTIindex(args.index_list)
realTIs, TIcounts = inputFastqs(args.input_fastq, refTIlen, refTIlist)

# add warnings if applicable
writewarnings(args.sample_id, refTInames, realTIs, TIcounts)

# Create a file for each TI in the config
for i in range(len(realTIs)):
    F = open(realTIs[i] + "_reads.txt", "w")
    F.write(
        "SuperloadedSample,"
        + args.sample_id
        + ","
        + refTInames[realTIs[i]]
        + ","
        + realTIs[i]
        + ","
        + str(TIcounts[i])
    )
    F.close()
