#!/usr/bin/env python3

import pyfastx
import argparse
import csv


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
        "-c",
        "--config",
        required=True,
        help="The csv-formatted pipeline config file to be used.",
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


# Imports and parses the TI config
def importConfig(config, sample):

    TIs = set()

    with open(config, "r") as f:
        readconfig = csv.DictReader(f)
        for row in readconfig:
            if row["fastq"] == sample:
                TIs.add(row["sequence"])
                # Create reads file for output later if this is a present TI
                readsfile = open(row["sequence"] + "_reads.txt", "w", newline="")
                csvwriter = csv.DictWriter(readsfile, fieldnames=row, lineterminator="")
                csvwriter.writerow(row)
                readsfile.close()

    return TIs


# Checks the TIs found in the fastq against the TIs listed in the config
def checkTIs(realTIs, configTIs, refTIs):

    errorlist = []
    for i in range(len(realTIs)):
        if realTIs[i] not in configTIs:
            errorlist.append(
                "\t[ERROR] TI "
                + refTIs[realTIs[i]]
                + " - "
                + realTIs[i]
                + " found in fastq but not found in config file. Please check the config file.\n"
            )
    for TI in configTIs:
        if TI not in realTIs:
            errorlist.append(
                "\t[ERROR] TI "
                + refTIs[TI]
                + " - "
                + TI
                + " found in config file but not found in fastq. Please check the config file.\n"
            )
    return errorlist


# Writes out error file if errors exist
def writeerrorfile(sample, errorlist, refTInames, realTIs, TIcounts):
    filename = sample + "_error.txt"
    F = open(filename, "w")
    F.write(sample + "\n")
    for i in range(len(errorlist)):
        F.write(errorlist[i])
    F.write("\n\tAll TIs found in " + sample + " fastq file\n")
    F.write("\tTI name,Sequence,Paired Read Count\n")
    for i in range(len(realTIs)):
        F.write(
            "\t"
            + refTInames[realTIs[i]]
            + ","
            + realTIs[i]
            + ","
            + str(TIcounts[i])
            + "\n"
        )
    F.write("\n")
    F.close()

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
configTIs = importConfig(args.config, args.sample_id)


# Check TIs
errors = checkTIs(realTIs, configTIs, refTInames)
if len(errors) > 0:
    writeerrorfile(args.sample_id, errors, refTInames, realTIs, TIcounts)

# add warnings if applicable
writewarnings(args.sample_id, refTInames, realTIs, TIcounts)

# Create a file for each TI in the config
for i in range(len(realTIs)):
    if realTIs[i] in configTIs:
        F = open(realTIs[i] + "_reads.txt", "a")
        F.write("," + str(TIcounts[i]))
        F.close()
