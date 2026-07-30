#!/usr/bin/env python3

import argparse
import os

parser = argparse.ArgumentParser(description="Compile deduplication statistics.")

# Defining argument group 'Options'
options = parser.add_argument_group("Options")
options.add_argument(
    "-s", "--sample-id", required=True, help="The sample ID of the compiled sample."
)
options.add_argument(
    "-i",
    "--input-directory",
    required=True,
    help="Input directory containing dedup_stats.txt files.",
)


# Make a list of all dedup stats filenames
def grabDedupFiles(inputdir):
    dedupfiles = []
    dirfiles = os.listdir(inputdir)
    for i in range(len(dirfiles)):
        if dirfiles[i].find("dedup_stats.txt") != -1:
            dedupfiles.append(inputdir + dirfiles[i])

    return dedupfiles


# Parse the mark duplicate output files and create a compiled file
def parseAndSum(sampleid, dedupfiles):
    # Open the output files
    dedupdata = []
    for i in range(len(dedupfiles)):
        with open(dedupfiles[i], "r") as f:
            data = f.readline()
            # Grab the ## METRICS section which is used in
            # the reports and aggregate metrics
            while data.find("## METRICS") == -1:
                data = f.readline()
            header1 = data
            header2 = f.readline()
            data = f.readline()
            dedupdata.append(data.strip().split("\t"))

    outdata = [sampleid]
    for i in range(1, len(dedupdata[0])):
        # Sum together all columns except percent duplication
        if i != 8:
            totalsum = 0
            for j in range(len(dedupdata)):
                totalsum += int(dedupdata[j][i])
            outdata.append(totalsum)
        else:
            # Percent duplicates = (2*READ_PAIR_DUPLICATES+UNPAIRED_READ_DUPLICATES)/
            # (2*READ_PAIRS_EXAMINED+UNPAIRED_READS_EXAMINED)
            percdup = ((outdata[6] * 2) + outdata[5]) / ((outdata[2] * 2) + outdata[1])
            outdata.append(percdup)

    # Write the compiled numbers tp a new file
    outfile = sampleid + ".dedup_stats.txt"
    with open(outfile, "w") as f:
        f.write(header1)
        f.write(header2)
        outdata = [str(x) for x in outdata]
        f.write("\t".join(outdata) + "\n")


# Executing script
if __name__ == "__main__":
    # Establish inputs
    args = parser.parse_args()
    dedupfiles = grabDedupFiles(args.input_directory)

    # Parse the mark duplicate output files and create a compiled file
    parseAndSum(args.sample_id, dedupfiles)
