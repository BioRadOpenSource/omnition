#!/usr/bin/env python3

import argparse
import os
import csv


parser = argparse.ArgumentParser(
    description="Validate sample TI and fastq configuration."
)

# Defining argument group 'Options'
options = parser.add_argument_group("Options")
options.add_argument(
    "-c",
    "--config",
    required=True,
    help="The csv-formatted pipeline config file to be used.",
)
options.add_argument(
    "-l",
    "--sample-list",
    required=True,
    help="Text file containing list of sample IDs to be analyzed.",
)
options.add_argument(
    "-i",
    "--index-list",
    required=True,
    help="Text file containing list of reference tagmentation indexes.",
)
options.add_argument(
    "-o",
    "--override",
    required=True,
    help="Override parameter to ignore fastqs not being used in the config",
)


# Imports and parses the TI config
def importConfig(config):

    samples = list()
    fastqs = set()
    TIs = set()
    TIfastqPairs = list()
    missing_fastq = list()

    with open(config, "r") as f:
        readconfig = csv.DictReader(f)
        for row in readconfig:
            if "sample" not in row.keys():
                print(
                    "[ERROR] The config file should be a CSV formatted with three columns: sample, fastq, and ti",
                    ", ".join(missing_fastq),
                )
                exit(1)
            if "fastq" not in row.keys():
                print(
                    "[ERROR] The config file should be a CSV formatted with three columns: sample, fastq, and ti",
                    ", ".join(missing_fastq),
                )
                exit(1)
            if "ti" not in row.keys():
                print(
                    "[ERROR] The config file should be a CSV formatted with three columns: sample, fastq, and ti",
                    ", ".join(missing_fastq),
                )
                exit(1)
            samples.append(row)
            fastqs.add(row["fastq"])
            TIs.add(row["ti"])
            TIfastqPairs.append(row["fastq"] + "-" + row["ti"])

    return (samples, fastqs, TIs, TIfastqPairs)


# Imports and parses the reference index list
def importTIindex(index):

    TIdict = dict()
    TIlist = list()

    with open(index, "r") as f:
        data = f.readline()
        while data:
            data = data.strip().split(",")
            TIlist.append(data[0])
            TIdict[data[0]] = data[1]
            data = f.readline()

    return (TIdict, TIlist)


# Function for importing list of provided sample IDs
def importSampleList(sample_list_arg):

    with open(sample_list_arg, "r") as f:
        sample_string = f.readline()
        sample_list = (
            sample_string.strip()
            .replace("[", "")
            .replace("]", "")
            .replace(" ", "")
            .split(",")
        )

    return sample_list


# Function for checking if the fastqs in the config file are present and account for all fastqs provided
def checkSamples(config_fastqs, sample_list, override):

    # Check if any fastqs in the config file are missing
    missing_fastq = list(config_fastqs - set(sample_list))
    if len(missing_fastq) != 0:
        print(
            "[ERROR] The following fastqs were specified in the config but not found:",
            ", ".join(missing_fastq),
        )
        exit(1)

    # Check if any fastqs present are not in the config file
    unaccounted_fastq = list(set(sample_list) - config_fastqs)
    if override != "true":
        if len(unaccounted_fastq) != 0:
            print(
                "[ERROR] The following fastqs were not assigned to a sample in the config:",
                ", ".join(unaccounted_fastq),
            )
            exit(1)


# Function for checking if each TI fastq pair is used only once
def checkTIfastqPairs(TIfastqPairs):

    check = set(TIfastqPairs)

    if len(check) != len(TIfastqPairs):
        multiples = []
        for item in check:
            if TIfastqPairs.count(item) > 1:
                multiples.append(item)
        print(
            "[ERROR] The following TI and fastq pairs were used more than once in the config:",
            ", ".join(multiples),
        )
        exit(1)


# Function for checkng if the TIs listed in the config are valid TIs
def checkTIindex(TIs, TIlist):

    errors = []

    for index in TIs:
        if index not in TIlist:
            errors.append(index)

    if len(errors) != 0:
        print("[ERROR] The following TIs are not valid TIs:", ", ".join(errors))
        exit(1)


# Create an output csv config file with the sequence of each TI added as a new column
def outputCSV(outfile, csvdict, TIdict):
    with open(outfile, "w", newline="") as file:
        fieldnames = ["sample", "fastq", "ti", "sequence"]
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        for row in csvdict:
            row["sequence"] = TIdict[row["ti"]]
            writer.writerow(row)


# Executing script
if __name__ == "__main__":

    # Establish inputs
    args = parser.parse_args()
    outfile = os.path.basename(args.config).split(".")[0] + ".validated.csv"

    # Import and parse inputs
    cfsamples, cffastqs, cfTIs, cfTIfastqPairs = importConfig(args.config)
    sample_list = importSampleList(args.sample_list)
    indexdict, indexlist = importTIindex(args.index_list)

    # Sanity Checks
    checkSamples(cffastqs, sample_list, args.override)
    checkTIindex(cfTIs, indexlist)
    checkTIfastqPairs(cfTIfastqPairs)

    # Create outputs
    outputCSV(outfile, cfsamples, indexdict)
