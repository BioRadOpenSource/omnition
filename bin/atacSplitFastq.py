#!/usr/bin/env python3

import pyfastx
import argparse
import os


# This function defines the input arguments.
def getargs():
    parser = argparse.ArgumentParser(
        description="Script to pull reads for a TI from a fastq file."
    )

    options = parser.add_argument_group("Options")
    options.add_argument(
        "-r1", "--read1-fastq", required=True, help="Input R1 fastq file"
    )
    options.add_argument(
        "-r2", "--read2-fastq", required=True, help="Input R2 fastq file"
    )
    options.add_argument(
        "-i", "--index", required=True, type=str, help="Tagmentation index"
    )

    return parser.parse_args()


# Reads fastq file line by line and writes out reads matching the index
def parseFastq(filename, outfile, index):

    F = open(outfile, "w")
    for name, seq, qual in pyfastx.Fastq(filename, build_index=False, full_name=True):
        myTI = name.split("_")[0][(len(index) * -1) :]
        if myTI == index:
            F.write("@" + name + "\n")
            F.write(seq + "\n")
            F.write("+\n")
            F.write(qual + "\n")

    F.close()


args = getargs()

# Set up variables
sample_id = os.path.basename(
    args.read1_fastq.replace("_R1_001.complete_debarcoded.fastq.gz", "")
)
outfile1 = sample_id + "-" + args.index + "_R1.complete_debarcoded.split.fastq"
outfile2 = sample_id + "-" + args.index + "_R2.complete_debarcoded.split.fastq"

# Parse R1 and R2
parseFastq(args.read1_fastq, outfile1, args.index)
parseFastq(args.read2_fastq, outfile2, args.index)
