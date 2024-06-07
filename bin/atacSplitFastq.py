#!/usr/bin/env python3

import pyfastx
import click
import ray
from ray.exceptions import ObjectStoreFullError
import sys


# This function defines the input arguments.
@click.command()
@click.option("--read1-fastq", "-r1", help="Input R1 fastq file")
@click.option("--read2-fastq", "-r2", help="Input R1 fastq file")
@click.option("--sample-id", "-s", help="Path of cell barcode allowlist CSV file.")
@click.option("--counts-file", "-c", help="Path of cell barcode allowlist CSV file.")
@click.option(
    "--cpus", "-p", default=1, help="The number of CPUs for parallel processing."
)
@click.option(
    "--override_errors", "-o", is_flag=True,
    help="When present, FASTQ-TIs not in index file are ignored."
)
def main(read1_fastq, read2_fastq, sample_id, counts_file, override_errors, cpus):
    # Read in indexes
    indexes = parse_counts_file(counts_file, sample_id)

    # Parse R1 and R2 in parallel
    ray.init(num_cpus=int(cpus))
    try:
        ray.get(
            [
                parse_fastq.remote(sample_id, readfile, indexes, override_errors)
                for readfile in [read1_fastq, read2_fastq]
            ]
        )
    except ObjectStoreFullError:
        ray.shutdown()
        print("OOM error: ObjectStoreFullError")
        sys.exit(137)
    except Exception as err:
        ray.shutdown()
        print(f"Unexpected {err=}, {type(err)=}")
        sys.exit(1)
    ray.shutdown()


def parse_counts_file(filename: str, sampleID: str) -> [str]:
    """
    Parses a counts file and returns a list of indexes for a given sample ID.

    Args:
        filename (str): The path to the counts file.
        sampleID (str): The ID of the sample to extract indexes for.

    Returns:
        list: A list of indexes for the given sample ID.

    """
    indexes = []
    with open(filename, "r") as f:
        for data in f:
            columns = data.split(",")
            if len(columns) > 1:
                if columns[1] == sampleID:
                    indexes.append(columns[3])
    return indexes


@ray.remote
def parse_fastq(sample_id: str, filename: str, indexes: [str],
                override_errors: bool) -> None:
    """
    Split a FASTQ file into multiple files based on the index sequence.

    Args:
        sample_id (str): The sample ID to use in the output file names.
        filename (str): The path to the input FASTQ file.
        indexes (list of str): The list of index sequences to split the file by.
        override_errors (bool): Whether to skip sequences with unknown indexes.

    Returns:
        None
    """
    if filename.find("_R2_001.complete_debarcoded.fastq.gz") == -1:
        outfile = sample_id + "-{}_R1.complete_debarcoded.split.fastq"
    else:
        outfile = sample_id + "-{}_R2.complete_debarcoded.split.fastq"
    # Open a new file for each index
    outfiles = [open(outfile.format(index), "w") for index in indexes]
    # Iterate line by line through the fastq
    for name, seq, qual in pyfastx.Fastq(filename, build_index=False, full_name=True):
        # Pull the index sequence
        myTI = name.split("_")[0][(len(indexes[0]) * -1):]
        # If the index is not in the index file, skip it
        if override_errors and myTI not in indexes:
            continue
        # Write to the correct file
        outfiles[indexes.index(myTI)].write("@" + name + "\n")
        outfiles[indexes.index(myTI)].write(seq + "\n")
        outfiles[indexes.index(myTI)].write("+\n")
        outfiles[indexes.index(myTI)].write(qual + "\n")
    # Close all open files
    for i in outfiles:
        i.close()


if __name__ == "__main__":
    main()
