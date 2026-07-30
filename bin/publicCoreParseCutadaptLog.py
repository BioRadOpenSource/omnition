#!/usr/bin/env python3

from __future__ import annotations
import click
import polars as pl


@click.command()
@click.option(
    "--input-log",
    "-i",
    required=True,
    type=click.Path(exists=True),
    help="Input log file"
)
@click.option(
    "--sample-id",
    "-s",
    required=True,
    type=str,
    help="Name of sample"
)
def main(
    input_log,
    sample_id
):
    """
    Parse the trimming counts of each adapter from the cutadapt log file.
    Creates a CSV file with the counts of each adapter trimmed from each read.

    Args:
        input_log: Path to the cutadapt log file
        sample_id:  Name of the sample
    """
    adapter_counts = parse_cutadapt_log(input_log, sample_id)
    adapter_counts.write_csv(f"{sample_id}_cutadapt_trim_read_counts.csv")


def parse_cutadapt_log(input_log: str, sample_id: str) -> dict:
    """
    Parse the trimming counts of each adapter from the cutadapt log file

    Input read line example:
    Total read pairs processed:            148,738

    Output read line example:
    Pairs written (passing filters):       147,172 (98.9%)

    Adapter lines example:
    === Second read: Adapter Trueseq_P5_Forward ===

    Sequence: ACGACGCTCTTCCGATCT; Type: regular 5'; Length: 18; Trimmed: 73 times

    Args:
        input_log (str): Path to the cutadapt log file
        sample_id (str): Name of the sample

    Returns:
        pl.DataFrame: DataFrame of adapter count summary
    """
    adapter_counts = {
        "sample": [],
        "process": [],
        "metric": [],
        "value": []
    }
    with open(input_log, "r") as f:
        for line in f:
            # Pull out the total input and output read counts
            if "Total read pairs processed:" in line:
                adapter_counts["metric"].append("input")
                adapter_counts["value"].append(line.split(" ")[-1].strip().replace(",", ""))
            if "Pairs written (passing filters):" in line:
                adapter_counts["metric"].append("output")
                adapter_counts["value"].append(line.split(" ")[-2].replace(",", ""))
            # For adapters, we need to know which read they were trimmed from and which end
            if "Adapter" in line:
                start = line
                _ = next(f)
                data = next(f)
                adapter, count = parse_adapter(start, data)
                adapter_counts["metric"].append(adapter)
                adapter_counts["value"].append(count)
    adapter_counts["sample"] = [sample_id] * len(adapter_counts["metric"])
    adapter_counts["process"] = ["cutadapt_trim"] * len(adapter_counts["metric"])
    return pl.DataFrame(adapter_counts)


def parse_adapter(header: str, data: str) -> tuple:
    """
    Parse the adapter name and count of times it was trimmed

    Adapter header line example:
    === Second read: Adapter Trueseq_P5_Forward ===

    Adapter data line example:
    Sequence: ACGACGCTCTTCCGATCT; Type: regular 5'; Length: 18; Trimmed: 73 times

    Args:
        header (str): Header line of the adapter
        data (str): Data line of the adapter

    Returns:
        adapter (str): Name of the adapter
        count (str): Count of times the adapter was trimmed
    """
    if "First read:" in header:
        read = "R1"
    else:
        read = "R2"
    if "5';" in data:
        end = "5"
    else:
        end = "3"
    # Extract the adapter name and the count of times it was trimmed
    # See above for example formatting of the cutadapt log
    name = header.split("Adapter ")[1].split(" ")[0].replace(",", "")
    count = data.split("Trimmed: ")[1].split(" ")[0].replace(",", "")
    adapter = f"{name}_{read}_{end}"
    return adapter, count


if __name__ == "__main__":
    main()
