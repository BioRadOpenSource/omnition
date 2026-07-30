#!/usr/bin/env python3

import click
import pyfastx
import polars as pl


@click.command()
@click.option(
    "--fastq",
    "-f",
    type=click.Path(),
    required=True,
    help="FASTQ file."
)
@click.option(
    "--output",
    "-o",
    type=click.Path(),
    required=True,
    help="output file."
)
@click.option(
    "--sample-id",
    "-s",
    type=str,
    required=True,
    help="Sample ID."
)
def main(
    fastq,
    output,
    sample_id,
):
    """
    Generate a list of bead IDs from a FASTQ file and save it to a CSV file.

    Args:
        fastq (str): Path to the input FASTQ file.
        output (str): Path to the output CSV file.
        sample_id (str): Sample ID for the bead list.
    """
    barcodes = get_barcodes_from_fastq(sample_id, fastq)

    # Create a DataFrame with the barcodes
    df = pl.DataFrame({"barcode": barcodes})

    # Write the DataFrame to a CSV file
    df.write_csv(output, include_header=False, separator="\t")

    print(f"Bead list saved to {output}")


def get_barcodes_from_fastq(sample_id: str, fastq: str) -> list[str]:
    """
    Open up a FASTQ file for parsing and store the barcodes.

    Args:
        sample_id (str): The sample ID to use in the output file names.
        filename (str): The path to the input FASTQ file.
    Returns:
        List[str]: A list of barcodes extracted from the FASTQ file.
    """
    barcodes = set()

    # Iterate line by line through the fastq
    for name, seq, qual in pyfastx.Fastq(fastq, build_index=False, full_name=False):
        # First field is the cell barcode, 2nd is UMI, we only want the cell barcode
        barcode = "".join(name.split("_")[0])
        # Store the barcode in the list
        barcodes.add(barcode)

    return list(barcodes)


if __name__ == "__main__":
    main()
