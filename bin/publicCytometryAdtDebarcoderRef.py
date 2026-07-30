#!/usr/bin/env python3
"""
Script to convert CSV file with barcode sequences to JSON config format for use with the debarcoder.

CSV format expected:
SEQUENCE,DESCRIPTION
CTGGGCAATTACTCG,Hu.CD19
CTCATTGTAACTCCT,Hu.CD3
...

Output: JSON config file for barcode processing with debarcoder.
"""

import csv
import json
import click
import sys
from pathlib import Path


def read_csv_barcodes(csv_file):
    """
    Read barcodes from CSV file and return list of sequences.

    Args:
        csv_file (str): Path to CSV file

    Returns:
        list: List of barcode sequences from first column
    """
    barcodes = []

    try:
        with open(csv_file, 'r', newline='') as file:
            # Use csv.reader to handle potential commas in descriptions
            reader = csv.reader(file)

            for row_num, row in enumerate(reader, 1):
                if len(row) < 2:
                    print(f"Warning: Row {row_num} has fewer than 2 columns, skipping: {row}")
                    continue

                sequence = row[0].strip()
                if sequence:  # Skip empty sequences
                    barcodes.append(sequence)
                else:
                    print(f"Warning: Empty sequence in row {row_num}")

    except FileNotFoundError:
        print(f"Error: File '{csv_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        sys.exit(1)

    return barcodes


def create_json_config(barcodes):
    """
    Create JSON configuration with barcode sequences.

    Args:
        barcodes (list): List of barcode sequences

    Returns:
        dict: JSON configuration dictionary
    """
    config = {
        "primary_fastq": 0,
        "sequences": [
            {
                "type": "Barcode",
                "max_dist": 1,
                "values": barcodes
            },
            {
                "type": "Umi",
                "length": 1
            },
            {
                "type": "Const",
                "max_dist": 1,
                "value": "AAAAAAAAAA"
            }
        ]
    }

    return config


@click.command()
@click.argument('input_csv', type=click.Path(exists=True, path_type=Path))
@click.option('-o', '--output', type=click.Path(path_type=Path),
              help='Output JSON file (default: same name as input with .json extension)')
@click.option('--pretty', is_flag=True,
              help='Pretty print JSON with indentation')
@click.help_option('-h', '--help')
def main(input_csv, output, pretty):
    """
    Convert CSV barcode file to JSON config format for use with the debarcoder.

    INPUT_CSV: Input CSV file with barcode sequences in first column

    \b
    CSV format expected:
        SEQUENCE,DESCRIPTION
        CTGGGCAATTACTCG,Hu.CD19
        CTCATTGTAACTCCT,Hu.CD3
        ...

    \b
    Examples:
        python adt_debarcoder_ref.py input.csv -o output.json
        python adt_debarcoder_ref.py barcodes.csv
        python adt_debarcoder_ref.py antibodies.csv --pretty
    """

    # Determine output filename
    if output:
        output_file = output
    else:
        output_file = input_csv.with_suffix('.json')

    # Read barcodes from CSV
    click.echo(f"Reading barcodes from: {input_csv}")
    barcodes = read_csv_barcodes(str(input_csv))

    if not barcodes:
        click.echo("Error: No valid barcodes found in CSV file", err=True)
        sys.exit(1)

    click.echo(f"Found {len(barcodes)} barcodes")

    # Create JSON configuration
    config = create_json_config(barcodes)

    # Write JSON file
    try:
        with open(output_file, 'w') as file:
            if pretty:
                json.dump(config, file, indent=2)
            else:
                json.dump(config, file, indent="\t")

        click.echo(f"JSON config written to: {output_file}")
        click.echo(f"Barcode sequences included: {len(barcodes)}")

    except Exception as e:
        click.echo(f"Error writing JSON file: {e}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
