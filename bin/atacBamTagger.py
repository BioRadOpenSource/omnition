#!/usr/bin/env python3

import pysam
import click


@click.command()
@click.option(
    "--cell-barcode-tag", "-ct", default="XC", help="Tag to use for the cell barcode."
)
@click.option(
    "--cell-barcode-position",
    "-cp",
    default=0,
    help="In an underscore delimited read name, the zero-based field that contains the cell barcode.",
)

def main(cell_barcode_tag,cell_barcode_position):
    infile = pysam.AlignmentFile("-", "rb")
    outfile = pysam.AlignmentFile("-", "w", template=infile)
    cell_barcode_position = 0
    cell_barcode_tag = 'XB'

    for line in infile:
        rname = line.query_name.split("_")

        barcode = rname[cell_barcode_position]
        readname = rname[cell_barcode_position + 1 :]

        line.set_tag(cell_barcode_tag, barcode)

        line.query_name = readname[0]
        outfile.write(line)

if __name__ == "__main__":
    main()
