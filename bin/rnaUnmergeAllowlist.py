#!/usr/bin/env python3

import click


@click.command()
@click.option(
    "--allowlist_file",
    "-a",
    help="Allowlist file",
)
@click.option(
    "--barcode_translate_file",
    "-b",
    help="Bead and drop barcode translate file.",
)
def main(allowlist_file, barcode_translate_file):
    # Get barcode translate
    translate_dict = read_barcode_translate(barcode_translate_file)

    # Get allowlist
    allowlist = read_allowlist(allowlist_file)

    # Convert cell barcodes to bead barcodes
    unmerged_allowlist = convert_barcodes(allowlist, translate_dict)

    # Write outputs
    write_unmerged_allowlist(allowlist_file, unmerged_allowlist)


def read_barcode_translate(bt_file):
    """
    Iterate over barcode translate file and read
    into a dictonary where drop_barcodes are keys
    """
    # Initialize output dictonary
    translate_dict = {}

    # Iterate over file and update output dictonary
    # Read in barcode translate file with drop barcodes as keys
    with open(bt_file) as f:
        for line in f:
            bead_barcode, drop_barcode, nBeads = line.strip().split(",")
            # Add a set for every drop barcode to store all bead barcodes
            if drop_barcode not in translate_dict.keys():
                translate_dict[drop_barcode] = set()
            translate_dict[drop_barcode].add(bead_barcode)
    # Return compiled dict
    return translate_dict


def read_allowlist(allowlist_file):
    """
    Read in allowlist drop barcodes into a set
    """
    # Initialize output set
    allowlist = set()

    # Read allowlist drop barcoedes into set
    with open(allowlist_file) as f:
        for line in f:
            allowlist.add(line.strip())
    return allowlist


def convert_barcodes(allowlist, translate_dict):
    """
    Convert the cell barcodes from the allowlist into bead
    barcodes.
    """
    # Initialize output set
    unmerged_allowlist = set()

    # Convert drop barcodes to bead barcodes
    for barcode in allowlist:
        if barcode not in translate_dict.keys():
            unmerged_allowlist.add(barcode)
        else:
            unmerged_allowlist.update(translate_dict[barcode])
    # Return new set of barcodes
    return unmerged_allowlist


def write_unmerged_allowlist(allowlist_file, unmerged_allowlist):
    """
    Write out bead barcode allowlist
    """
    outfile = allowlist_file.replace("allowlist", "unmerged_allowlist")
    with open(outfile, "w") as f:
        for bead in unmerged_allowlist:
            f.write("{}\n".format(bead))


if __name__ == "__main__":
    main()
