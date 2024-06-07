#!/usr/bin/env python3

import click


@click.command()
@click.option(
    "--duplicates_file",
    "-d",
    help="Duplicate counts file",
)
@click.option(
    "--barcode_translate_file",
    "-b",
    help="Bead and drop barcode translate file.",
)
def main(duplicates_file, barcode_translate_file):
    # Read barcode translate
    translate_dict = read_barcode_translate(barcode_translate_file)

    # Analyze the duplicate barcodes
    duplicates_dict = analyze_duplicates(translate_dict, duplicates_file)

    # Write outputs
    write_merged_duplicates(duplicates_file, duplicates_dict)


def read_barcode_translate(file):
    """
    Read the barcode translate file for later use.
    """
    # Initialize output dictionary
    output_dict = {}

    # Read in barcode translate file setting bead barcode as keys
    with open(file) as f:
        for line in f:
            bead_barcode, drop_barcode, nBeads = line.strip().split(",")
            output_dict[bead_barcode] = drop_barcode

    # Return the translation dictionary
    return output_dict


def analyze_duplicates(barcode_translate, duplicates_file):
    """
    Read through duplicate file and convert to drop barcode
    """
    # Initialize output dictionary
    duplicates_dict = {}

    # Convert entries in duplicates file to drop barcode
    with open(duplicates_file) as f:
        for line in f:
            bead_barcode, gene, count = line.strip().split(",")
            # If the bead barcode is not translated, keep as is
            if bead_barcode not in barcode_translate.keys():
                # Add a dictionary for every bead barcode
                if bead_barcode not in duplicates_dict.keys():
                    duplicates_dict[bead_barcode] = {}
                # Add a counter for every gene
                if gene not in duplicates_dict[bead_barcode].keys():
                    duplicates_dict[bead_barcode][gene] = 0
                # Add read counts for every gene
                duplicates_dict[bead_barcode][gene] += int(count)
            # Save gene counts under drop barcode
            else:
                # Add a dictionary for every drop barcode
                if barcode_translate[bead_barcode] not in duplicates_dict.keys():
                    duplicates_dict[barcode_translate[bead_barcode]] = {}
                # Add a counter for every gene
                if gene not in duplicates_dict[barcode_translate[bead_barcode]].keys():
                    duplicates_dict[barcode_translate[bead_barcode]][gene] = 0
                # Add read counts for every gene
                duplicates_dict[barcode_translate[bead_barcode]][gene] += int(count)

    # Return the new dictonary
    return duplicates_dict


def write_merged_duplicates(duplicates_file, duplicates_dict):
    """
    Write out translated duplicate counts file
    """
    outfile = duplicates_file.replace("duplicate", "merged_duplicate")
    with open(outfile, "w") as f:
        for barcode in duplicates_dict.keys():
            for gene in duplicates_dict[barcode].keys():
                f.write(
                    "{},{},{}\n".format(barcode, gene, duplicates_dict[barcode][gene])
                )


if __name__ == "__main__":
    main()
