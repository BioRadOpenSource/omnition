#!/usr/bin/env python3

import click
import polars as pl


@click.command()
@click.option(
    "--filtered_bead_file",
    "-b",
    help="List of barcodes that do not have a cell_barcode in the DO file",
)
@click.option(
    "--barcode_translate_file",
    "-bt",
    help="Barcode Translate file made with cell barcodes "
    "found in the DO edgelist file.",
)
@click.option(
    "--sample_id",
    "-s",
    help="Sample be analyzed",
)
def main(filtered_bead_file, barcode_translate_file, sample_id):
    # Initialize sample Id
    sampleId = sample_id.replace("_", "")

    # Read in needed files
    bead_list = build_barcode_list(filtered_bead_file, sampleId)
    barcode_df = read_barcode_translate(barcode_translate_file)

    # Get original droplet count
    old_barcode_count = barcode_df.select(pl.col("DropBarcode")).n_unique()

    # Get the number of additional droplets from inclusion of unmerged beads
    additional_barcode_count = bead_list.select(pl.col("BeadBarcode")).n_unique()

    # Get updated total barcode count
    updated_barcode_count = old_barcode_count + additional_barcode_count

    # Get new padding value based on total barcodes
    padding_length = len(str(updated_barcode_count))

    # Generate new Cell Barcodes for barcodes without one in the edge file
    new_barcode_translate = format_bead_list(
        bead_list, old_barcode_count, updated_barcode_count, padding_length
    )

    reformated_barcode_translate = reformat_barcode_translate(
        barcode_df, padding_length
    )

    # Concat the reformated original to the new translate with unmerged beads included
    full_barcode_translate = pl.concat(
        [reformated_barcode_translate, new_barcode_translate]
    )

    full_barcode_translate.write_csv(f"{sample_id}_barcodeTranslate.tsv", sep=",")


def build_barcode_list(bead_file, sampleid):
    """
    Read in unmerged bead barcode list file and format as df
    """
    bead_barcodes = (
        pl.read_csv(bead_file, has_header=False).sort(by="column_1").to_series()
    )
    data = {
        "BeadBarcode": bead_barcodes,
        "sampleId": [sampleid] * len(bead_barcodes),
        "nBeads": [1] * len(bead_barcodes),
    }
    bead_list = pl.DataFrame(data)
    return bead_list


def read_barcode_translate(file):
    """
    Read barcode lookup table into polars dataframe
    Assumes columns are bead barcode, cell barcode, n beads in partition
    """
    barcode_df = pl.read_csv(file, sep=",", has_header=True)
    return barcode_df


def format_bead_list(
    bead_list, old_barcode_count, updated_barcode_count, padding_length
):
    """
    Format bead list to barcode_translate format
    """
    bead_list = bead_list.with_column(
        pl.arange(old_barcode_count + 1, updated_barcode_count + 1, eager=True)
        .cast(pl.Utf8)
        .str.zfill(padding_length)
        .alias("idNum")
    )

    bead_list = (
        bead_list.sort(by="BeadBarcode")
        .with_columns(
            pl.format(
                "{}BC{}N{}", pl.col("sampleId"), pl.col("idNum"), pl.col("nBeads")
            ).alias("DropBarcode"),
        )
        .select(["BeadBarcode", "DropBarcode", "nBeads"])
    )

    return bead_list


def reformat_barcode_translate(barcode_df, padding_length):
    """
    Update 0 padding of droplet IDs
    """
    barcode_df = barcode_df.with_columns(
        [
            barcode_df.select(pl.col("DropBarcode"))
            .to_series()
            .str.extract(pattern="BC([0-9]*)N")
            .str.zfill(padding_length)
            .alias("idNum"),
            barcode_df.select(pl.col("DropBarcode"))
            .to_series()
            .str.extract(pattern="^(.*)BC[0-9]*N")
            .alias("sampleId"),
        ]
    )

    barcode_df = (
        barcode_df.with_columns(
            pl.format(
                "{}BC{}N{}", pl.col("sampleId"), pl.col("idNum"), pl.col("nBeads")
            ).alias("DropBarcode"),
        )
        .select(["BeadBarcode", "DropBarcode", "nBeads"])
        .sort(["DropBarcode", "BeadBarcode"])
    )
    return barcode_df


if __name__ == "__main__":
    main()
