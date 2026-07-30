#!/usr/bin/env python3

import polars as pl
import click
from scipy.sparse import csr_matrix
import gzip


@click.command()
@click.option(
    "--barcodes",
    "-b",
    help="barcode list file",
    required=True
)
@click.option(
    "--adts",
    "-a",
    help="gene list file",
    required=True
)
@click.option(
    "--matrix",
    "-m",
    help="sparse matrix file",
    required=True
)
@click.option(
    "--translate",
    "-t",
    help="barcode translation file",
    required=True
)
@click.option(
    "--sample",
    "-s",
    help="Sample ID string to use for labelling outputs.",
    required=True
)
def main(matrix, adts, barcodes, translate, sample):
    sparse_mtx, adts, barcodes, translate = read_input_files(matrix, adts, barcodes, translate)
    full_mtx = make_full_matrix(adts, barcodes, sparse_mtx)
    droplet_mtx = merge_droplet_barcodes(full_mtx, translate)
    transposed_mtx = transpose_mtx(droplet_mtx)
    old_create_and_save_sparse_mtx(droplet_mtx)
    create_and_save_sparse_mtx(transposed_mtx, sample)


def read_input_files(matrix, adts, barcodes, translate):
    """
    Read in the input files for the sparse matrix, adts, barcodes, and the barcode translate file
    """
    sparse_mtx = pl.read_csv(matrix, separator=' ', skip_rows=3, new_columns=["barcode", "adt", "count"])
    adts = pl.read_csv(adts, has_header=False, new_columns=["adt_seq"])
    barcodes = pl.read_csv(barcodes, has_header=False, new_columns=["BeadBarcode"])
    translate = pl.read_csv(translate, separator=',')

    return sparse_mtx, adts, barcodes, translate


def make_full_matrix(adts, barcodes, sparse_mtx):
    """
    Recreate full matrix from sparse matrix files
    """
    # Add row to translate between sparse matrix and axis labels
    adts = adts.with_columns((pl.arange(0, adts.height) + 1).alias("adt"))
    barcodes = barcodes.with_columns((pl.arange(0, barcodes.height) + 1).alias("barcode"))

    # Join the dataframes on the adt and barcode columns
    sparse_mtx = sparse_mtx.join(adts, on="adt")
    sparse_mtx = sparse_mtx.join(barcodes, on="barcode")

    # Pivot the sparse matrix to a full matrix and fill null values with 0
    full_mtx = sparse_mtx.pivot(index='BeadBarcode', on='adt_seq', values='count').fill_null(strategy="zero")
    return full_mtx


def merge_droplet_barcodes(full_mtx, translate):
    """
    Merge the full matrix with the barcode translation file and group by the droplet barcode column
    """
    # Merge the dataframes on the barcode column
    combined_mtx = full_mtx.join(translate, on="BeadBarcode")
    # Group by the droplet barcode column and sum the counts from individual beads
    droplet_mtx = combined_mtx.group_by("DropBarcode").sum().drop("BeadBarcode", "nBeads")

    return droplet_mtx


def transpose_mtx(droplet_mtx):
    # transpose the mtx
    transposed_mtx = droplet_mtx.transpose(include_header=True)

    # Extract the first row to use as the new header
    new_header = list(transposed_mtx.row(0))

    # Modify the first element of the new header
    new_header[0] = "ADTs"

    # Remove the first row from the DataFrame
    transposed_mtx = transposed_mtx[1:]

    # Set the new header for the DataFrame
    transposed_mtx.columns = new_header

    return transposed_mtx


def old_create_and_save_sparse_mtx(droplet_mtx):
    """
    Create and save the sparse matrix file
    """
    # Save row/column values to text files
    drop_barcodes = droplet_mtx['DropBarcode'].to_list()
    adt_list = droplet_mtx.columns
    adt_list.remove('DropBarcode')

    # Write the row indices to a text file
    with open("merged.barcodes.txt", "w") as file:
        for barcode in drop_barcodes:
            file.write(f"{barcode}\n")

    # Write the column headers to a text file
    with open("merged.genes.txt", "w") as file:
        for adt in adt_list:
            file.write(f"{adt}\n")

    # Make the sparse matrix
    sparse_mtx = csr_matrix(droplet_mtx.drop('DropBarcode').to_numpy())

    # Get the row and column indices and index them at 0 for downstream compatibility
    rows, cols = sparse_mtx.nonzero()
    rows += 1
    cols += 1
    # Get the data values
    data = sparse_mtx.data

    # Format data for sparse matrix file
    combined = pl.DataFrame({
        'Column1': rows,
        'Column2': cols,
        'Column3': data
    })
    new_row = pl.DataFrame({
        'Column1': [combined['Column1'].max()],
        'Column2': [combined['Column2'].max()],
        'Column3': [len(combined)]
    })
    combined = pl.concat([new_row, combined], how="vertical_relaxed")

    # Write expected first line to the text file
    with open('merged.mtx', 'w') as f:
        f.write('%%MatrixMarket matrix coordinate real general\n%\n%                                               \n')
        # Append the DataFrame to the file
        combined.write_csv(f, separator=' ', include_header=False)


def create_and_save_sparse_mtx(transposed_mtx, sample):
    """
    Create and save the sparse matrix file
    """
    # Save row/column values to text files
    adt_list = transposed_mtx['ADTs'].to_list()
    drop_barcodes = transposed_mtx.columns
    drop_barcodes.remove('ADTs')

    # Write the row indices to a text file
    with open(sample + ".unfiltered.barcodes.tsv", "w") as file:
        for barcode in drop_barcodes:
            file.write(f"{barcode}\n")

    # Write the column headers to a text file
    with open(sample + ".unfiltered.genes.tsv", "w") as file:
        for adt in adt_list:
            file.write(f"{adt}\n")

    # Make the sparse matrix
    sparse_mtx = csr_matrix(transposed_mtx.drop('ADTs').cast({pl.String: pl.Int64}).to_numpy())

    # Get the row and column indices and index them at 0 for downstream compatibility
    rows, cols = sparse_mtx.nonzero()
    rows += 1
    cols += 1
    # Get the data values
    data = sparse_mtx.data

    # Format data for sparse matrix file
    combined = pl.DataFrame({
        'Column1': rows,
        'Column2': cols,
        'Column3': data
    })
    new_row = pl.DataFrame({
        'Column1': [combined['Column1'].max()],
        'Column2': [combined['Column2'].max()],
        'Column3': [len(combined)]
    })
    combined = pl.concat([new_row, combined], how="vertical_relaxed")

    # Write expected first line to the gzipped file
    with gzip.open(sample + '.unfiltered.mtx.gz', 'wt') as f:
        f.write('%%MatrixMarket matrix coordinate real general\n%\n%                                               \n')
        # Convert the DataFrame to a CSV string
        csv_string = combined.write_csv(separator=' ', include_header=False)
        # Write the CSV string to the gzipped file
        f.write(csv_string)


if __name__ == "__main__":
    main()
