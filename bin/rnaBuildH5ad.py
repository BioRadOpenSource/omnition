#!/usr/bin/env python3

import anndata as ad
import pandas as pd
import click


@click.command()
@click.option("--matrix", "-m", help="A .mtx(.gz) file.")
@click.option(
    "--barcodes", "-b", help="A .tsv(.gz) of barcodes corresponding to the matrix."
)
@click.option(
    "--features", "-f", help="A .tsv(.gz) file of features corresponding to the matrix."
)
@click.option(
    "--metadata", "-d", help="A .csv(.gz) file of metadata for barcodes in the matrix."
)
@click.option(
    "--embeddings",
    "-e",
    help="A .csv(.gz) of metadata created by Seurat, includes dim "
    "reduction embeddings.",
)
@click.option("--name", "-n", default="Sample", help="Sample name identifier.")
def main(matrix, barcodes, features, metadata, embeddings, name):
    h5 = one_sample(matrix, barcodes, features, metadata, embeddings, name)

    # write to h5
    h5.write("{}.h5ad".format(name))


def one_sample(matrix, barcodes, features, metadata, embeddings, name):
    # read matrix and row/colnames
    adata = ad.read_mtx(matrix).T
    barc = pd.read_csv(barcodes, header=None, sep="\t")
    features = pd.read_csv(features, header=None, sep="\t")

    # read embeddings and metadata
    embed = pd.read_csv(embeddings, header=0, sep=",")
    meta = pd.read_csv(metadata, header=0, sep=",")

    # filter metadata down to include only filtered barcodes (= those in .mtx)
    meta = meta[meta.barcode.isin(barc[0])]

    # join so that we're working on one df
    meta = meta.set_index("barcode").join(embed.set_index("barcode"), how="left")

    # start building h5ad
    # set variable features to gene symbols
    adata.var = features.drop(0, axis=1)

    # set the variable index to the Ensembl ID and naming the variable index
    adata.var.index = features[0]
    adata.var.index.name = "ensembl_id"

    # naming the variable column
    adata.var.columns = ["gene"]

    # set the observation index to the barcodes
    adata.obs.index = barc[0]

    # add an observation for each barcode, not saving the rest of the metadata
    adata.obs = meta.drop(list(meta.columns), axis=1)

    return adata


def add_embedding(h5ad, meta, embedding_prefix):
    h5ad.obsm["X_" + embedding_prefix.strip("_")] = meta.filter(
        regex=embedding_prefix
    ).values
    return h5ad


if __name__ == "__main__":
    main()
