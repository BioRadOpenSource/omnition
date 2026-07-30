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
    h5, assay_type = one_sample(matrix, barcodes, features, metadata, embeddings, name)

    # write to h5
    if assay_type == "RNA":
        ext = "_rna.h5ad"
    elif assay_type == "ADT":
        ext = "_cyto.h5ad"
    else:
        raise ValueError(f"Unknown assay_type: {assay_type}")
    h5.write(f"{name}{ext}")


def one_sample(matrix, barcodes, features, metadata, embeddings, name):
    # read matrix and row/colnames
    adata = ad.read_mtx(matrix).T
    barc = pd.read_csv(barcodes, header=None, sep="\t")
    features = pd.read_csv(features, header=None, sep="\t")

    # read embeddings and metadata
    embed = pd.read_csv(embeddings, header=0, sep=",")
    # Auto-detect assay type by looking for RNA or ADT substrings in column names
    embed_columns_str = ' '.join(embed.columns.tolist())
    has_adt = 'ADT' in embed_columns_str
    has_rna = 'RNA' in embed_columns_str

    # Determine assay type
    if has_rna and not has_adt:
        assay_type = "RNA"
    elif has_adt and not has_rna:
        assay_type = "ADT"
    else:
        # Default to RNA if no clear indicators
        assay_type = "RNA"
        print("Warning: Unable to determine assay. Defaulting to RNA.")

    print(f"Auto-detected assay type: {assay_type}")

    # Handle optional metadata file
    if metadata is not None:
        meta = pd.read_csv(metadata, header=0, sep=",")
        # filter metadata down to include only filtered barcodes (= those in .mtx)
        meta = meta[meta.barcode.isin(barc[0])]
        # join so that we're working on one df
        meta = meta.set_index("barcode").join(embed.set_index("barcode"), how="left")
    else:
        # If no metadata file provided, just use embeddings
        meta = embed.set_index("barcode")

    # start building h5ad
    # Handle features file structure - check if there are feature names beyond IDs
    if features.shape[1] > 1:
        # Features file has both IDs and names (typical for RNA data)
        adata.var = features.drop(0, axis=1)
        adata.var.index = features[0].astype(str)
        if assay_type == "ADT":
            adata.var.index.name = "antibody_id"
            adata.var.columns = ["antibody"]
        elif assay_type == "RNA":
            adata.var.index.name = "ensembl_id"
            adata.var.columns = ["gene"]
    else:
        # Features file only has one column (typical for ADT data with antibody names)
        feature_names = features[0].astype(str)
        if assay_type == "ADT":
            # For ADT: antibody names like "Hu.CD4" are both ID and display name
            adata.var = pd.DataFrame({"antibody": feature_names.values},
                                     index=pd.Index(feature_names.values, name="antibody_id"))
        elif assay_type == "RNA":
            # For RNA: gene IDs are used as both ID and display name
            adata.var = pd.DataFrame({"gene": feature_names.values},
                                     index=pd.Index(feature_names.values, name="ensembl_id"))

    # set the observation index to the barcodes
    adata.obs.index = barc[0]

    # add metadata as observations for each barcode
    # Convert all columns to strings to avoid H5AD writing errors
    meta = meta.astype(str)
    adata.obs = meta

    return adata, assay_type


def add_embedding(h5ad, meta, embedding_prefix):
    h5ad.obsm["X_" + embedding_prefix.strip("_")] = meta.filter(
        regex=embedding_prefix
    ).values
    return h5ad


if __name__ == "__main__":
    main()
