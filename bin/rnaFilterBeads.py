#!/usr/bin/env python3

import click
import polars as pl


@click.command()
@click.option(
    "--edge_file",
    "-e",
    type=click.Path(exists=True),
    required=True,
    help="Edgelist file to filter.",
)
@click.option(
    "--barcodes",
    "-b",
    type=click.Path(exists=True),
    required=True,
    help="Barcodes with which to filter.",
)
@click.option(
    "--sample_id",
    "-s",
    type=str,
    required=False,
    help="Sample ID for formatting output files.",
)
def main(edge_file, barcodes, sample_id):
    # read allowlist
    allowlist = read_allowlist(barcodes, True)

    # read edgelist
    edgelist = read_edgelist(edge_file)

    # filter edges
    filtered_edgelist = filter_edgelist(edgelist.lazy(), allowlist.lazy())

    # write output
    filtered_edgelist.write_csv(f"{sample_id}_edges.tsv", has_header=False, sep="\t")


def read_allowlist(allowlist_file, merged):
    """
    Read allowlist file into a dataframe
    """
    if merged:
        allowlist = pl.read_csv(
            allowlist_file, sep="\t", has_header=False, new_columns=["droplet"]
        ).with_row_count(name="rank", offset=1)
    else:
        allowlist = pl.read_csv(
            allowlist_file, sep="\t", has_header=False, new_columns=["bead"]
        ).with_row_count(name="rank", offset=1)
    return allowlist


def read_edgelist(edge_file):
    """
    Read edgelist file into a dataframe, removing duplicate rows.
    Order such that bead1 and umi1 are alphabetically first.
    """
    cols = [0, 1, 2, 3]
    new_cols = ["bead1", "umi1", "bead2", "umi2"]
    df = (
        pl.read_csv(
            edge_file,
            has_header=False,
            sep="\t",
            columns=cols,
            new_columns=new_cols,
        )
        .unique()
        .lazy()
    )
    df = df.with_columns(
        [
            # Reorder beads alphabetically
            pl.when(pl.col("bead1") < pl.col("bead2"))
            .then(pl.col("bead1") + "_" + pl.col("bead2"))
            .otherwise(pl.col("bead2") + "_" + pl.col("bead1"))
            .alias("bead1_bead2"),
            # Reorder umis alphabetically
            pl.when(
                (pl.col("bead1") == pl.col("bead2")) & (pl.col("umi1") < pl.col("umi2"))
            )
            .then(pl.col("umi1") + "_" + pl.col("umi2"))
            .when(
                (pl.col("bead1") == pl.col("bead2")) & (pl.col("umi1") > pl.col("umi2"))
            )
            .then(pl.col("umi2") + "_" + pl.col("umi1"))
            .when(
                (pl.col("bead1") == pl.col("bead2"))
                & (pl.col("umi1") == pl.col("umi2"))
            )
            .then(pl.col("umi1") + "_" + pl.col("umi2"))
            .when(pl.col("bead1") < pl.col("bead2"))
            .then(pl.col("umi1") + "_" + pl.col("umi2"))
            .otherwise(pl.col("umi2") + "_" + pl.col("umi1"))
            .alias("umi1_umi2"),
        ]
    )
    expr = [
        pl.col("bead1_bead2").alias("edge_id"),
        pl.col("umi1_umi2").alias("edge_umi"),
        pl.col("bead1_bead2").str.split("_").arr.get(0).alias("bead1_reorder"),
        pl.col("bead1_bead2").str.split("_").arr.get(1).alias("bead2_reorder"),
        pl.col("umi1_umi2").str.split("_").arr.get(0).alias("umi1_reorder"),
        pl.col("umi1_umi2").str.split("_").arr.get(1).alias("umi2_reorder"),
    ]
    df = df.with_columns(expr)
    df = df.select(
        [
            "edge_id",
            "edge_umi",
            "bead1_reorder",
            "umi1_reorder",
            "bead2_reorder",
            "umi2_reorder",
        ]
    ).rename(
        {
            "bead1_reorder": "bead1",
            "umi1_reorder": "umi1",
            "bead2_reorder": "bead2",
            "umi2_reorder": "umi2",
        }
    )
    return df


def filter_edgelist(edgelist, allowlist):
    """
    Filter edgelist by inner joins & filter bead1 and bead2.
    """
    out = edgelist.join(
        allowlist.drop("rank"), left_on="bead1", right_on="droplet", how="inner"
    ).join(allowlist.drop("rank"), left_on="bead2", right_on="droplet", how="inner")
    out = out.collect()
    return out


if __name__ == "__main__":
    main()
