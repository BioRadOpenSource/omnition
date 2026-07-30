#!/usr/bin/env python3

import click
import polars as pl


@click.command()
@click.option(
    "--raw_edges",
    "-r",
    type=click.Path(exists=True),
    required=True,
    help="List of raw input edges.",
)
@click.option(
    "--corrected_edges",
    "-c",
    type=click.Path(exists=True),
    required=True,
    help="List of corrected edges.",
)
@click.option(
    "--barcodes",
    "-b",
    type=click.Path(exists=True),
    required=True,
    help="Barcode allowlist.",
)
@click.option(
    "--sample_id",
    "-s",
    type=str,
    required=False,
    help="Sample ID for formatting output files.",
)
def main(raw_edges, corrected_edges, barcodes, sample_id):
    """
    Collect the metadata for raw edges,
    corrected edges, and above knee edges

    Args:
        raw_edges (str): Path to raw edges file.
        corrected_edges (str): Path to corrected edges file.
        barcodes (str): Path to barcode allowlist.
        sample_id (str): Sample ID for formatting output files.
    """
    # read allowlist
    allowlist = read_allowlist(barcodes)

    # read edgelists
    raw_edgelist = read_edgelist(raw_edges)
    corrected_edgelist = read_edgelist(corrected_edges)

    # Filter above knee edges
    filtered_edgelist = filter_edgelist(corrected_edgelist, allowlist.lazy())

    # count edges
    count_dict = {}
    count_dict = count_edge_types(raw_edgelist.collect(), "raw", count_dict)
    count_dict = count_edge_types(corrected_edgelist.collect(), "corrected", count_dict)
    count_dict = count_edge_types(filtered_edgelist, "above_knee", count_dict)

    # write output
    final_summary_df = create_final_summary(count_dict, sample_id)
    final_summary_df.write_csv(f"{sample_id}_edge_metadata.csv")


def read_allowlist(allowlist_file: str) -> pl.DataFrame:
    """
    Read allowlist file into a dataframe

    Args:
        allowlist_file (str): The path to the allowlist file.

    Returns:
        pl.DataFrame: The dataframe with allowlist
    """
    allowlist = pl.read_csv(
        allowlist_file,
        separator="\t",
        has_header=False,
        new_columns=["bead"],
        schema_overrides={"bead": pl.Utf8}
    ).with_row_index(name="rank", offset=1)
    return allowlist


def read_edgelist(edge_file):
    """
    Read edgelist file into a dataframe, removing duplicate rows.
    Order such that bead1 and umi1 are alphabetically first.
    """
    cols = [0, 1, 2, 3]
    new_cols = ["bead1", "umi1", "bead2", "umi2"]
    set_schema = {"bead1": pl.Utf8, "umi1": pl.Utf8, "bead2": pl.Utf8, "umi2": pl.Utf8}
    df = (
        pl.read_csv(
            edge_file,
            has_header=False,
            separator="\t",
            columns=cols,
            new_columns=new_cols,
            schema_overrides=set_schema,
        ).lazy()
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
        pl.col("bead1_bead2").str.split("_").list.get(0).alias("bead1_reorder"),
        pl.col("bead1_bead2").str.split("_").list.get(1).alias("bead2_reorder"),
        pl.col("umi1_umi2").str.split("_").list.get(0).alias("umi1_reorder"),
        pl.col("umi1_umi2").str.split("_").list.get(1).alias("umi2_reorder"),
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


def count_edge_types(edgelist, edge_type, count_dict):
    """
    Count the total, unique, redundant, and abmiguous edges.
    Total = all edges including duplicates
    Unique = unique edges
    Ambiguous = edges with N bases
    Final = Unique bead pair edges with no Ns
    Redundant = unique edges - ambiguous edges - final edges

    Args:
        edgelist (pl.DataFrame): Edgelist dataframe.
        edge_type (str): Type of edge.
        count_dict (dict): Dictionary to store counts.

    Returns:
        dict: Dictionary with updated counts
    """
    if edge_type == "raw":
        count_dict["total_edges"] = edgelist.shape[0]
    elif edge_type == "above_knee":
        count_dict[f"{edge_type}_edges"] = edgelist.shape[0]
    count_dict[f"{edge_type}_unique_edges"] = edgelist.n_unique()
    count_dict[f"{edge_type}_ambiguous_umis"] = edgelist.filter(
        (pl.col("edge_id").str.contains("N")
         | pl.col("edge_umi").str.contains("N"))
    ).n_unique()
    final_edges = edgelist.filter(
        (~pl.col("edge_id").str.contains("N"))
        & (~pl.col("edge_umi").str.contains("N"))
    ).select("edge_id").n_unique()
    count_dict[f"{edge_type}_redundant_edges"] = count_dict[f"{edge_type}_unique_edges"] - count_dict[f"{edge_type}_ambiguous_umis"] - final_edges

    return count_dict


def filter_edgelist(edgelist, allowlist):
    """
    Filter edgelist by inner joins & filter bead1 and bead2.
    """
    out = edgelist.join(
        allowlist.drop("rank"), left_on="bead1", right_on="bead", how="inner"
    ).join(allowlist.drop("rank"), left_on="bead2", right_on="bead", how="inner")
    out = out.collect()
    return out


def create_final_summary(count_dict, sample_id):
    """
    Create final metric summary data frame

    Args:
        count_dict (dict): Dictionary with counts.
        sample_id (str): Sample ID.

    Returns:
        pl.DataFrame: Final summary dataframe.
    """
    count_dict["percent_raw_unique_deconvolution_reads"] = round(
        count_dict["raw_unique_edges"] / count_dict["total_edges"] * 100, 2
    )
    count_dict["percent_corrected_unique_deconvolution_reads"] = round(
        count_dict["corrected_unique_edges"] / count_dict["total_edges"] * 100, 2
    )
    count_dict["percent_above_knee_unique_deconvolution_reads"] = round(
        count_dict["above_knee_unique_edges"] / count_dict["total_edges"] * 100, 2
    )
    final_metrics_df = pl.DataFrame(
        {
            "metric": [metric for metric in count_dict.keys()],
            "value": [count for count in count_dict.values()]
        },
        strict=False
    )
    final_metrics_df = final_metrics_df.with_columns(
        [
            pl.lit(sample_id).alias("sample"),
            pl.lit("edge_metadata").alias("process")
        ]
    ).select(["sample", "process", "metric", "value"])
    return final_metrics_df


if __name__ == "__main__":
    main()
