#!/usr/bin/env python3

import click
import polars as pl
import os
from glob import glob


@click.command()
@click.option(
    "--input-dir",
    "-i",
    required=True,
    type=click.Path(exists=True),
    help="Directory containing the input files.",
)
@click.option(
    "--sample-id",
    "-s",
    required=True,
    type=str,
    help="Sample ID to be used in the output file.",
)
@click.option(
    "--read-counts-suffix",
    "-f",
    default="_read_counts.csv",
    show_default=True,
    type=str,
    help="Suffix for the read count files.",
)
def main(input_dir, sample_id, read_counts_suffix):
    """
    Merge the read counts from the reannotate fragments and
    reannotate bam output files into a single CSV file.

    Args:
        input_dir (str): Directory containing the input files.
        sample_id (str): Sample ID to be used in the output file.
        read_counts_suffix (str): Suffix for the read count files.
    """
    all_counts = collect_counts(input_dir, read_counts_suffix)
    summarized_counts = summarize_counts(all_counts, sample_id)
    summarized_counts.write_csv(
        f"{input_dir}/{sample_id}_merged_reannotate_read_counts.csv"
    )


def collect_counts(input_dir: str, read_counts_suffix: str) -> pl.DataFrame:
    """
    Collect read counts from fragments and BAM files.

    Args:
        input_dir (str): Directory containing the input files.
        read_counts_suffix (str): Suffix for the read count files.

    Returns:
        pl.DataFrame: DataFrame containing all read counts.
    """
    # list files in the input directory that have the specified suffixes
    count_files = glob(os.path.join(input_dir, f"*{read_counts_suffix}"))
    count_dfs = [
        pl.read_csv(file) for file in count_files
    ]

    return pl.concat(count_dfs, how="vertical")


def sum_counts(df: pl.DataFrame, process: str, metric: str) -> pl.DataFrame:
    """
    Sum the counts for a specific process and metric.

    Args:
        df (pl.DataFrame): DataFrame containing read counts.
        process (str): Process name to filter by.
        metric (str): Metric name to filter by.

    Returns:
        pl.DataFrame: DataFrame with summed counts.
    """
    return df.filter(
        (pl.col("process").str.starts_with(process)) & (pl.col("metric") == metric)
    )["value"].sum()


def summarize_counts(all_counts: pl.DataFrame, sample_id: str) -> pl.DataFrame:
    """
    Summarize read counts by sample ID.

    Args:
        all_counts (pl.DataFrame): DataFrame containing all read counts.
        sample_id (str): Sample ID to be used in the output file.

    Returns:
        pl.DataFrame: DataFrame with summarized read counts.
    """
    processes = ["assemble_fragments", "reannotate_fragments", "reannotate_bam"]
    output_metrics = {
        "sample": [],
        "process": [],
        "metric": [],
        "value": []
    }
    for process in processes:
        for metric in ["input", "output"]:
            output_metrics["sample"].append(sample_id)
            output_metrics["process"].append(process)
            output_metrics["metric"].append(metric)
            value = sum_counts(all_counts, process, metric)
            output_metrics["value"].append(value)

    return pl.DataFrame(output_metrics)


if __name__ == "__main__":
    main()
