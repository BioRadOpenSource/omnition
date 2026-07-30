#!/usr/bin/env python3

import click
import os
import statistics
from typing import Dict, List, Tuple
import polars as pl


@click.command()
@click.option(
    "--sample-id",
    "-s",
    required=True,
    type=str,
    help="The sample ID of the compiled sample."
)
@click.option(
    "--input-directory",
    "-i",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Input directory containing deconvolution files.",
)
@click.option(
    "--qc-stats-suffix",
    "-q",
    default=".QCstats.csv",
    required=False,
    show_default=True,
    type=str,
    help="Suffix for the QC stats files."
)
@click.option(
    "--fastq-ti-counts-file",
    "-f",
    default="fastqTIreadcounts.csv",
    required=False,
    show_default=True,
    type=str,
    help="Name of the fastqTIreadcounts file."
)
def main(sample_id, input_directory, qc_stats_suffix, fastq_ti_counts_file):
    """
    Collects metrics from the QC stats files and fastqTIreadcounts.csv files
    for the deconvolution process and outputs a summary CSV file.

    Args:
        sample_id (str): The sample ID for which to collect metrics.
        input_directory (str): The directory containing the QC stats files and fastqTIreadcounts.csv.
        qc_stats_suffix (str): The suffix for the QC stats files.
        fastq_ti_counts_file (str): The name of the fastqTIreadcounts file.
    """
    # Initialize the metrics dictionary and sample metrics
    metrics_dict = {
        "sample": [],
        "process": [],
        "metric": [],
        "value": [],
    }
    qc_stats = []
    # Gather the fastq+TI stats
    ti_read_counts = read_ti_counts(sample_id, input_directory, fastq_ti_counts_file)
    names = find_fastq_ti_names(input_directory, qc_stats_suffix)
    for i in range(len(names)):
        metrics_dict, qc_stats_df = collect_fastq_ti_metrics(
            metrics_dict, names[i], input_directory, ti_read_counts[names[i]], qc_stats_suffix
        )
        qc_stats.append(qc_stats_df)

    # Make stats for the sample as a whole
    metrics_dict = collect_sample_metrics(
        sample_id, metrics_dict, ti_read_counts, qc_stats
    )

    metrics_df = pl.DataFrame(metrics_dict, strict=False)
    metrics_df.write_csv(
        input_directory + sample_id + ".beadFiltSummary.csv"
    )


def find_fastq_ti_names(input_dir: str, qc_stats_suffix: str) -> List[str]:
    """
    Makes a list of all fastq+TI pair names from the input directory.

    Args:
        input_dir (str): The directory containing the QC stats files.
        qc_stats_suffix (str): The suffix for the QC stats files.

    Returns:
        list: A sorted list of fastq+TI pair names.
    """
    names = []
    dirfiles = os.listdir(input_dir)
    for i in range(len(dirfiles)):
        if dirfiles[i].find(qc_stats_suffix) != -1:
            names.append(dirfiles[i].replace(qc_stats_suffix, ""))
    names.sort()

    return names


def read_ti_counts(sample: str, input_dir: str, fastq_ti_counts: str) -> Dict[str, int]:
    """
    Reads the fastq_ti_counts file and returns a dictionary
    mapping fastq+TI pairs to their read counts for the specified sample.

    Args:
        sample (str): The sample ID to filter the read counts.
        input_dir (str): The directory containing the fastq_ti_counts file.
        fastq_ti_counts (str): The name of the fastq_ti_counts file.

    Returns:
        dict: A dictionary mapping fastq+TI pairs to their read counts.
    """
    fastq_ti_counts_df = pl.read_csv(input_dir + fastq_ti_counts)
    fastq_ti_counts_dict = fastq_ti_counts_df.with_columns(
        (pl.col("fastq") + "-" + pl.col("sequence")).alias("fastqTI")
    ).filter(pl.col("sample") == sample).select("fastqTI", "count").rows_by_key(key=["fastqTI"])

    for fastq_ti, count in fastq_ti_counts_dict.items():
        fastq_ti_counts_dict[fastq_ti] = count[0][0]
    return fastq_ti_counts_dict


def collect_fastq_ti_metrics(metrics_dict: Dict[str, List], name: str, input_directory: str, reads: int, qc_stats_suffix: str) -> Tuple[Dict[str, List], pl.DataFrame]:
    """
    Collects metrics for a single fastq+TI pair and updates the metrics dictionary.

    Args:
        metrics_dict (dict): The dictionary to store collected metrics.
        name (str): The name of the fastq+TI pair.
        input_directory (str): The directory containing the QC stats files.
        reads (int): The number of reads for the fastq+TI pair.
        qc_stats_suffix (str): The suffix for the QC stats files.

    Returns:
        tuple: A tuple containing the updated metrics dictionary and the QC stats DataFrame.
    """
    qc_stats = pl.read_csv(input_directory + name + qc_stats_suffix, infer_schema_length=None)
    cells = qc_stats.shape[0]
    metrics_dict["sample"].append(name)
    metrics_dict["process"].append("deconvolution")
    metrics_dict["metric"].append("mean_valid_read_pairs_per_cell")
    metrics_dict["value"].append(str(round(reads / cells, 1)))

    return metrics_dict, qc_stats


def collect_sample_metrics(sample_id: str, metrics_dict: Dict[str, List], reads_per_index: Dict[str, int], qc_stats: List[pl.DataFrame]) -> Dict[str, List]:
    """
    Collects and calculates sample-level metrics from the provided data.

    Args:
        sample_id (str): The sample ID for which to collect metrics.
        metrics_dict (dict): The dictionary to store collected metrics.
        reads_per_index (dict): A dictionary mapping fastq+TI pairs to their read counts.
        qc_stats (list): A list of DataFrames containing QC stats for each fastq+TI pair.

    Returns:
        dict: The updated metrics dictionary with sample-level metrics.
    """
    sample_metrics = dict()
    sample_metrics["total_reads"] = sum(reads_per_index.values())
    sample_metrics["median_reads_per_index"] = round(
        statistics.median(reads_per_index.values()), 1
    )
    all_qc_stats = pl.concat(qc_stats, how="vertical")
    sample_metrics["total_cells"] = all_qc_stats.shape[0]
    sample_metrics["mean_valid_read_pairs_per_cell"] = round(
        sample_metrics["total_reads"] / sample_metrics["total_cells"], 1
    )
    sample_metrics["median_unique_nuclear_fragments_per_cell"] = all_qc_stats.select(
        pl.col("uniqueNuclearFrags")
    ).median().item()
    sample_metrics["median_estimated_library_size_per_cell"] = all_qc_stats.select(
        pl.col("librarySize")
    ).median().item()

    for metric, value in sample_metrics.items():
        metrics_dict["sample"].append(sample_id)
        metrics_dict["process"].append("summary")
        metrics_dict["metric"].append(metric)
        metrics_dict["value"].append(value)
    return metrics_dict


if __name__ == "__main__":
    main()
