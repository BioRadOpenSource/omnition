#!/usr/bin/env python3

import gzip
import statistics
import click
import csv
from typing import Tuple


@click.command()
@click.option(
    "--tss-window-size",
    "-w",
    default=4001,
    required=False,
    show_default=True,
    type=int,
    help="The value to use for tss window size, currently interpreted "
    "as the total window size, including central TSS base",
)
@click.option(
    "--tss-matrix",
    "-m",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="The .tss_data_matrix.gz file that is going to have its \
    enrichment score calculated.",
)
def main(tss_window_size, tss_matrix):
    """
    Calculate the TSS enrichment for each TSS window in the provided
    tss_data_matrix.gz file. The TSS enrichment is defined as the
    average coverage across all TSS windows, normalized by the average
    coverage in the first and last 100 bases of the TSS window.
    The TSS enrichment score is the maximum value of the normalized
    coverage across all bases in the TSS window.

    Args:
        tss_window_size (int): The total size of the TSS window, including
            the central TSS base. Default is 4001.
        tss_matrix (str): Path to the .tss_data_matrix.gz file containing
            TSS coverage data.
    """
    count, rolling_total = calculate_total_tss_coverage(
        tss_matrix, tss_window_size
    )

    norm, tss_enrichment_score = normalize_tss_coverage(
        rolling_total, count, tss_window_size
    )

    sample_name = tss_matrix.replace(".tss_data_matrix.gz", "")
    with open("{}.tss_enrichment.csv".format(sample_name), "w") as f:
        csv_out = csv.writer(f)
        csv_out.writerows([norm[index]] for index in range(0, len(norm)))

    with open("{}.tss_metric.csv".format(sample_name), "w") as f:
        f.write("sample,process,metric,value\n")
        f.write(f"{sample_name},tss_enrichment,tss_score,{tss_enrichment_score}\n")


def calculate_total_tss_coverage(tss_matrix: str, tss_window_size: int) -> Tuple[int, list]:
    """
    Calculate the total TSS coverage from the provided TSS matrix.

    Args:
        tss_matrix (str): Path to the .tss_data_matrix.gz file.
        tss_window_size (int): The total size of the TSS window, including
            the central TSS base.

    Returns:
        count (int): The total number of TSS windows with coverage.
        rolling_total (list): A list containing the total coverage at each
            relative base pair position around the TSS.
    """
    # create a vector initialized with tss window sized numbers of 0s
    rolling_total = [0] * (
        tss_window_size
    )

    # iterate over TSS windows; for those that have coverage,
    # add to the rolling total of coverage at each base
    count = 0
    with gzip.open(tss_matrix, "rt") as f:
        for line in f:
            vals = [
                int(i) for i in line.strip().split("\t")
            ]  # strip newlines, split by tab, convert to int. is a
            #  vector of coverage values for bp around tss
            if sum(vals) > 0:  # check if there's any coverage in the tss window
                count += 1  # total number of gene tss with any coverage
                # add new tss coverage data to all previous data to get
                for i, val in enumerate(vals):
                    rolling_total[i] += val
                # global coverage values at each relative bp
    return count, rolling_total


def normalize_tss_coverage(rolling_total: list, count: int, tss_window_size: int) -> Tuple[list, float]:
    """
    Normalize the TSS coverage by the count of TSS windows.

    Args:
        rolling_total (list): A list containing the total coverage at each
            relative base pair position around the TSS.
        count (int): The total number of TSS windows with coverage.
        tss_window_size (int): The total size of the TSS window, including

    Returns:
        norm (list): A normalized list of coverage values.
        tss_enrichment_score (float): The maximum normalized coverage value
    """
    if count > 0:
        # compute the average coverage for each base
        div = [x / count for x in rolling_total]

        # caluclate noise, defined as the average for the first and last 100 bases
        avg_noise = statistics.mean(div[:100] + div[-100:])

        # normalize by average noise
        norm = [x / avg_noise for x in div]

    else:
        norm = [0] * (tss_window_size)
    # Calculate the tss enrichment score by taking the maximum value
    tss_enrichment_score = round(max(norm), 1)
    return norm, tss_enrichment_score


if __name__ == "__main__":
    main()
