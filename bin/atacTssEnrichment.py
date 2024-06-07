#!/usr/bin/env python3

import gzip
import statistics
import click
import csv


@click.command()
@click.option(
    "--tss-window-size",
    "-w",
    default=4001,
    help="The value to use for tss window size, currently interpreted "
    "as the total window size, including central TSS base",
)
@click.option(
    "--tss-matrix",
    "-m",
    help="The .tss_data_matrix.gz file that is going to have its \
    enrichment score calculated.",
)
def main(tss_window_size, tss_matrix):
    rolling_total = [0] * (
        tss_window_size
    )  # where 4001 = the tss window size + 1 base for the actual TSS,
    # creates vector initialized with tss window sized numbers of 0s

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
                rolling_total = [
                    a + b for a, b in zip(rolling_total, vals)
                ]  # add new tss coverage data to all previous data to get
                # global coverage values at each relative bp

    if count > 0:
        # compute the average coverage for each base
        div = [x / count for x in rolling_total]

        # caluclate noise, defined as the average for the first and last 100 bases
        avg_noise = statistics.mean(div[:100] + div[-100:])

        # normalize by average noise
        norm = [x / avg_noise for x in div]

    else:
        norm = [0] * (tss_window_size)

    sample_name = tss_matrix.replace(".tss_data_matrix.gz", "")
    with open("{}.tss_enrichment.csv".format(sample_name), "w") as f:
        csv_out = csv.writer(f)
        csv_out.writerows([norm[index]] for index in range(0, len(norm)))


if __name__ == "__main__":
    main()
