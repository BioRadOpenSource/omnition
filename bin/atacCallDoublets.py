#!/usr/bin/env python3

import click
import polars as pl
import os
import numpy as np
from rpy2.robjects import vectors
from rpy2.robjects.packages import STAP
import math
import rle
import re
from collections import Counter
import gzip
import shutil


@click.command()
@click.option(
    "--csv_dir",
    "-d",
    help="directory of csv (previously .rds) files",
)
@click.option(
    "--n_bc",
    "-c",
    help="File for the number of reads supporting each barcode.",
)
@click.option(
    "--hq_bc",
    "-q",
    help="file for the HQ barcodes that were nominated.",
)
@click.option(
    "--params_file",
    "-p",
    help="deconvolution params file.",
)
@click.option(
    "--min_jaccard_frag",
    "-m",
    default=0.0,
    help="Must be provided as int.",
)
@click.option(
    "--name",
    "-n",
    help="name prefix for file naming convention.",
)
@click.option(
    "--one_to_one",
    "-o",
    default="False",
    help="Boolean arguement for keeping bead : drop conversion 1 to 1.",
)
@click.option(
    "--catac_assay",
    "-b",
    help="Boolean value",
)
@click.option(
    "--ti_len",
    "-l",
    help="Must be provided as int.",
)
def main(
    csv_dir,
    n_bc,
    hq_bc,
    params_file,
    min_jaccard_frag,
    name,
    one_to_one,
    catac_assay,
    ti_len,
):
    # get the csv files from the the dir
    overlap_df = load_overlap_df(csv_dir)

    # Only consider merging when Tn5 is the same
    catac_assay = eval(catac_assay)
    if catac_assay:
        overlap_df = substr_right(overlap_df, "barc1", int(ti_len))
        overlap_df = substr_right(overlap_df, "barc2", int(ti_len))
        overlap_df = overlap_df.filter(
            pl.col("match_barc1") == pl.col("match_barc2")
        ).select([pl.col("barc1"), pl.col("barc2"), pl.col("n_both")])

    # Import number of barcodes
    allowlist_bc = read_hq_bc_file(hq_bc)
    quantification_df = read_n_bc_file(n_bc, allowlist_bc["bc"].to_list())
    count_df = create_count_df(quantification_df.lazy())
    implicated_df = create_implicated_df(overlap_df, count_df)

    min_jaccard_frag = float(min_jaccard_frag)
    # Call knee if we need to
    if min_jaccard_frag == 0.0:
        print("Computing jaccard index for bead merging via a knee call--")
        jaccard_results = get_density_threshold(
            implicated_df.select([pl.col("jaccard_frag")]).head(1000000),
            "jaccard",
            logTransform=True,
        )
        min_jaccard_frag = jaccard_results[0]
        jaccard_called_frag = jaccard_results[1]
    else:
        jaccard_called_frag = min_jaccard_frag

    # Appending column stating whether merged or n
    implicated_df = implicated_df.with_column(
        pl.col("jaccard_frag").apply(lambda x: x > min_jaccard_frag).alias("merged")
    )
    print("PROGRESS: Merging bead barcodes into droplet barcodes.")

    # Filter out barcodes that will not be merged
    barcode_filtered_df = implicated_df.filter(pl.col("merged"))

    # Prepare barcode translate df
    barcode_translate_df = quantification_df.with_column(
        pl.lit("").alias("droplet_barcode")
    )

    # merge those beads
    barcode_translate_df = group_beads(
        quantification_df,
        barcode_filtered_df,
        barcode_translate_df,
        eval(one_to_one),
        name,
        ti_len,
        catac_assay,
    )

    # Output ------------------------------------------------------------------
    with open(params_file, "a") as f:
        f.write("jaccard_threshold_nosafety," + str(jaccard_called_frag) + "\n")
        # Append to deconvolution parameters
        f.write("jaccard_threshold," + str(min_jaccard_frag) + "\n")

    # update for outputs
    implicated_df = implicated_df.select(
        [
            pl.all().exclude("merged"),
            pl.col("merged")
            .cast(str)
            .apply(lambda x: "TRUE" if x == "true" else "FALSE")
            .alias("merged"),
        ]
    )

    #  Write output files based on presence/absence of TIs
    if not catac_assay:
        implicated_file = name + ".implicatedBarcodes.csv"
        translate_file = name + ".barcodeTranslate.tsv"

        # Export the implicated barcode table df
        implicated_df.write_csv(implicated_file, sep=",")

        # Compress output file
        with open(implicated_file, "rb") as f_in:
            with gzip.open(implicated_file + ".gz", "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        # Export the barcode translate df
        barcode_translate_df.select(
            [pl.col("bead_barcode"), pl.col("droplet_barcode")]
        ).write_csv(translate_file, sep="\t", has_header=False)
    else:
        fastq_tis = [
            f.split(".")[0]
            for f in os.listdir(csv_dir)
            if f.endswith(".deconvolutionParams.orig.csv")
        ]
        for fastq_ti in fastq_tis:
            ti = fastq_ti[fastq_ti.find("-") + 1:]
            # Export the implicated barcodes on a ti level
            implicated_dfti = implicated_df.filter(pl.col("barc1").str.contains(ti))
            implicated_dfti.write_csv(fastq_ti + ".implicatedBarcodes.csv", sep=",")

            # Compress output file
            with open(fastq_ti + ".implicatedBarcodes.csv", "rb") as f_in:
                with gzip.open(fastq_ti + ".implicatedBarcodes.csv.gz", "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)

            # Export the implicated barcodes
            barcode_translate_final_df = barcode_translate_df.select(
                [pl.col("bead_barcode"), pl.col("droplet_barcode")]
            )
            barcode_translate_ti_final_df = barcode_translate_final_df.filter(
                pl.col("bead_barcode").str.contains(ti)
            )
            barcode_translate_ti_final_df.write_csv(
                fastq_ti + ".barcodeTranslate.tsv", sep="\t", has_header=False
            )


def load_overlap_df(csv_dir):
    """
    Function to load the csv files from the current csv directory and create overlap_df.
    """
    files = [
        pl.read_csv(
            csv_dir + f, dtypes={"barc1": str, "barc2": str, "n_both": pl.Int64}
        )
        for f in os.listdir(csv_dir)
        if f.endswith("_overlapCount.csv.gz")
    ]
    df_vertical_concat = pl.concat(
        files,
        how="vertical",
    )

    return df_vertical_concat


def substr_right(input_df, column, ti_len):
    """
    Substring the right porton of a string specified by the ti_len parameter
    """
    result = input_df.select(
        [pl.all(), pl.col(column).apply(lambda x: x[-ti_len:]).alias("match_" + column)]
    )
    return result


def read_hq_bc_file(hq_bc):
    """
    Function to read in the high quality barcode file.
    """
    df = pl.read_csv(hq_bc, has_header=False, new_columns=["bc"])
    return df


def read_n_bc_file(n_bc, allowlist_bc):
    """
    Function to read in the number of barcodes file.
    """
    quantification_df = pl.read_csv(
        n_bc, has_header=False, new_columns=["bead_barcode", "count"], sep=","
    )
    quantification_df = quantification_df.filter(
        pl.col("bead_barcode").is_in(allowlist_bc)
    )
    quantification_df = quantification_df.sort("count", reverse=True)
    return quantification_df


def create_count_df(quantification_df):
    """
    Function to create the count df from the quantification_df
    """
    count_df = quantification_df.select(
        [pl.col("bead_barcode"), pl.col("count").apply(lambda x: x * 2)]
    )
    return count_df


def create_implicated_df(overlap_df, count_df):
    """
    Function to create the implicated df.
    """
    implicated_df = overlap_df.groupby(["barc1", "barc2"]).agg(
        pl.col("n_both").sum().alias("N_both")
    )
    implicated_df = implicated_df.join(
        count_df, left_on="barc1", right_on="bead_barcode"
    )
    implicated_df.columns = ["barc1", "barc2", "N_both", "N_barc1"]
    implicated_df = implicated_df.join(
        count_df, left_on="barc2", right_on="bead_barcode"
    )
    implicated_df.columns = ["barc1", "barc2", "N_both", "N_barc1", "N_barc2"]
    implicated_df = implicated_df.with_column(
        (
            pl.col("N_both")
            / (pl.col("N_barc1") + pl.col("N_barc2") - pl.col("N_both") + 0.05)
        ).alias("jaccard_frag")
    )
    implicated_df = implicated_df.select(
        [
            pl.all().exclude("jaccard_frag"),
            pl.col("jaccard_frag").apply(lambda x: round(x, 5)),
        ]
    )
    implicated_df = implicated_df.filter(
        pl.col("jaccard_frag") > 0
    )  # Remove pairs with score of 0
    implicated_df = implicated_df.sort(
        "jaccard_frag", reverse=True
    )  # Arrange from highest to lowest Jaccard index
    return implicated_df


def get_density_threshold(df, call_type, logTransform=True):
    """
    Function for calling jaccard knee based on local minima
    """
    # Initialize using some reasonable value and filter anything below
    threshold = get_mode(df) * 0.001
    filtered_df = df.filter(pl.col("jaccard_frag") > threshold)

    # Parameterize the log transformation to work with non-count data
    # May get confusing eventually but for now, seemingly a decent hack
    if logTransform:
        filtered_df = filtered_df.select(
            [
                pl.all(),
                pl.col("jaccard_frag")
                .apply(lambda x: round(math.log10(x), 6))
                .alias("log_counts"),
            ]
        )
    else:
        filtered_df = filtered_df.select(
            [pl.all(), pl.col("jaccard_frag").alias("log_counts")]
        )

    # Calculate the density using a gaussian kernel
    filtered_df_vec = vectors.FloatVector(filtered_df["log_counts"].to_list())
    v_dens_string = """vector_density_calc <- function(filtered_list) {
                    xx_values <- 10000
                    vector_density <-
                    density(
                        filtered_list,
                        bw = 0.1,
                        kernel = "gaussian",
                        n = xx_values,
                        from = min(filtered_list),
                        to = max(filtered_list)
                    )
                return(vector_density)
            }
        """
    r_pkg = STAP(v_dens_string, "r_pkg")

    vector_density = r_pkg.vector_density_calc(filtered_df_vec)
    vector_density_dict = dict(zip(vector_density.names, list(vector_density)))
    vector_density_y_list = list(vector_density_dict["y"])
    local_mins = get_local_minima(vector_density_y_list)

    # If a minima was called at the very start or end of the distribution, remove it
    local_mins = [
        val for val in local_mins if (val != 1) and (val != len(vector_density_y_list))
    ]

    local_min = find_local_min_in_list(
        local_mins[::-1], logTransform, filtered_df_vec, vector_density_dict
    )
    if local_min != 0:
        if logTransform:
            threshold = 10 ** vector_density_dict["x"][local_min]
        else:
            threshold = vector_density_dict["x"][local_min]
        print("Setting knee threshold to: ", threshold)
    else:
        print("No reliable knee found-- setting threshold to 0")
        threshold = 0
        local_min = 1
        local_mins = 1

    safety = 0
    # Safe guard for Jaccard Index failure
    if call_type == "jaccard" and (threshold > 0.5 or threshold < 0.000001):
        print("No reliable knee found-- setting threshold to 0.005")
        safety = 0.005

    # Safe guard for knee counts failure
    if call_type == "bead" and (threshold > 100000 or threshold < 100):
        print("No reliable knee found-- setting threshold to 500")
        safety = 500

    # Safety is with the guard rails; threshold is what the knee calls
    if not safety > 0:
        safety = threshold

    return (safety, threshold)


def get_mode(df):
    """
    Function to determine mode of an input df column
    """
    x = df["jaccard_frag"].to_list()
    ux = df["jaccard_frag"].unique().to_list()
    match = [ux.index(ind) if ind in ux else None for ind in x]
    tabulate = (np.bincount(match)).tolist()
    index = tabulate.index(max(tabulate))
    return ux[index]


def get_local_minima(x):
    """
    Get the local min
    """
    # Find the indices where we go from negative diff's to positive
    # if this returns less then 0 then we are descending
    y = np.diff([math.inf] + x) < 0

    # Find number of downward and upward steps and identify index of inflection point
    y = np.cumsum(rle.encode(list(y))[1])

    # If we only get three values, then this becomes a problem
    # Getting TRUE,FALSE,TRUE will keep the extreme values and remove the true minimum
    if len(y) > 3:
        y = [int(y[int(x)]) for x in np.arange(start=0, stop=len(y), step=2)]
    else:
        y = y.tolist()

    # Seth's modification for removing duplicated elements at the beginning
    if x[0] == x[1]:
        y = y[-1]
        # handling edge case
        if isinstance(y, int):
            y = [y]

    return y


def find_local_min_in_list(x, logTransform, filtered_df_vec, vector_density_dict):
    """
    Creating a new function to perform the find
    """
    # Make sure that the selected min includes at least 20% of barcodes
    # and that the difference between the min and the max differences
    # by some appreciable amount
    # both in terms of absolute difference (changed 0.5
    # to 0.05) AND relative difference

    if logTransform:
        abs_difference = 0.5
    else:
        abs_difference = 0.05
    xx_values = 10000

    found = 0
    for item in x:
        if item >= (0.2 * xx_values) and (
            (max(list(filtered_df_vec)) - vector_density_dict["x"][item])
            > abs_difference
            or (vector_density_dict["x"][item] < max(list(filtered_df_vec)) / 2)
        ):
            found = item  # return the first time this criteria is met
            break
    return found


def group_beads(
    quantification_df,
    barcode_filtered_df,
    barcode_translate_df,
    one_to_one,
    name,
    ti_len,
    catac_assay,
):
    """
    Function to group the beads together and create the droplets
    """
    # Guess at how many places will be needed for barcode index
    # when creating droplet barcode names
    # E.g., "<SAMPLE>_BC25_N03" requires 2 places for "BC25"
    droplet_barcode_places = math.ceil((math.log10(quantification_df.shape[0])))

    # Initializing index
    idx = 1

    # Loop through each row of the barcode quantification for beads on the allowlist
    while quantification_df.shape[0] > 0:
        # Pull first barcode in df
        barcode = quantification_df[0, 0]
        # Initializing merge list
        barcode_merge = [barcode]

        # If one-to-one is specified, don't merge barcodes
        # Otherwise, create list of barcodes to merge
        if not one_to_one:
            # Find barcodes that have overlapping fragments with barcode
            related_barcodes1 = barcode_filtered_df.filter(pl.col("barc1") == barcode)
            related_barcodes2 = barcode_filtered_df.filter(pl.col("barc2") == barcode)
            barcode_filtered_df = barcode_filtered_df.filter(
                pl.col("barc1") != barcode
            )  # to drop the ones where it was true
            barcode_filtered_df = barcode_filtered_df.filter(
                pl.col("barc2") != barcode
            )  # to drop the ones where it was true
            # Create list of barcodes to merge based on similarity
            if len(related_barcodes1) > 0:
                barcode_merge.extend(related_barcodes1["barc2"].to_list())
            if len(related_barcodes2) > 0:
                barcode_merge.extend(related_barcodes2["barc1"].to_list())

        # Format droplet barcode name depending on presence/absence of TIs
        if not catac_assay:
            droplet_name = name + "_BC" + f"{idx:0>{droplet_barcode_places}}" + "_N"
        else:
            droplet_name = (
                name
                + "_Tn5-"
                + substr_right(
                    pl.DataFrame(barcode_merge, columns=["bc"]), "bc", int(ti_len)
                )[0, 1]
                + "_BC"
                + f"{idx:0>{droplet_barcode_places}}"
                + "_N"
            )

        # Get the implicated barcodes droplet barcodes (if they exist)
        old_droplets = barcode_translate_df.filter(
            pl.col("bead_barcode").is_in(barcode_merge)
        )["droplet_barcode"].to_list()
        # want to only count barcodes that are already a part of droplets
        # so remove empty strings from list
        non_empty_old_droplets = Counter([item for item in old_droplets if item])
        droplet_size = 0
        intersection_size = 0
        for key, value in non_empty_old_droplets.items():
            droplet_size += int(re.search(r"N(\w+)", key[len(name + "_BC"):]).group(1))
            intersection_size += value
        # subtract from the droplet size the number of bc that are already in it,
        # add the new number of implicated bc
        droplet_size = droplet_size - intersection_size + len(barcode_merge)
        barcode_translate_df = barcode_translate_df.with_column(
            pl.when(
                pl.col("droplet_barcode").is_in(list(non_empty_old_droplets.keys()))
            )
            .then(droplet_name + "%02d" % droplet_size)
            .when(
                (pl.col("bead_barcode").is_in(barcode_merge))
                & (pl.col("droplet_barcode") == "")
            )
            .then(droplet_name + "%02d" % droplet_size)
            .otherwise(pl.col("droplet_barcode"))
            .alias("droplet_barcode")
        )

        # Remove barcodes that have been merged
        quantification_df = quantification_df.filter(
            ~pl.col("bead_barcode").is_in(barcode_merge)
        )
        # Increment index for next iteration
        idx += 1

    return barcode_translate_df


if __name__ == "__main__":
    main()
