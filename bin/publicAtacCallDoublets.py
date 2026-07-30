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
from typing import Tuple
from collections import Counter
import gzip
import shutil


@click.command()
@click.option(
    "--csv_dir",
    "-d",
    required=True,
    type=click.Path(exists=True),
    help="directory of csv (previously .rds) files",
)
@click.option(
    "--n_bc",
    "-c",
    required=True,
    type=click.Path(exists=True),
    help="File for the number of reads supporting each barcode.",
)
@click.option(
    "--hq_bc",
    "-q",
    required=True,
    type=click.Path(exists=True),
    help="file for the HQ barcodes that were nominated.",
)
@click.option(
    "--input_params_suffix",
    "-s",
    type=str,
    required=False,
    show_default=True,
    help="Suffix of the input deconvolution params file.",
    default=".deconvolutionParams.orig.csv",
)
@click.option(
    "--min_jaccard_frag",
    "-m",
    required=True,
    type=float,
    help="Must be provided as int.",
)
@click.option(
    "--name",
    "-n",
    required=True,
    type=str,
    help="name prefix for file naming convention.",
)
@click.option(
    "--one_to_one",
    "-o",
    is_flag=True,
    required=False,
    type=bool,
    help="Boolean arguement for keeping bead : drop conversion 1 to 1.",
)
@click.option(
    "--catac_assay",
    "-b",
    is_flag=True,
    required=False,
    type=bool,
    help="Boolean value",
)
@click.option(
    "--ti_len",
    "-l",
    required=True,
    type=int,
    help="Must be provided as int.",
)
def main(
    csv_dir,
    n_bc,
    hq_bc,
    input_params_suffix,
    min_jaccard_frag,
    name,
    one_to_one,
    catac_assay,
    ti_len,
):
    # get the csv files from the the dir
    overlap_df = load_overlap_df(csv_dir)

    # Only consider merging when Tn5 is the same
    if catac_assay:
        overlap_df = substr_right(overlap_df, "barc1", int(ti_len))
        overlap_df = substr_right(overlap_df, "barc2", int(ti_len))
        overlap_df = overlap_df.filter(
            pl.col("match_barc1") == pl.col("match_barc2")
        ).select([pl.col("barc1"), pl.col("barc2"), pl.col("n_both")])

    # Import number of barcodes
    allowlist_bc = read_hq_bc_file(hq_bc)
    unfiltered_quant_df, quantification_df = read_n_bc_file(n_bc, allowlist_bc["bc"].to_list())
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
    implicated_df = implicated_df.with_columns(
        (pl.col("jaccard_frag") > min_jaccard_frag).alias("merged")
    )
    print("PROGRESS: Merging bead barcodes into droplet barcodes.")

    # Filter out barcodes that will not be merged
    barcode_filtered_df = implicated_df.filter(pl.col("merged"))

    # Prepare barcode translate df
    barcode_translate_df = quantification_df.with_columns(
        pl.lit("").alias("droplet_barcode")
    )

    # merge those beads
    barcode_translate_df = group_beads(
        quantification_df,
        barcode_filtered_df,
        barcode_translate_df,
        one_to_one,
        name,
        ti_len,
        catac_assay,
    )

    # update for outputs
    implicated_df = implicated_df.select(
        [
            pl.all().exclude("merged"),
            pl.col("merged")
            .cast(str)
            .map_elements(lambda x: "TRUE" if x == "true" else "FALSE", return_dtype=pl.String)
            .alias("merged"),
        ]
    )

    #  Write output files based on presence/absence of TIs
    if not catac_assay:
        implicated_file = name + ".implicatedBarcodes.csv"
        translate_file = name + ".barcodeTranslate.tsv"

        # Export the implicated barcode table df
        implicated_df.write_csv(implicated_file, separator=",")

        # Compress output file
        with open(implicated_file, "rb") as f_in:
            with gzip.open(implicated_file + ".gz", "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        # Export the barcode translate df
        barcode_translate_df.select(
            [pl.col("bead_barcode"), pl.col("droplet_barcode")]
        ).write_csv(translate_file, separator="\t", include_header=False)

        # Create the new deconvolution params file
        bead_threshold, updated_params_df = create_updated_deconvolution_params(
            jaccard_called_frag,
            min_jaccard_frag,
            name,
            input_params_suffix,
        )
        updated_params_df.write_csv(
            name + ".deconvolutionParams.csv",
            separator=",",
            include_header=False,
        )

        # Generate metrics summary output
        metrics_df = compile_deconvolution_metrics(
            name,
            unfiltered_quant_df,
            quantification_df,
            bead_threshold,
            implicated_df,
            min_jaccard_frag,
            barcode_translate_df
        )
    else:
        fastq_tis = [
            f.split(".")[0]
            for f in os.listdir(csv_dir)
            if f.endswith(".deconvolutionParams.orig.csv")
        ]
        metrics_ti_dfs = []
        for fastq_ti in fastq_tis:
            ti = fastq_ti[fastq_ti.find("-") + 1:]
            # Export the implicated barcodes on a ti level
            implicated_dfti = implicated_df.filter(pl.col("barc1").str.ends_with(ti))
            implicated_dfti.write_csv(fastq_ti + ".implicatedBarcodes.csv", separator=",")

            # Compress output file
            with open(fastq_ti + ".implicatedBarcodes.csv", "rb") as f_in:
                with gzip.open(fastq_ti + ".implicatedBarcodes.csv.gz", "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)

            # Export the implicated barcodes
            barcode_translate_final_df = barcode_translate_df.select(
                [pl.col("bead_barcode"), pl.col("droplet_barcode")]
            )
            barcode_translate_ti_final_df = barcode_translate_final_df.filter(
                pl.col("bead_barcode").str.ends_with(ti)
            )
            barcode_translate_ti_final_df.write_csv(
                fastq_ti + ".barcodeTranslate.tsv", separator="\t", include_header=False
            )

            # Create the new deconvolution params file
            bead_threshold, updated_params_df = create_updated_deconvolution_params(
                jaccard_called_frag,
                min_jaccard_frag,
                fastq_ti,
                input_params_suffix,
            )
            updated_params_df.write_csv(
                fastq_ti + ".deconvolutionParams.csv",
                separator=",",
                include_header=False,
            )

            # Generate metrics summary output
            unfiltered_quant_df_ti = unfiltered_quant_df.filter(
                pl.col("bead_barcode").str.ends_with(ti)
            )
            quantification_df_ti = quantification_df.filter(
                pl.col("bead_barcode").str.ends_with(ti)
            )
            metrics_ti_dfs.append(compile_deconvolution_metrics(
                fastq_ti,
                unfiltered_quant_df_ti,
                quantification_df_ti,
                bead_threshold,
                implicated_dfti,
                min_jaccard_frag,
                barcode_translate_ti_final_df
            ))
        metrics_df = pl.concat(metrics_ti_dfs, how="vertical")
    metrics_df.write_csv(
        name + ".deconvolution_metrics.csv",
        separator=",",
        include_header=True,
    )


def load_overlap_df(csv_dir: str) -> pl.DataFrame:
    """
    Function to load the csv files from the current csv directory and create overlap_df.

    Args:
        csv_dir (str): Directory containing the CSV files.
    Returns:
        pl.DataFrame: A DataFrame containing the concatenated data from all CSV files.
    """
    files = [
        pl.read_csv(
            csv_dir + f, schema_overrides={"barc1": str, "barc2": str, "n_both": pl.Int64}
        )
        for f in os.listdir(csv_dir)
        if f.endswith("_overlapCount.csv.gz")
    ]
    df_vertical_concat = pl.concat(
        files,
        how="vertical",
    )

    return df_vertical_concat


def substr_right(input_df: pl.DataFrame, column: str, ti_len: int) -> pl.DataFrame:
    """
    Substring the right porton of a string specified by the ti_len parameter

    Args:
        input_df (pl.DataFrame): Input DataFrame containing the column to be processed.
        column (str): Name of the column to be processed.
        ti_len (int): Length of the substring to extract from the right side of the string.
    Returns:
        pl.DataFrame: A DataFrame with an additional column containing the right substring of the specified column.
    """
    result = input_df.select(
        [pl.all(), pl.col(column).map_elements(lambda x: x[-ti_len:], return_dtype=pl.Utf8).alias("match_" + column)]
    )
    return result


def read_hq_bc_file(hq_bc: str) -> pl.DataFrame:
    """
    Function to read in the high quality barcode file.

    Args:
        hq_bc (str): Path to the high quality barcode file.
    Returns:
        pl.DataFrame: A DataFrame containing the high quality barcodes.
    """
    df = pl.read_csv(hq_bc, has_header=False, new_columns=["bc"])
    return df


def read_n_bc_file(n_bc: str, allowlist_bc: list) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """
    Function to read in the number of barcodes file.

    Args:
        n_bc (str): Path to the number of barcodes file.
        allowlist_bc (list): List of bead barcodes that are allowed.
    Returns:
        - unfiltered_quantification_df: DataFrame with all bead barcodes and their counts.
        - quantification_df: Filtered DataFrame with bead barcodes that are in the allowlist.
    """
    unfiltered_quantification_df = pl.read_csv(
        n_bc, has_header=False, new_columns=["bead_barcode", "count"], separator=","
    )
    quantification_df = unfiltered_quantification_df.filter(
        pl.col("bead_barcode").is_in(allowlist_bc)
    )
    quantification_df = quantification_df.sort("count", descending=True)
    return unfiltered_quantification_df, quantification_df


def create_count_df(quantification_df: pl.DataFrame) -> pl.DataFrame:
    """
    Function to create the count df from the quantification_df

    Args:
        quantification_df (pl.DataFrame): DataFrame containing the bead barcodes and their counts.
    Returns:
        pl.DataFrame: A DataFrame containing the bead barcodes and their counts multiplied by 2.
    """
    count_df = quantification_df.select(
        [pl.col("bead_barcode"), (pl.col("count") * 2)]
    )
    return count_df


def create_implicated_df(overlap_df: pl.DataFrame, count_df: pl.DataFrame) -> pl.DataFrame:
    """
    Function to create the implicated df.

    Args:
        overlap_df (pl.DataFrame): DataFrame containing the overlap information between barcodes.
        count_df (pl.DataFrame): DataFrame containing the bead barcodes and their counts.
    Returns:
        pl.DataFrame: A DataFrame containing pairs of barcodes with their Jaccard index and counts.
    """
    implicated_df = overlap_df.group_by(["barc1", "barc2"]).agg(
        pl.col("n_both").sum().alias("N_both")
    )
    implicated_df = implicated_df.join(
        count_df.collect(), left_on="barc1", right_on="bead_barcode"
    )
    implicated_df.columns = ["barc1", "barc2", "N_both", "N_barc1"]
    implicated_df = implicated_df.join(
        count_df.collect(), left_on="barc2", right_on="bead_barcode"
    )
    implicated_df.columns = ["barc1", "barc2", "N_both", "N_barc1", "N_barc2"]
    implicated_df = implicated_df.with_columns(
        (
            pl.col("N_both")
            / (pl.col("N_barc1") + pl.col("N_barc2") - pl.col("N_both") + 0.05)
        ).alias("jaccard_frag")
    )
    implicated_df = implicated_df.select(
        [
            pl.all().exclude("jaccard_frag"),
            pl.col("jaccard_frag").map_elements(lambda x: round(x, 5), return_dtype=pl.Float64),
        ]
    )
    implicated_df = implicated_df.filter(
        pl.col("jaccard_frag") > 0
    )  # Remove pairs with score of 0
    implicated_df = implicated_df.sort(
        "jaccard_frag", descending=True
    )  # Arrange from highest to lowest Jaccard index
    return implicated_df


def get_density_threshold(df: pl.DataFrame, call_type: str, logTransform: bool = False) -> Tuple[float, float]:
    """
    Function for calling jaccard knee based on local minima

    Args:
        df (pl.DataFrame): DataFrame containing the Jaccard index values.
        call_type (str): Type of call being made, either "jaccard" or "bead".
        logTransform (bool): Whether to apply log transformation to the Jaccard index values.
    Returns:
        Tuple[float, float]: A tuple containing the safety threshold and the knee threshold.
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
                .map_elements(lambda x: round(math.log10(x), 6), return_dtype=pl.Float64)
                .alias("log_counts"),
            ]
        )
    else:
        filtered_df = filtered_df.select(
            [pl.all(), pl.col("jaccard_frag").alias("log_counts")]
        )

    # Calculate the density using a gaussian kernel NEEDS TO BE STRIPPED OUT FROM R CODE
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


def get_mode(df: pl.DataFrame) -> float:
    """
    Function to determine mode of an input df column

    Args:
        df (pl.DataFrame): DataFrame containing the Jaccard index values.
    Returns:
        float: The mode of the Jaccard index values.
    """
    x = df["jaccard_frag"].to_list()
    ux = df["jaccard_frag"].unique().to_list()
    match = [ux.index(ind) if ind in ux else None for ind in x]
    tabulate = (np.bincount(match)).tolist()
    index = tabulate.index(max(tabulate))
    return ux[index]


def get_local_minima(x: list) -> list:
    """
    Get the local minima of a list of values

    Args:
        x (list): A list of numerical values.
    Returns:
        list: A list of indices where local minima occur.
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


def find_local_min_in_list(x: list, logTransform: bool, filtered_df_vec: vectors.FloatVector, vector_density_dict: dict) -> int:
    """
    Creating a new function to perform the find

    Args:
        x (list): List of indices where local minima occur.
        logTransform (bool): Whether the Jaccard index values are log-transformed.
        filtered_df_vec (vectors.FloatVector): Vector of Jaccard index values.
        vector_density_dict (dict): Dictionary containing the density values.
    Returns:
        int: The index of the first local minimum that meets the criteria.
    """
    # Make sure that the selected min includes at least 20% of barcodes
    # and that the difference between the min and the max differences
    # by some appreciable amount
    # both in terms of absolute difference (Caleb changed 0.5
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
    quantification_df: pl.DataFrame,
    barcode_filtered_df: pl.DataFrame,
    barcode_translate_df: pl.DataFrame,
    one_to_one: bool,
    name: str,
    ti_len: int,
    catac_assay: bool = False
):
    """
    Function to group the beads together and create the droplets

    Args:
        quantification_df (pl.DataFrame): DataFrame containing the bead barcodes and their counts.
        barcode_filtered_df (pl.DataFrame): DataFrame containing the barcodes that will be merged.
        barcode_translate_df (pl.DataFrame): DataFrame to store the mapping of bead barcodes to droplet barcodes.
        one_to_one (bool): Whether to keep a one-to-one mapping between beads and droplets.
        name (str): Name prefix for the output files.
        ti_len (int): Length of the Tn5 sequence if applicable.
        catac_assay (bool): Whether this is a cATAC assay, affecting droplet naming conventions.
    Returns:
        pl.DataFrame: DataFrame containing the mapping of bead barcodes to droplet barcodes.
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
                + barcode_merge[0][:int(ti_len)]
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
        barcode_translate_df = barcode_translate_df.with_columns(
            pl.when(
                pl.col("droplet_barcode").is_in(list(non_empty_old_droplets.keys()))
            )
            .then(pl.lit(droplet_name + "%02d" % droplet_size))
            .when(
                (pl.col("bead_barcode").is_in(barcode_merge))
                & (pl.col("droplet_barcode") == "")
            )
            .then(pl.lit(droplet_name + "%02d" % droplet_size))
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


def create_updated_deconvolution_params(
    jaccard_called_frag: float,
    min_jaccard_frag: float,
    name: str,
    input_params_suffix: str
) -> Tuple[float, pl.DataFrame]:
    """
    Function to write the updated deconvolution params file.

    Args:
        jaccard_called_frag (float): The Jaccard index value called from the knee.
        min_jaccard_frag (float): The minimum Jaccard index value.
        name (str): Name prefix for the output files.
        input_params_suffix (str): Suffix of the input deconvolution params file.
    Returns:
        Tuple[float, pl.DataFrame]: A tuple containing the bead threshold and the updated parameters DataFrame.
    """
    input_params_df = pl.read_csv(
        name + input_params_suffix,
        has_header=False,
        new_columns=["param", "value"],
    )
    bead_threshold = input_params_df.filter(
        pl.col("param") == "bead_threshold"
    )["value"].to_list()[0]
    bead_threshold = float(bead_threshold)
    # Add new rows to the params df
    output_params_df = input_params_df.cast({"value": pl.Utf8}).extend(
        pl.DataFrame({
            "param": ["jaccard_threshold_nosafety", "jaccard_threshold"],
            "value": [str(jaccard_called_frag), str(min_jaccard_frag)],
        }, strict=False)
    )
    return bead_threshold, output_params_df


def compile_deconvolution_metrics(
        sample: str,
        unfiltered_quant_df: pl.DataFrame,
        quantification_df: pl.DataFrame,
        bead_threshold: float,
        implicated_df: pl.DataFrame,
        jaccard_threshold: float,
        barcode_translate_df: pl.DataFrame
) -> pl.DataFrame:
    """
    Function to compile deconvolution metrics.

    Args:
        sample (str): Sample name.
        unfiltered_quant_df (pl.DataFrame): DataFrame with all bead barcodes and their counts.
        quantification_df (pl.DataFrame): Filtered DataFrame with bead barcodes that are in the allowlist.
        bead_threshold (float): The threshold for bead counts.
        implicated_df (pl.DataFrame): DataFrame containing pairs of barcodes with their Jaccard index and counts.
        jaccard_threshold (float): The Jaccard index threshold used for merging.
        barcode_translate_df (pl.DataFrame): DataFrame containing the mapping of bead barcodes to droplet barcodes.
    Returns:
        pl.DataFrame: A DataFrame containing the compiled metrics for deconvolution.
    """
    metrics_dict = {
        "total_barcodes_observed": unfiltered_quant_df.shape[0],
        "beads_above_knee": quantification_df.shape[0],
        "fragment_thresh_at_bead_knee": bead_threshold,
        "median_frags_per_above_knee_bead": quantification_df["count"].median(),
        "bead_pairs_above_jaccard_knee": implicated_df.filter(pl.col("merged") == "TRUE").shape[0],
        "jaccard_thresh_at_knee": jaccard_threshold,
        "total_cells": barcode_translate_df.select(pl.col("droplet_barcode")).unique().shape[0],
    }
    return pl.DataFrame({
        "sample": [sample] * len(metrics_dict),
        "process": ["deconvolution"] * len(metrics_dict),
        "metric": list(metrics_dict.keys()),
        "value": list(metrics_dict.values()),
    }, strict=False)


if __name__ == "__main__":
    main()
