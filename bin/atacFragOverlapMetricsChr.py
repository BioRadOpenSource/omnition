#!/usr/bin/env python3

import argparse
import csv
import gzip
import pandas as pd
from collections import Counter, defaultdict


# This function defines the input arguments.
def getargs():
    parser = argparse.ArgumentParser(
        description="Script for a given annotated chromosome fragments file saves a "
                    ".csv.gz file for each pair of barcodes given how many times "
                    "they share an identical fragment in this given chromosome."
    )

    options = parser.add_argument_group("Options")
    options.add_argument("-nc", "--ncthreshold", required=True, help="nc threshold")
    options.add_argument(
        "-hq",
        "--hqbeadsfile",
        required=True,
        help="hq beads file/barcode allow list file",
    )
    options.add_argument(
        "-be", "--bedpefile", required=True, help="annotated fragment file, bedpe"
    )
    options.add_argument(
        "-csv", "--csvfile", required=True, help="output .csv file name"
    )
    options.add_argument("-c", "--countfile", required=True, help="nc count file, .tsv")
    options.add_argument(
        "-rt", "--regularizethreshold", required=True, help="regularize threshold"
    )
    options.add_argument("-m", "--mergemethod", required=True, help="merge method")
    options.add_argument("-r", "--rounding", required=True, help="rounding value")

    return parser.parse_args()


def import_frags_export_NC():
    # Import fragments
    frags = gzip.open(args.bedpefile, "rb")
    frags_tsv = pd.read_csv(
        frags,
        delimiter="\t",
        header=None,
        names=[
            "read_name",
            "chr",
            "start",
            "end",
            "orientation",
            "support",
            "bead_barcode",
        ],
    )

    HQbeads = pd.read_csv(args.hqbeadsfile, header=None, names=["beads"])

    # Filter 1 for eligible barcodes
    frags_filt1 = frags_tsv[frags_tsv.bead_barcode.isin(HQbeads.beads)]

    # Quantify NC + export
    nccounts = (
        frags_filt1.groupby(["start", "end"], sort=False)
        .size()
        .reset_index(name="n_distinct_barcodes")
    )
    frags_filt1 = pd.merge(
        frags_filt1,
        nccounts,
        how="left",
        left_on=["start", "end"],
        right_on=["start", "end"],
    )
    del nccounts

    # Filter out high NC values + first pass PCR duplicate removal at the bead level
    nc_value = frags_filt1[frags_filt1["n_distinct_barcodes"] <= int(args.ncthreshold)]
    frags_filt2 = (
        nc_value.groupby(
            ["chr", "start", "end", "support", "bead_barcode", "orientation"],
            sort=False,
        )
        .size()
        .reset_index(name="PCRdupCount")
    )

    return frags_filt2


def make_insert_overlap(element_in_table, frags_filt2):
    """
    :description: Pull out barcode pairs/counts for retianed fragments
    :param element_in_table: Which column to use for the """\
    """insertion site (start, end, or cutsite)
    :param frags_filt2: dataframe containing bead barcodes,"""\
    """ insertion sites, and chromosome.
    :return: dictionary containing the implicated pair counts
    """

    # Find find overlaps of Tn5 insertions
    inserts_df = frags_filt2[
        ["chr", element_in_table, element_in_table, "bead_barcode"]
    ].copy()
    inserts_df.columns = ["chr", "start", "end", "bead_barcode"]

    # group by the insertion start site, and get all
    # combinations that have fragments overlapping
    groups = (
        inserts_df.groupby("start", sort=False)["bead_barcode"]
        .agg(lambda x: x.tolist())
        .reset_index(name="hits")
    )
    # merge the list of "hits"/barcodes for each row back into the inserts_df
    inserts_df = pd.merge(
        inserts_df, groups, how="left", left_on=["start"], right_on=["start"]
    )

    del groups  # remove groups dataframe that is no longer needed

    # create dictionary associating the current row barcode to the list
    # of barcodes with the same insertion start site
    inserts_df["barc"] = inserts_df.apply(
        lambda x: {x["bead_barcode"]: x["hits"]}, axis=1
    )

    # create a list of each barcode dictionary from each row
    bcs = list(inserts_df["barc"])

    implicatedPaircount = defaultdict(int)

    # iterate over the list of dictionaries
    for i in range(len(bcs)):
        for key, value in bcs[i].items():
            barcodes = list(
                filter(lambda x: x != key, value)
            )  # filter out reads that map to the same barcode
            # concatenate the barcodes to create a single key since
            # the barcodes are treated as an implicated pair
            counter = Counter(
                list(
                    map(
                        lambda x, y: x + "," + y,  # count the implicated pairs
                        list(
                            map(lambda x: key if key > x else x, barcodes)
                        ),  # update order of bc1
                        list(map(lambda x: x if key > x else key, barcodes)),
                    )
                )
            )  # update order of bc2
            for k, v in counter.items():
                implicatedPaircount[k] += v

    return dict(implicatedPaircount)


args = getargs()

frags_filt2 = import_frags_export_NC()

if len(frags_filt2.index) > 0:
    # handle rounding of inserts sites; rounding to -1,-2,-3 depending
    # on if its 10, 100, or 1000 respectively
    if int(args.rounding) != 0:
        rounding = (
            -1
            if int(args.rounding) == 10
            else (-2 if int(args.rounding) == 100 else -3)
        )
        frags_filt2 = frags_filt2.apply(
            lambda x: round(x, rounding) if x.name in ["start", "end"] else x
        )
    # create a column in frags_filt2 that has the coordinate of the relevant insertion
    if args.mergemethod == "r1":
        frags_filt2["cutsite"] = frags_filt2.apply(
            lambda x: x["end"] if x["orientation"] == "r2r1" else x["start"], axis=1
        )
        implicatedPairs = make_insert_overlap("cutsite", frags_filt2)
    elif args.mergemethod == "r2":
        frags_filt2["cutsite"] = frags_filt2.apply(
            lambda x: x["end"] if x["orientation"] == "r1r2" else x["start"], axis=1
        )
        implicatedPairs = make_insert_overlap("cutsite", frags_filt2)
    else:
        implicatedPaircount_start = make_insert_overlap("start", frags_filt2)
        implicatedPaircount_end = make_insert_overlap("end", frags_filt2)
        implicatedPairs_start_end = [implicatedPaircount_start, implicatedPaircount_end]

        # create a single dictionary, adding values from Counter dictionaries
        # with the same key (the implicated pair)
        implicatedPairs = defaultdict(int)
        for x in implicatedPairs_start_end:
            for k, v in x.items():
                implicatedPairs[k] += v

        implicatedPairs = dict(implicatedPairs)

    # Divide each count by (following the previous  data.table .(n_both = .N/2) syntax)
    for key, val in implicatedPairs.items():
        implicatedPairs[key] = int(val / 2)

    # filter out n_both values that are < regularize threshold
    implicatedPairs = dict(
        filter(lambda x: x[1] >= int(args.regularizethreshold), implicatedPairs.items())
    )

else:
    implicatedPairs = {}

# Export
try:
    with gzip.open((args.csvfile), "wt") as f:
        writer = csv.writer(f)
        writer.writerow(["barc1", "barc2", "n_both"])
        for k, v in implicatedPairs.items():
            barc = k.split(",")
            writer.writerow([barc[0], barc[1], v])
except IOError:
    print("I/O error")
