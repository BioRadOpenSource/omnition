#!/usr/bin/env python3

import click
from umi_tools import UMIClusterer
import math
import ray
from ray.exceptions import ObjectStoreFullError
import sys


@click.command()
@click.option(
    "--edge_file",
    "-e",
    help="Edge list file.",
)
@click.option(
    "--sample",
    "-s",
    help="Sample ID string to use for labelling outputs.",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="Print progress to the screen.",
)
@click.option(
    "--umi_threshold",
    "-u",
    default=1,
    help="Hamming distance threshold for correcting UMIs. If the distance between \
    two UMIs is less than or equal to this value, they will be corrected.",
)
@click.option("--cpus", "-c", default=1, help="The number of cpus to be used.")
def main(edge_file, sample, verbose, umi_threshold, cpus):
    if verbose:
        print("Reading edgelist...")

    # load the bead dicts and read in the edge list
    edges = make_edge_dict(edge_file)

    ray.init(num_cpus=cpus)
    updated_counts = []

    if umi_threshold == 0:
        c_umi_edges = edges
        updated_counts.append(0)
        if verbose:
            print("Skipping Correct UMIs...")
    else:
        if verbose:
            print("Correcting UMIs...")
        try:
            # Split the edge dictionary into equal sized chunks
            # based on the number of cpus supplied.
            # Run each chunk on its own thread through correct_umis
            corrected_umis = ray.get(
                [
                    correct_umis.remote(chunk, umi_threshold)
                    for chunk in split_dict(edges, cpus)
                ]
            )
        except ObjectStoreFullError:
            ray.shutdown()
            print("OOM error: ObjectStoreFullError")
            sys.exit(137)
        except Exception as err:
            ray.shutdown()
            print(f"Unexpected {err=}, {type(err)=}")
            sys.exit(1)
        # Combine the edge dictionaries and corrected read counts
        c_umi_edges = {}
        corrected_umi_count = 0
        for i in range(len(corrected_umis)):
            c_umi_edges.update(corrected_umis[i][0])
            corrected_umi_count += corrected_umis[i][1]

        if verbose:
            print("{} UMIs corrected".format(corrected_umi_count))
        updated_counts.append(corrected_umi_count)

    # Write outputs
    if verbose:
        print("Writing outputs...")
    write_corrected_edges(c_umi_edges, sample)
    write_counts_file(updated_counts, sample)


""" Begin Input data parsers """


# Edgelist is formatted:
# bead barcode1, umi1, bead barcode2, umi2
# Read these into a dictionary
def make_edge_dict(edge_file):
    bead_dict = {}
    with open(edge_file) as f:
        for line in f:
            b1, u1, b2, u2 = line.strip().split("\t")
            # Sort beads and umis lexographically
            # so we don't end up with duplicates
            bead1, bead2, umi1, umi2 = order_bead_umi_pair(b1, u1, b2, u2)
            # clusterer only takes bytes, not strings
            # setting up umi_pair with _ so we can
            # split it back into two sequences later
            umi_pair = bytes(umi1 + "_" + umi2, "utf-8")
            # Create dictionary for bead1 if not available
            if bead1 not in bead_dict.keys():
                bead_dict[bead1] = {}
            # Create dictionary for bead1 if not available
            if bead2 not in bead_dict[bead1].keys():
                bead_dict[bead1][bead2] = {}
            # Set up count for umi_pair if not already done
            if umi_pair not in bead_dict[bead1][bead2].keys():
                bead_dict[bead1][bead2][umi_pair] = 0
            # Increment UMI count
            bead_dict[bead1][bead2][umi_pair] += 1
    return bead_dict


""" End data parsers """


""" Begin helper functions """


# Sort bead and UMI into lexicographical order.
def order_bead_umi_pair(bead1, umi1, bead2, umi2):
    if bead1 == bead2:
        if umi1 < umi2:
            return bead1, bead2, umi1, umi2
        else:
            return bead2, bead1, umi2, umi1
    elif bead1 < bead2:
        return bead1, bead2, umi1, umi2
    else:
        return bead2, bead1, umi2, umi1


# Function for spliting dictionary into
# equal sized chunks
def split_dict(d, chunk_count):
    # Create a list of the keys
    keys = list(d.keys())
    # Divide list length by number of
    # chunks and round up
    n = math.ceil(len(keys) / chunk_count)
    # create dictionary chunks of size n
    for i in range(0, len(keys), n):
        yield {k: d[k] for k in keys[i: i + n]}


# Function for spliting dictionary into
# equal sized chunks
def split_list(large_list, chunk_count):
    # Divide list length by number of
    # chunks and round up
    n = math.ceil(len(large_list) / chunk_count)
    # create dictionary chunks of size n
    for i in range(0, len(large_list), n):
        yield large_list[i: i + n]


""" End helper functions """


""" Begin sequence correction functions """


# Function to correct BEAD UMI sequences based
# on a set Hamming distance threshold
# bead_dict is a dictionary of dictionaries
# first layer: key = bead, value = dictionary
# second layer: key = bead2, value = dictionary
# third layer: key = Bead UMI/Bead UMI pair, value = read count
@ray.remote
def correct_umis(bead_dict, umi_threshold):
    # Set up count of corrected reads
    corrected = 0
    # Establish clusterer class
    clusterer = UMIClusterer(cluster_method="directional")
    # Iterate over beads in the dictionary
    for bead in bead_dict.keys():
        # Iterate over either the second bead or the do
        # depending on the format
        for second in bead_dict[bead].keys():
            # If you have more than one UMI, start correction
            if len(bead_dict[bead][second].keys()) > 1:
                # Establish empty final UMI dictionary
                new_umis = {}
                # Cluster UMIs based on threshold and read count
                # Clusters are yielded as a list of lists
                # UMIs in their own list alone are not corrected
                # UMIs in a list with other UMIs are clustered
                # The first UMI in the list is the correct UMI
                clustered_umis = clusterer(
                    bead_dict[bead][second], threshold=umi_threshold
                )
                for cluster in clustered_umis:
                    # Establish a new UMI count
                    new_umis[cluster[0]] = 0
                    # Iterate over UMIs in the cluster and add their counts
                    # to the new total
                    for umi in cluster:
                        new_umis[cluster[0]] += bead_dict[bead][second][umi]
                        if umi != cluster[0]:
                            corrected += bead_dict[bead][second][umi]
                # Update the bead dictionary with the new UMIs
                bead_dict[bead][second] = new_umis
    return (bead_dict, corrected)


""" End sequence correction functions """


""" Begin output functions """


# Write a new edges file with the same number of reads
# as we started with, just now corrected sequences
def write_corrected_edges(edges_dict, sample_name):
    with open("{}_corrected_edges.tsv".format(sample_name), "w") as f:
        # Iterate over beads in the edges dict
        for bead in edges_dict.keys():
            # Iterate over bead2s
            for bead2 in edges_dict[bead].keys():
                # For each umi pair, split it back into two UMIs
                for umi_pair in edges_dict[bead][bead2].keys():
                    umi1, umi2 = umi_pair.decode("utf-8").split("_")
                    # Write out the read as many times as we have read counts
                    for _ in range(0, edges_dict[bead][bead2][umi_pair]):
                        f.write("{}\t{}\t{}\t{}\n".format(bead, umi1, bead2, umi2))


# Write a log file of the read counts corrected
def write_counts_file(updated_counts, sample_name):
    with open("{}_corrected_edge_counts.csv".format(sample_name), "w") as f:
        f.write("sample,process,read,metric,value\n")
        f.write(
            "{},correct_edges,NA,UMI_corrected_reads,{}\n".format(
                sample_name, updated_counts[-1]
            )
        )


""" End output functions"""

if __name__ == "__main__":
    main()
