#!/usr/bin/env python3

import click
from dataclasses import dataclass
import csv


@click.command()
@click.option(
    "--edge_file",
    "-e",
    help="Edge list file",
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
    help="Number of UMIs required for a bead-bead or bead-DO connection "
    "to be considered real.",
)
@click.option(
    "--one_to_one",
    "-oo",
    is_flag=True,
    help="When this flag is present, each bead is mapped to a unique cell barcode. This"
    " effectively negates all bead merging.",
)
@click.option(
    "--bead_file",
    "-bf",
    help="Barcode file containing all bead barcodes.",
)
def main(edge_file, sample, verbose, umi_threshold, one_to_one, bead_file):
    if verbose:
        print("Reading edgelist...")

    # load the bead dicts and read in the edge list
    edges = read_edgelist(edge_file)
    edges = prune_edges(edges, umi_threshold)
    bead_dict = format_edges(edges)
    full_bead_list = read_bead_file(bead_file)

    # Create droplets
    if verbose:
        print("Building droplets...")
    if one_to_one:
        if verbose:
            print("One-to-one specified, mapping each bead to a unique droplet...")
        droplets = build_one_to_one_droplets(bead_dict)
    else:
        droplets = build_droplets(bead_dict, verbose)

    # Compare master bead list with droplet beads
    final_bead_list = compare_bead_files(full_bead_list, droplets)

    # Write outputs
    if verbose:
        print("Writing outputs...")
    write_final_bead_output(final_bead_list)
    write_barcode_translate(droplets, sample, one_to_one)


""" Begin Input data parsers """


def read_edgelist(edge_file):
    """
    Fill a dictionary with edge list contents.
    This creates a map of each R1 barcode to its R2 barcodes
    Filter UMIs that are ambiguous (=has an N)
    Filter self-connections
    """
    bead_dict = {}
    with open(edge_file) as f:
        for line in f:
            b1, u1, b2, u2 = line.strip().split("\t")
            if not b1 == b2:
                if "N" not in u1 and "N" not in u2:
                    # create a string for the UMI
                    bead1, bead2, umi = order_barcode_umi(b1, u1, b2, u2)
                    if bead1 not in bead_dict.keys():
                        bead_dict[bead1] = {}
                        bead_dict = add_new_connection(bead_dict, bead1, bead2, umi)
                    else:
                        # check if bead2 has already been seen
                        bead_connections = bead_dict[bead1].keys()
                        # if bead2 has already been seen
                        if bead2 in bead_connections:
                            # if bead1 has already seen bead2,
                            # only add a count if the umi is new
                            if umi not in bead_dict[bead1][bead2].umis:
                                bead_dict = add_new_umi_connection(
                                    bead_dict, bead1, bead2, umi
                                )
                        # if bead2 has not been seen, add a connection
                        else:
                            bead_dict = add_new_connection(bead_dict, bead1, bead2, umi)
                    if bead2 not in bead_dict.keys():
                        bead_dict[bead2] = {}
                        bead_dict = add_new_connection(bead_dict, bead2, bead1, umi)
                    else:
                        # check if bead1 has already been seen
                        bead_connections = bead_dict[bead2].keys()
                        # if bead1 has already been seen
                        if bead1 in bead_connections:
                            if umi not in bead_dict[bead2][bead1].umis:
                                bead_dict = add_new_umi_connection(
                                    bead_dict, bead2, bead1, umi
                                )
                        else:
                            bead_dict = add_new_connection(bead_dict, bead2, bead1, umi)
    return bead_dict


def read_bead_file(bead_file):
    """
    Read in the bead file list made upstream
    that contains a list of all beads. This will be
    used to filter to bead barcodes not found in the
    DO edge file
    """
    bead_set = set()
    with open(bead_file) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            bead_set.add(line[0])
    return bead_set


""" End data parsers """

""" Begin Dataclasses """


@dataclass
class Connection:
    """
    Dataclass to store bead-bead UMI count
    """

    __slots__ = ["umis", "count"]
    umis: set()
    count: int


@dataclass
class Droplet:
    """
    Dataclass to store beads for a partition.
    """

    __slots__ = ["beads"]
    beads: set()


""" Begin droplet builders """


def build_droplets(bead_dict, verbose):
    """
    Build droplets by connecting beads directly.
    How:
    1. Select a pair of beads from the bead:bead map.
    2. Look for other beads that "saw" the value in the bead:bead map.
    3. If found, add those beads and any new beads to the droplet.
    4. Repeat this process until no new beads are added to
    the droplet, meaning we have found all beads that
    the sequencing data say were together in the same droplet.
    5. Remove all beads and oligos in our droplet from the bead
    and oligo map dictionaries.
    6. Return these updated dictionaries to the global environment
    and use them for the next iteration.
    6. Repeat until no beads are present in the bead map, meaning
    all beads have been assigned to a droplet.
    """
    # create empty dict
    droplet_dict = {}
    # Start at 0 for the droplet IDs
    droplet_id = 0
    # Iterate while we have beads that have yet to be assigned
    # to a droplet
    while len(bead_dict.keys()) > 0:
        # Select the next bead barcode in the bead map
        barc = next(iter(bead_dict))
        # Assign a new droplet ID
        droplet_id += 1
        # Add the bead's barcode and oligos to a new Droplet
        bead_dict[barc].add(barc)
        droplet_dict.update({droplet_id: Droplet(bead_dict[barc])})
        # Provide an update, if desired
        if verbose:
            print(
                "Number of partitions built {}, Number of remaining beads: {}".format(
                    droplet_id, len(bead_dict.keys())
                ),
                end="\r",
                flush=True,
            )
        # Build a droplet using logic defined above
        droplet_dict, bead_dict_forward = expand_droplet(
            droplet_dict, bead_dict, droplet_id
        )
    return droplet_dict


def build_one_to_one_droplets(bead_dict):
    """
    Create a Droplet dict where each bead in the bead dict
    is placed into a unique droplet.
    """
    # create empty dict
    droplet_dict = {}
    droplet_id = 0
    for bead in bead_dict:
        droplet_id += 1
        droplet = Droplet(bead_dict[bead])
        droplet_dict.update({droplet_id: droplet})
    return droplet_dict


""" End droplet builders """

""" Begin helper functions """


def add_new_connection(bead_dict, k, v, umi):
    """
    Add a new Connection to a dictionary.
    k and v represent one bead and another bead
    """
    bead_dict[k].update({v: Connection(set(), 0)})
    bead_dict[k][v].umis.add(umi)
    bead_dict[k][v].count += 1
    return bead_dict


def add_new_umi_connection(bead_dict, k, v, umi):
    """
    Add a new UMI to an existing Connection.
    k and v can represent one bead and another bead
    """
    bead_dict[k][v].umis.add(umi)
    bead_dict[k][v].count += 1
    return bead_dict


def order_barcode_umi(bead1, umi1, bead2, umi2):
    """
    Sort bead and UMI into lexicographical order.
    """
    if bead1 < bead2:
        return bead1, bead2, f"{umi1}{umi2}"
    else:
        return bead2, bead1, f"{umi2}{umi1}"


""" End helper functions """


def prune_edges(bead_dict, umi_threshold):
    """
    Filters a dictionary produced by read_edgelist or
    read_do_edgelist to remove Connections with fewer than
    umi_threshold.
    """
    for bead in list(bead_dict.keys()):
        # iterate over bead-target
        for target in list(bead_dict[bead].keys()):
            # drop Connections with < umi_threshold umis
            if bead_dict[bead][target].count < umi_threshold:
                bead_dict[bead].pop(target)
    return bead_dict


def format_edges(out_dict):
    """
    Convert a dictionary of Connections to simple key:{values} dict
    """
    forward = {k: set(v.keys()) for k, v in out_dict.items()}
    return forward


def compare_bead_files(bead_list, droplet_dict):
    """
    Remove all instances of beads that map to edgelist from the master
    bead file.
    """
    for droplet in droplet_dict:
        for bead in droplet_dict[droplet].beads:
            if bead in bead_list:
                bead_list.remove(bead)
    return bead_list


""" Begin Droplet expansions functions """


def expand_droplet(droplet_dict, bead_dict, droplet_id):
    """
    Logic for this approach is described in build_droplets().
    """
    # Define the number of beads and deconvolution oligos in the partition
    # before potentially adding new items
    start_size = len(droplet_dict[droplet_id].beads)
    # Set counts for total beads to infinity in order
    # to get the while loop below to execute at least one time
    end_size = float("inf")
    # Begin iterating, only stop when no further additions to the
    # droplet happen
    while end_size > start_size:
        # Count the number of beads at the start of this iteration
        start_size = len(droplet_dict[droplet_id].beads)
        # Get the beads that are currently in the Droplet
        begin_beads = droplet_dict.get(droplet_id).beads.copy()
        # Iterate over each bead in the droplet to see if it brings any new beads
        for bead in begin_beads:
            # Get the beads seen by our bead that's already in the droplet
            seen_beads = bead_dict.get(bead)
            # If the bead saw no other beads, skip it
            if seen_beads is None:
                continue
            # If the bead saw new beads, add them to the droplet
            # Sets don't allow duplication, so no worries about repeats
            for seen_bead in seen_beads:
                droplet_dict[droplet_id].beads.add(seen_bead)
        # After iterating over all beads, see if the size has changed
        # If the droplet hasn't grown, it's a stable size and the loop breaks
        # If the droplet HAS grown, we'll repeat the loop
        end_size = len(droplet_dict[droplet_id].beads)
    # Once the droplet is stable, remove the beads in this droplet from the input dict
    for bead in droplet_dict[droplet_id].beads.copy():
        if bead in bead_dict.keys():
            del bead_dict[bead]
    return (droplet_dict, bead_dict)


""" End Droplet expansions functions """


""" Begin barcode formatting and output """


def format_droplet_barcode(sample_name, digits, droplet_id, total_beads):
    """
    Format the barcode name to match format that is the
    de facto standard
    """
    # Remove dashes and underscores
    fmt_sample_name = sample_name.replace("-", "").replace("_", "")

    # Pad droplet ID with leading zeros
    padded_id = format(droplet_id, "0{}".format(digits))

    # Slap it all together
    barcode = "{}BC{}N{}".format(fmt_sample_name, padded_id, total_beads)

    return barcode


def write_barcode_translate(droplet_dict, sample_name, one_to_one):
    """
    Write out a barcode translation file matching the
    de facto standard
    """
    n_digits = len(str(len(droplet_dict.keys())))
    with open("{}.barcodeTranslate.tsv".format(sample_name), "w") as f:
        f.write("{},{},{}\n".format("BeadBarcode", "DropBarcode", "nBeads"))
        for droplet in droplet_dict:
            if one_to_one:
                total_beads = 1
            else:
                total_beads = len(droplet_dict[droplet].beads)
            drop_barcode = format_droplet_barcode(
                sample_name, n_digits, droplet, total_beads
            )
            for bead in droplet_dict[droplet].beads:
                f.write("{},{},{}\n".format(bead, drop_barcode, total_beads))


def write_final_bead_output(beads):
    """
    Write the now filtered bead list to an output
    for use with generating cell_barcodes
    """
    with open("filtered_bead_list.csv", "w", newline="") as bead_file:
        for bead in beads:
            bead_file.write(f"{bead}\n")


""" End barcode formatting and output """

if __name__ == "__main__":
    main()
