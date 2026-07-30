#!/usr/bin/env python3

import pysam
import click
import ray
import glob
import sys
from ray.exceptions import ObjectStoreFullError
from coreBamCounter import getChromosomesWithAlignments


@click.command()
@click.option("--bamfile", "-b", help="A .bam file with an index.")
@click.option(
    "--cell-barcode-tag", "-ct", default="XC", help="Tag to use for the cell barcode."
)
@click.option(
    "--cell-barcode-position",
    "-cp",
    default=0,
    help="In an underscore delimited read name, the zero-based field that "
    "contains the cell barcode.",
)
@click.option("--umi-tag", "-ut", default="XZ", help="Tag to use for the UMI.")
@click.option(
    "--umi-position",
    "-up",
    default=1,
    help="In an underscore delimited read name, the zero-based "
    "field that contains the UMI.",
)
@click.option(
    "--tmp-dir", "-t", default="./tmp", help="Temp directory for split contig bams."
)
@click.option(
    "--out-dir", "-o", default="./", help="Output directory for final bam file."
)
@click.option(
    "--cpus", "-c", default=1, help="The number of CPUs for parallel processing."
)
@click.option(
    "--lookup", "-l", default=None, help="Barcode lookup table for bead merging."
)
def main(
    bamfile,
    cell_barcode_tag,
    cell_barcode_position,
    umi_tag,
    umi_position,
    tmp_dir,
    out_dir,
    cpus,
    lookup,
):
    # build parallel pool
    contigs = getChromosomesWithAlignments(bamfile)
    ray.init(num_cpus=cpus)

    if lookup is not None:
        # create lookup and update cell count for later use
        lookup_dict = build_barcode_lookup(lookup)

        try:
            ray.get(
                [
                    annotate_cell_barcodes.remote(
                        bamfile, contig, tmp_dir, cell_barcode_tag, lookup_dict
                    )
                    for contig in contigs
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

    else:
        try:
            ray.get(
                [
                    process_one_contig.remote(
                        bamfile,
                        contig,
                        tmp_dir,
                        cell_barcode_tag,
                        cell_barcode_position,
                        umi_tag,
                        umi_position,
                    )
                    for contig in contigs
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

    cat_args = glob.glob(tmp_dir + "/*.bam")
    with open(tmp_dir + "/samtools_cat.txt", 'w') as p:
        p.write(" ".join(cat_args))
    ray.shutdown()


def get_contigs_with_reads(bamfile):
    """Get the contigs from an indexed bam file that have alignments"""

    idxstats = pysam.idxstats(bamfile).strip().split("\n")
    chroms = []

    for line in idxstats:
        contig, len, nmapped, nunmapped = line.split("\t")
        if int(nmapped) > 0:
            chroms.append(contig)

    return chroms


def build_barcode_lookup(lookup):
    """
    Read barcode lookup table into dictionary
    Barcode lookup is formatted upstream in rnaMergeBeads.R
    Assumes columns are bead barcode, cell barcode, n beads in partition
    """
    # Initalize output varaibles
    output_count_dict = {}
    # Step through lookup table
    with open(lookup) as f:
        for line in f:
            (bead, cell, beads) = line.split(",")
            output_count_dict[bead] = cell
    return output_count_dict


@ray.remote
def process_one_contig(
    bamfile,
    contig,
    tmp_dir,
    cell_barcode_tag,
    cell_barcode_position,
    umi_tag,
    umi_position,
):
    """
    Iterates over a contig and moves barcodes and UMI from read name to desired tags
    Assumes the non-barcode/umi portion of the read name follows immediately
    """
    # fetch the contig

    bam = pysam.AlignmentFile(bamfile, "rb")

    Itr = bam.fetch(str(contig), multiple_iterators=True)

    # create a file to write outputs to
    splitBam = pysam.AlignmentFile(f"{tmp_dir}/{str(contig)}.bam", "wbu", template=bam)
    # open file for writing bead barcodes
    with open("all_beads.tsv", "a") as cell_barcode_output:
        # start processing reads
        for read in Itr:
            # split the read name by underscore
            rname = read.query_name.split("_")

            # slice relevant fields
            barcode = rname[cell_barcode_position]
            umi = rname[umi_position]
            readname = rname[max(umi_position, cell_barcode_position) + 1:]

            # set tags as desired
            read.set_tag(cell_barcode_tag, barcode)
            read.set_tag(umi_tag, umi)

            # set readname
            read.query_name = readname[0]

            # write alignment
            splitBam.write(read)
            cell_barcode_output.write(f"{barcode}\n")

    splitBam.close()


@ray.remote
def annotate_cell_barcodes(bamfile, contig, tmp_dir, cell_barcode_tag, lookup_dict):
    """
    Iterates over a contig
    Moves existing cell barcode to new tag
    Adds merged cell barcode to original cell barcode tag
    If bead barcode is not in lookup_dict (i.e. there was no deconvolution
    oligo associated with that bead), the bead barcode is used
    as the cell barcode.
    """
    # fetch the contig

    bam = pysam.AlignmentFile(bamfile, "rb")

    Itr = bam.fetch(str(contig), multiple_iterators=True)

    # Create a file to write outputs to
    splitBam = pysam.AlignmentFile(
        "{}/{}.bam".format(tmp_dir, str(contig)), "wbu", template=bam
    )

    # start processing reads

    for read in Itr:
        # check for bead barcode, move to new tag if found
        if read.has_tag(cell_barcode_tag):
            # get cbc
            old_cbc = read.get_tag(cell_barcode_tag)

            # look up cbc
            new_cbc = lookup_dict.get(old_cbc)

            if new_cbc is not None:
                read.set_tag(cell_barcode_tag, new_cbc)
                read.set_tag("XB", old_cbc)
            else:
                read.set_tag("XB", old_cbc)

        # write alignment
        splitBam.write(read)

    splitBam.close()


if __name__ == "__main__":
    main()
