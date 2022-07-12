#!/usr/bin/env python3
"""
Gets fragment counts from a BAM file in BED regions
"""
import bx.intervals as bxi
import click
import pysam
import os
from collections import defaultdict, OrderedDict
from scipy.sparse import lil_matrix, spmatrix, hstack
from scipy.io import mmwrite
from dataclasses import dataclass
from typing import DefaultDict, Dict, Tuple
import numpy as np
import pandas as pd
import ray


@dataclass
class CellMetaData:
    """
    Hold metadata associated with a given cellid.
    """

    __slots__ = ["index", "total_reads", "reads_in_peaks"]
    index: int
    total_reads: int
    reads_in_peaks: int

    def __iadd__(self, new):
        if self.index != new.index:
            raise ValueError("Mismatched indexes")
        new_total = self.total_reads + new.total_reads
        new_rip = self.reads_in_peaks + new.reads_in_peaks
        return CellMetaData(self.index, new_total, new_rip)


def build_intervals(
    bed_file: str, contig: str
) -> Tuple[DefaultDict[str, bxi.Intersecter], Dict[str, int]]:
    """
    Create the interval trees, one per chromosome.
    
    Note that the return type here is just Dict due to limitations in Pythons typing
       
    intersecters is a dictionary that has one key per chromosome, and an interval tree underneath it.
    interval_index_map is maintains an ordering of interval to index as originally found in the bed file
    """
    intersecters: DefaultDict[str, bxi.Intersecter] = defaultdict(bxi.Intersecter)
    interval_index_map: OrderedDict[str, int] = OrderedDict()

    with open(bed_file, "r") as bed_fh:
        count = 0
        for line in bed_fh:
            fields = line.strip().split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom == contig:
                name = "{}:{}-{}".format(chrom, start, end)
                interval = bxi.Interval(start, end, name)
                intersecters[chrom].add_interval(interval)
                interval_index_map[name] = count
                count += 1
    return intersecters, interval_index_map


def get_chromosomes_with_alignments(bam: str) -> list:
    """
    Simple function for getting the chromosomes from a BAM that have alignments
    """
    idx = pysam.idxstats(bam).splitlines()
    chroms2use = []
    for chr in range(len(idx)):
        idxstat = idx[chr].split("\t")
        if int(idxstat[2]) + int(idxstat[3]) > 0:
            chroms2use.append(idxstat[0])
    return chroms2use


def get_chromosomes_with_peaks(bed_file) -> list:
    """
    Simple function for getting the chromosomes with a peak called from a BED file
    """
    chromsWithPeaks = []
    with open(bed_file, "r") as bed_fh:
        for line in bed_fh:
            fields = line.strip().split("\t")
            chrom = fields[0]
            if chrom not in chromsWithPeaks:
                chromsWithPeaks.append(chrom)
    return chromsWithPeaks


def list_extract(lst: list, index: int) -> list:
    return [item[index] for item in lst]


@ray.remote
def detect_overlaps(
    bam_path: str,
    contig: str,
    chromvar_compat: bool,
    goodbc: list,
    peak_bed: str,
    verbose: bool,
) -> Tuple[spmatrix, Dict[str, CellMetaData], list]:
    """
    Iterate over SAM record and check for overlaps in the read fragments against the intervals from the bed file.
    Note that interval_index_map and the returned dict are OrderedDicts
    The read fragment refers to the start of the read1 to the end of the read2.
    """
    # set up interval tree
    intersecters, interval_index_map = build_intervals(peak_bed, contig)

    # set up index map
    cellid_index_map: OrderedDict[str, CellMetaData] = OrderedDict()

    # add entry to index map for each barcode on the allowlist
    for barcode in goodbc:
        cellid_index_map[barcode] = CellMetaData(barcode, 0, 0)

    # change to a set
    barcodes = set(goodbc)

    # create list for looking up cell IDs and zip to dictionary for faster lookups
    cellid_index_list = list(cellid_index_map.keys())
    cellid_indexes = list(range(0, len(cellid_index_list)))
    barcode_lookup = dict(zip(cellid_index_list, cellid_indexes))

    # create empty sparse matrix of dimensions barcode-by-peak
    matrix: spmatrix = lil_matrix(
        (
            len(goodbc),
            len(list(filter(lambda x: contig in x, list(interval_index_map.keys())))),
        ),
        dtype=np.uint,
    )

    # open alignment file for reading
    samfile = pysam.AlignmentFile(bam_path, "rb")

    # fetch the contig in question
    Itr = samfile.fetch(str(contig))

    # iterate over reads
    if verbose:
        print("Analyzing {}".format(contig))
        count = 0

    for read in Itr:
        if verbose:
            count += 1
            if count % 1000000 == 0:
                print("Processed {} reads on {}".format(count, contig))
        if read.has_tag("DB") and read.is_proper_pair and not read.is_reverse:
            cell_tag = read.get_tag("DB")
            if cell_tag in barcodes:
                fragment_start = read.reference_start
                fragment_end = fragment_start + abs(
                    read.template_length
                )  # - 1 # since pysam is converting to 0-based start pos, we don't need this... I think, output matches chromvar  # The -1 needs some research, just copying chromVar? ^ research
                if contig in intersecters:  # The chromosome is in peaks
                    # start = time.process_time()
                    intersections = intersecters[contig].find(
                        fragment_start, fragment_end
                    )
                    # print("Time for intersection {}".format(time.process_time() - start))
                    # Filter to produce exact output as chromvar
                    # Don't count fragments that totally encompass a region, only count if htey have a start or end within the region
                    # This should be removed at a future point
                    if chromvar_compat:
                        intersections = [
                            interval
                            for interval in intersections
                            if fragment_start >= interval.start
                            or fragment_end <= interval.end
                        ]
                    if intersections:  # Our fragment intersects peaks
                        for interval in intersections:
                            # Bump count in each peak that was intersected
                            # start1 = time.process_time()
                            matrix[
                                barcode_lookup.get(cell_tag),
                                interval_index_map[interval.value],
                            ] += 1
                            # print("Time for storing results {}".format(time.process_time() - start1))
                            # count the read as in a peak
                            if (
                                chromvar_compat
                            ):  # Chromvar appears to count add each time it sees a read in a peak, even that means double couting it
                                cellid_index_map[cell_tag].reads_in_peaks += 1
                        if (
                            not chromvar_compat
                        ):  # This is where a read in peak should be counted
                            cellid_index_map[cell_tag].reads_in_peaks += 1
                        cellid_index_map[cell_tag].total_reads += 1
                    else:  # read not in peak
                        cellid_index_map[cell_tag].total_reads += 1
                else:  # read not in peak
                    cellid_index_map[cell_tag].total_reads += 1
    if verbose:
        print("Completed {}".format(contig))
    return matrix.tocsr(), cellid_index_map, list(interval_index_map.keys()), contig


@click.command()
@click.option(
    "--bam_file", "-b", type=click.Path(), required=True, help="Path to bam file"
)
@click.option(
    "--allowlist",
    "-w",
    type=click.Path(),
    required=True,
    help="allowlist of cell barcodes, one per line",
)
@click.option(
    "--peak_bed", "-p", type=click.Path(), required=True, help="The peaks bed file"
)
@click.option(
    "--index",
    "-i",
    default="MULTIPLEXED",
    help="Combinatorial index within which to create matrix.",
)
@click.option(
    "--out_dir",
    "-o",
    type=click.Path(),
    required=True,
    help="The directory to write output to",
)
@click.option(
    "--chromvar_compat",
    is_flag=True,
    required=False,
    help="Force output to be what chromvar would output - arguably wrong",
)
@click.option(
    "--cpus", "-c", default=1, help="Number of CPUs to use in parallel processing."
)
@click.option(
    "--sampleid",
    "-s",
    required=True,
    help="sampleId for the current sample to be prepended to read_counts.txt",
)
@click.option(
    "--verbose",
    is_flag=True,
    required=False,
    help="Prints an update ever 1M reads processed per contig.",
)
def main(bam_file, allowlist, peak_bed, out_dir, index, chromvar_compat, cpus, sampleid, verbose):
    """
    Build a counts matrix of cell tags (from the DB tag) in the sam file against the peaks in the bed file.

    The count matrix has rows that are the cellids from the DB tags, and colums that are each range from the bed file.
    """
    # Read allowlist; All barcodes in this file are good barcodes
    wl = pd.read_csv(allowlist, delimiter=",")

    # Filter barcodes to include only tagmentation index if necessary
    goodbc = wl.DropBarcode

    if not index == "MULTIPLEXED":
        goodbc = goodbc[goodbc.str.contains(index)]

    goodbc = goodbc.tolist()

    # Get intersection of contigs with peaks and contigs with alignments
    contigs = get_chromosomes_with_alignments(bam_file)
    peaks = get_chromosomes_with_peaks(peak_bed)
    contigs_to_use = list(set(contigs) & set(peaks))

    if verbose:
        print("Counting reads in peaks on {}".format(list(contigs_to_use)))

    # Set up parallel processing
    ray.init(num_cpus=cpus)

    # Process contigs in parallel
    results = ray.get(
        [
            detect_overlaps.remote(
                bam_file, contig, chromvar_compat, goodbc, peak_bed, verbose
            )
            for contig in contigs_to_use
        ]
    )

    # Horizontally stack (column bind) outputs; convert back to csr
    out_mtx = hstack(list_extract(results, 0)).tocsr()

    # Sort output; hstack does not reorder rows
    out_mtx.sort_indices()

    # Write matrix to disk
    mmwrite(os.path.join(out_dir, sampleid + ".count_matrix"), out_mtx, field="integer")

    # Write matrix as a dense matrix in a CSV - Kept for debug purposes
    # with open(os.path.join(out_dir, "dense_matrix.tsv"), "w") as dense_fh:
    #     for cell_id, metadata in cell_index_map.items():
    #         row = []
    #         for interval, j in interval_index_map.items():
    #             row.append("{}".format(str(matrix[metadata.index, j])))
    #         row = '\n'.join(row)
    #         dense_fh.write(f"{row}\n")

    # Get peak names --> sort these by contigs_to_use
    all_peaks = list_extract(results, 2)
    peak_names = [item for sublist in all_peaks for item in sublist]

    # Write column names to a file; columns are the peak names
    with open(os.path.join(out_dir, sampleid + ".column_names.txt"), "w") as column_fh:
        for peak in peak_names:
            column_fh.write("{}\n".format(peak))

    # Write rows to a file; rows are the allowlist barcodes
    with open(os.path.join(out_dir, sampleid + ".row_names.txt"), "w") as row_fh:
        for line in goodbc:
            row_fh.write("{}\n".format(line))

    # Merge all the cell metadata
    meta = list_extract(results, 1)
    merged_dict = meta[0]
    for i in range(1, len(meta)):
        for barcode in goodbc:
            merged_dict[barcode] += meta[i][barcode]

    # Write Cell info to a file
    with open(os.path.join(out_dir, sampleid + ".read_counts.txt"), "w") as read_counts_fh:
        read_counts_fh.write("cellid\ttotal_reads\treads_in_peaks\n")
        for cellid, meta_data in merged_dict.items():
            read_counts_fh.write(
                "{}\t{}\t{}\n".format(
                    cellid, meta_data.total_reads, meta_data.reads_in_peaks
                )
            )


if __name__ == "__main__":
    main()
