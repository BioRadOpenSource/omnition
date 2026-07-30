#!/usr/bin/env python3
# Helper script
import math
import pysam
from typing import List, Generator


def get_chromosomes_with_alignments(bam: str) -> List[str]:
    """
    Generate a list of chromosomes with alignments in a bam file.

    Args:
        bam (str): Path to the bam file.

    Returns:
        List[str]: List of chromosomes with alignments.
    """
    idx = pysam.idxstats(bam).splitlines()
    chroms_with_alignments = []
    for chr in range(len(idx)):
        idxstat = idx[chr].split("\t")
        if int(idxstat[2]) + int(idxstat[3]) > 0:
            chroms_with_alignments.append(idxstat[0])
    return chroms_with_alignments


def split_dict(d: dict, chunk_count: int) -> Generator[dict, None, None]:
    """
    Split a dictionary into equal sized chunks

    Args:
        d (dict): Dictionary to be split
        chunk_count (int): Number of chunks to create

    Returns:
        Generator[dict]: Generator of chunked dictionaries
    """
    # Create a list of the keys
    keys = list(d.keys())
    # Divide list length by number of
    # chunks and round up
    n = math.ceil(len(keys) / chunk_count)
    # create dictionary chunks of size n
    for i in range(0, len(keys), n):
        yield {k: d[k] for k in keys[i: i + n]}
