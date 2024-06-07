#!/usr/bin/env python3

import click
import pysam


# Get arguments
@click.command()
@click.option("--bam", "-b", help="Bam file to split.")
@click.option("--id", "-i", help="Sample ID.")
@click.option("--species1", "-s1", help="First species.")
@click.option("--species2", "-s2", help="Second Species")
def main(bam, id, species1, species2):

    # Open the bam file
    bamfile = read_bam(bam)

    # Create species 1 bam file to write to
    species_1_bam = make_bam(f"{species1}.mixed.filtered.bam", bamfile)

    # Create species 2 bam file to write to
    species_2_bam = make_bam(f"{species2}.mixed.filtered.bam", bamfile)

    # Sort reads by species
    sort_reads_by_species(bamfile, species1, species2, species_1_bam, species_2_bam)


def read_bam(bam_path: str) -> pysam.AlignmentFile:
    """
    Reads a bam file and returns a pysam.AlignmentFile object.

    Args:
    - bam_path (str): path to the bam file.

    Returns:
    - bam (pysam.AlignmentFile): a pysam.AlignmentFile object representing the bam file.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    return bam


def make_bam(bam_path: str, template_bam: pysam.AlignmentFile) -> pysam.AlignmentFile:
    """
    Creates a new bam file and returns a pysam.AlignmentFile object.

    Args:
    - bam_path (str): path to the bam file.
    - template_bam (pysam.AlignmentFile): a pysam.AlignmentFile object
        epresenting the bam file to use as a template.

    Returns:
    - bam (pysam.AlignmentFile): a pysam.AlignmentFile object representing the bam file.
    """
    bam = pysam.AlignmentFile(bam_path, "wb", template=template_bam)
    return bam


def sort_reads_by_species(
    input_bam: pysam.AlignmentFile,
    species1: str,
    species2: str,
    out_bam_1: pysam.AlignmentFile,
    out_bam_2: pysam.AlignmentFile,
):
    """
    Sorts reads in a bam file by species and writes them to two different bam files.

    Args:
    - input_bam (pysam.AlignmentFile): a pysam.AlignmentFile object
        representing the input bam file.
    - species1 (str): the name of the first species.
    - species2 (str): the name of the second species.
    - out_bam_1 (pysam.AlignmentFile): a pysam.AlignmentFile object
        representing the output bam file for species 1.
    - out_bam_2 (pysam.AlignmentFile): a pysam.AlignmentFile object
        representing the output bam file for species 2.

    Returns:
    - None
    """
    # Iterate through reads
    for read in input_bam.fetch(until_eof=True):
        # Get name of species
        name = read.reference_name

        # Write to different bam depeneding on what the species is
        if name is not None:
            if species1 in name:
                out_bam_1.write(read)
            elif species2 in name:
                out_bam_2.write(read)

    # Close the open Bam files
    input_bam.close()
    out_bam_1.close()
    out_bam_2.close()


if __name__ == "__main__":
    main()
