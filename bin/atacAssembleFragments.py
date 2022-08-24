#!/usr/bin/env python3

"""
Filter reads in BAM file and convert to fragments file
"""

import csv
import gzip
from math import inf
from pysam import AlignmentFile
import click

@click.command("assemble_fragments")
@click.option(
    "--input",
    "-i",
    help="Path to an indexed BAM file with duplicates marked."
)
@click.option(
    "--sample",
    "-s",
    help="Sample ID string to use for labelling outputs.",
)
@click.option(
    "--chromosome",
    "-c",
    help="Reference contig to use for labelling outputs.",
)
@click.option(
    "--mapping-quality",
    "-q",
    default=30,
    help="Minimum mapping quality needed to retain a read.",
)
@click.option(
    "--insert-size",
    "-z",
    default=2000,
    help="Maximum insert size allowed for fragments.",
)
@click.option(
    "--chromosome-length",
    "-l",
    default=inf,
    help="Chromosome length.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Chromosome length.",
)
def assemble_fragments(
  input: str,
  sample: str,
  chromosome: str,
  mapping_quality: int,
  insert_size: int,
  chromosome_length: int,
  verbose: bool):
    """ Filter reads in BAM file and convert to fragments file """

    # Define input file
    bamfile = AlignmentFile(input, "rb")
    if verbose:
        print(f"verbose: input file: {input}")

    # Define fragment file path
    fragment_path = f"{input.replace('.bam', '')}.frag.bedpe.gz"
    if verbose:
        print(f"verbose: fragment path: {fragment_path}")
    input_count, output_count = write_to_output_file(
      bamfile,
      fragment_path,
      mapping_quality,
      insert_size,
      chromosome_length)
    write_to_readcount_file(
      sample,
      chromosome,
      input_count,
      output_count,
      verbose)

def write_to_output_file(
  bamfile,
  fragment_path: str,
  mapping_quality: int,
  insert_size: int,
  chromosome_length: int):
    """ filter reads in BAM file and convert to fragments file """
    # Initialize read counters
    input_count, output_count = 0, 0
    # Write to fragment file
    with gzip.open(fragment_path, 'wb') as fragment_file:
        for read_name, read_dict in bedpe_generator(bamfile, input_count):
    
            input_count += 1

            # If read 1 start is 5' of read 2 start OR read 2 is 5' of read 1
            # Filter reads based on quality score and insert size
            orientation = get_orientation(read_dict, mapping_quality, insert_size)
            if orientation == "skip":
                continue
            ref_start, ref_end = get_reference_values_by_orientation(
              orientation, read_dict, chromosome_length)

            output_count += 1

            # Format the fragment entry as tab-separated string
            entries = [
              read_dict['r1_ref_name'],
              str(ref_start),
              str(ref_end),
              read_name,
              orientation,
              str(read_dict['support'])]
            entry = '\t'.join(entries)

            # Write the fragment entry to file
            fragment_file.write(f"{entry}\n".encode('utf-8'))

    return input_count, output_count


def is_read1_start_5prime_of_read2_start(read_dict, mapping_quality, insert_size):
    """ If read 1 start is 5' of read 2 start """
    return (read_dict['r1_ref_start'] <= read_dict['r2_ref_start'] and
      min(read_dict['r1_map_qual'], read_dict['r2_map_qual']) >= mapping_quality and
      (read_dict['r2_ref_end'] - read_dict['r1_ref_start']) <= insert_size)

def is_read2_start_5prime_of_read1_start(read_dict, mapping_quality, insert_size):
    """ If read 2 start is 5' of read 1 start """
    return (read_dict['r1_ref_start'] > read_dict['r2_ref_start'] and
        min(read_dict['r1_map_qual'], read_dict['r2_map_qual']) >= mapping_quality and
        (read_dict['r1_ref_end'] - read_dict['r2_ref_start']) <= insert_size)

def get_orientation(read_dict, mapping_quality, insert_size):
    """ get the string orientation """
    orientation = "skip"
    if is_read1_start_5prime_of_read2_start(read_dict, mapping_quality, insert_size):
        orientation = "r1r2"
    elif is_read2_start_5prime_of_read1_start(read_dict, mapping_quality, insert_size):
        orientation = "r2r1"
    return orientation

def get_read_count_list(
  sample,
  chromosome,
  input_count,
  output_count):
    """ get read count list """
    # Making read count list
    read_count_list = [
        ["sample","process","metric","count"],
        [sample, "assemble_fragments_" + chromosome, "input", input_count],
        [sample, "assemble_fragments_" + chromosome, "output", output_count]
    ]
    return read_count_list

def write_to_readcount_file(
  sample: str,
  chromosome: str,
  input_count: int,
  output_count: int,
  verbose: bool):
    """ write input/output counts to readcount file """
    # Define read count file path
    read_count_path = f"{sample}_{chromosome}_assemble_fragments_output_read_counts.csv"
    if verbose:
        print(f"verbose: write to {read_count_path}")
    read_count_list = get_read_count_list(
      sample,
      chromosome,
      input_count,
      output_count)

    # Write to read count file
    with open(read_count_path, "w", encoding="utf-8") as read_count_file:
        csv_out = csv.writer(read_count_file)
        csv_out.writerows(read_count_list)

def get_reference_values_by_orientation(orientation: str, read_dict, chromosome_length: int):
    """ Set fragment values based upon orientation """
    if orientation == "r1r2":
        ref_start = read_dict['r1_ref_start']
        ref_end = read_dict['r2_ref_end']
        # Modify reference position based on strand
        if read_dict['r1_is_fwd']:
            ref_start += 4
            ref_end += 4
        else:
            ref_start -= 5
            ref_end -= 5
    elif orientation == "r2r1":
        ref_start = read_dict['r2_ref_start']
        ref_end = read_dict['r1_ref_end']
        # Modify reference position based on strand
        if read_dict['r2_is_fwd']:
            ref_start += 4
            ref_end += 4
        else:
            ref_start -= 5
            ref_end -= 5
    ref_start = max(ref_start, 0)
    ref_end = min(ref_end, chromosome_length)
    return ref_start, ref_end

def update_paired_read_dict(read_name: str, read_dict, read):
    """ update paired read dict """
    if read.is_read1:
        read_dict[read_name]['r1_ref_name'] = read.reference_name
        read_dict[read_name]['r1_ref_start'] = read.reference_start
        read_dict[read_name]['r1_ref_end'] = read.reference_end
        read_dict[read_name]['r1_map_qual'] = read.mapping_quality
        read_dict[read_name]['r1_is_fwd'] = not read.is_reverse
        read_dict[read_name]['support'] = read.get_tag('DS') if read.has_tag('DS') else 1
    else:
        read_dict[read_name]['r2_ref_name'] = read.reference_name
        read_dict[read_name]['r2_ref_start'] = read.reference_start
        read_dict[read_name]['r2_ref_end'] = read.reference_end
        read_dict[read_name]['r2_map_qual'] = read.mapping_quality
        read_dict[read_name]['r2_is_fwd'] = not read.is_reverse

    return read_dict


def bedpe_generator(
    bamfile,
    input_count_int
):
    """ iterator for BAM records """

    # Initialize dictionary
    bedpe_dict = {}

    # Iterate through BAM records
    for read in bamfile:
        # Only keep properly paired primary alignments
        if (read.is_proper_pair
            and not read.is_unmapped
            and not read.is_secondary
            and not read.is_duplicate):
            # Increment read counter
            input_count_int += 1
            # Pull read name
            read_name = read.query_name
            # If this is the first read in the pair
            if read_name not in bedpe_dict:
                # Create a dictionary entry
                bedpe_dict[read_name] = {}
                # Add read data to dictionary
                bedpe_dict = update_paired_read_dict(read_name, bedpe_dict, read)
                # If the read's mate is already in the dictionary
            else:
                # Complete the dictionary record with the mate data
                bedpe_dict = update_paired_read_dict(read_name, bedpe_dict, read)
                # Output the read name and completed dictionary record
                yield read_name, bedpe_dict[read_name]
                # Delete the dictionary record after returning it
                del bedpe_dict[read_name]


def main():
    """ entry point for assemble_fragments """
    assemble_fragments()

if __name__ == "__main__":
    main()
