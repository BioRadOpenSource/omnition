#!/usr/bin/env python3

"""
Filter reads in BAM file and convert to fragments file
"""

import click
import csv
import gzip
import pysam




@click.command()
@click.option(
    "--input", 
    "-i", 
    help="Path to an indexed BAM file with duplicates marked.")
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




def main(input, mapping_quality, insert_size, sample, chromosome):
    
    # Define input file
    bamfile = pysam.AlignmentFile(input, "rb")

    # Initialize read counters
    input_count = 0
    output_count = 0

    # Define fragment file path
    fragment_path = "{}.frag.bedpe.gz".format(input.replace(".bam", ""))

    # Write to fragment file
    with gzip.open(fragment_path, 'wb') as fragment_file:
        for read_name, read_dict in bedpe_generator(bamfile, input_count):

            # Increment input counter
            input_count += 1

            # If read 1 start is 5' of read 2 start
            # Filter reads based on quality score and insert size
            if (read_dict['r1_ref_start'] <= read_dict['r2_ref_start']
                and min(read_dict['r1_map_qual'], read_dict['r2_map_qual']) 
                    >= mapping_quality
                and (read_dict['r2_ref_end'] - read_dict['r1_ref_start']) 
                    <= insert_size):
                        # Set fragment values
                        ref_name = read_dict['r1_ref_name']
                        ref_start = read_dict['r1_ref_start']
                        ref_end = read_dict['r2_ref_end']
                        query_name = read_name
                        orientation = "r1r2"
                        support = read_dict['support']
                        # Modifying reference position based on strand
                        if read_dict['r1_is_fwd']:
                            ref_start += 4
                            ref_end += 4
                        else:
                            ref_start -= 5
                            ref_end -= 5

            # If read 2 start is 5' of read 1 start
            # Filter reads based on quality score and insert size
            elif (read_dict['r1_ref_start'] > read_dict['r2_ref_start']
                and min(read_dict['r1_map_qual'], read_dict['r2_map_qual']) 
                    >= mapping_quality
                and (read_dict['r1_ref_end'] - read_dict['r2_ref_start']) 
                    <= insert_size):

                        # Set fragment values
                        ref_name = read_dict['r1_ref_name']
                        ref_start = read_dict['r2_ref_start']
                        ref_end = read_dict['r1_ref_end']
                        query_name = read_name
                        orientation = "r2r1"
                        support = read_dict['support']

                        # Modifying reference position based on strand
                        if read_dict['r2_is_fwd']:
                            ref_start += 4
                            ref_end += 4
                        else:
                            ref_start -= 5
                            ref_end -= 5

            else:
                continue

            # Increment output counter
            output_count += 1

            # Format the fragment entry in as tab-separate
            entry = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                ref_name,
                ref_start,
                ref_end,
                query_name,
                orientation,
                support
            ).encode()

            # Write the fragment entry to file
            fragment_file.write(entry)

    # Define read count file path
    read_count_path = "{0}_{1}_assemble_fragments_output_read_counts.csv".format(
        sample,
        chromosome
    )

    # Making read count list
    read_count_list = [
        ["sample","process","metric","count"],
        [sample, "assemble_fragments_" + chromosome, "input", input_count],
        [sample, "assemble_fragments_" + chromosome, "output", output_count]
    ]

    # Write to read count file
    with open(read_count_path, "w") as read_count_file:
        csv_out = csv.writer(read_count_file)
        csv_out.writerows(read_count_list)



def update_paired_read_dict(
    read_name, 
    read_dict, 
    read
):

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

    return(read_dict)


def bedpe_generator(
    bamfile,
    input_count_int
):

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
                    bedpe_dict = update_paired_read_dict(read_name, 
                        bedpe_dict, 
                        read)
                # If the read's mate is already in the dictionary
                else:
                    # Complete the dictionary record with the mate data
                    bedpe_dict = update_paired_read_dict(read_name, 
                        bedpe_dict, 
                        read)
                    # Output the read name and completed dictionary record
                    yield read_name, bedpe_dict[read_name]
                    # Delete the dictionary record after returning it
                    del bedpe_dict[read_name]




if __name__ == "__main__":
    main()
