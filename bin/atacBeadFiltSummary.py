#!/usr/bin/env python3

import argparse
import os
import statistics
import re

parser = argparse.ArgumentParser(
    description="Create a summary file of the bead filtration stats."
)

# Defining argument group 'Options'
options = parser.add_argument_group("Options")
options.add_argument(
    "-s", "--sample-id", required=True, help="The sample ID of the compiled sample."
)
options.add_argument(
    "-i",
    "--input-directory",
    required=True,
    help="Input directory containing deconvolution files.",
)


# Make a list of all fastq+TI pair names
def grabNames(inputdir):
    names = []
    dirfiles = os.listdir(inputdir)
    for i in range(len(dirfiles)):
        if dirfiles[i].find(".QCstats.csv") != -1:
            names.append(dirfiles[i].replace(".QCstats.csv", ""))
    names.sort()

    return names


# Get all sample files of a given type from a larger list of files
def getFileList(pattern, file_list):
    target_files = []
    for filename in file_list:
        reg_check = re.search(pattern, filename)
        if reg_check:
            target_files.append(reg_check[0])
    return target_files


# Get the threshold value for a given sample
# Take the median from the Unique fragment threshold at
# knee from all TIs associated with the sample
def getSampleThresh(inputdir):
    decon_params_pattern = r".*\.deconvolutionParams.csv$"
    dirfiles = os.listdir(inputdir)
    # Get just the deconvoultionParams files
    decon_params_files = getFileList(decon_params_pattern, dirfiles)

    # Put the bead thresholds from all the decon params files into a list
    all_thresholds = []
    for params_file in decon_params_files:
        with open(params_file, "r") as file:
            for line in file:
                row_data = line.strip().split(",")
                if row_data[0] == "bead_threshold":
                    all_thresholds.append(int(row_data[1]))

    # Take the median of the list
    sample_median = statistics.median(all_thresholds)
    return sample_median


# Read in the info from deconvolution params and find the fragment threshold
def readDeconvolutionParams(sample, inputdir):
    F = open(inputdir + sample + ".deconvolutionParams.csv")
    data = F.readlines()
    F.close()
    for i in range(len(data)):
        if data[i].find("bead_threshold,") != -1:
            thresh = data[i].strip().split(",")[1]

    return thresh


# Get the Median Unique Fragment Count Per Above Knee Bead
# Aggregate all above knee beads from all TIs and and find
# the median accross all the above knee beads
def readSampleBarcodeQuant(inputdir, thresh):
    barcode_quant_pattern = r".*\.barcodeQuantSimple.csv$"
    dirfiles = os.listdir(inputdir)
    # Get just the barcode quants files
    barcode_quant_files = getFileList(barcode_quant_pattern, dirfiles)

    above_knee_counts = []

    for quant_file in barcode_quant_files:
        with open(quant_file, "r") as file:
            for line in file:
                # The first item is the bead barcode, the second is the count
                row_data = line.strip().split(",")
                count = int(row_data[1])
                if count >= int(thresh):
                    above_knee_counts.append(count)

    all_TI_median = statistics.median(above_knee_counts)
    return all_TI_median


# Open the barcode quant file and calculate the median
# fragments in beads above the threshold
def readbarcodeQuant(sample, inputdir, thresh):
    F = open(inputdir + sample + ".barcodeQuantSimple.csv")
    data = F.readline()

    frags = []
    while data:
        count = int(data.strip().split(",")[1])
        # Only look at the beads above the threshold
        if count >= int(thresh):
            frags.append(count)
        data = F.readline()

    F.close()
    # Calculate the median fragment count
    medfrags = round(statistics.median(frags), 1)

    return medfrags


# Get Total Cells
# Sum total cells for each TI
# Get Median Unique Fragment Count Per Cell
# Aggregate all cells from all TIs and and find the median across all the cells
def readSampleQCStats(inputdir):
    QCstats_pattern = r".*\.QCstats.csv$"
    dirfiles = os.listdir(inputdir)
    # Get just the QCstats files
    QCstats_files = getFileList(QCstats_pattern, dirfiles)
    sample_cell_count = 0
    sample_frags = []
    for QCstats_file in QCstats_files:
        skip_line_count = 1
        with open(QCstats_file, "r") as file:
            for line in file:
                if skip_line_count > 0:
                    # Skip the header from each file
                    skip_line_count = skip_line_count - 1
                    # Read in the header
                    header = line.strip().split(",")
                    # Find the index of the unique nuclear fragments column
                    uniqnucfrags = header.index("uniqueNuclearFrags")
                else:
                    # Each line is a cell, so increment the cell count
                    sample_cell_count += 1
                    row_data = line.strip().split(",")
                    # Use the info we got from the header to get the frag data
                    sample_frags.append(int(row_data[uniqnucfrags]))
    sample_median_frags = statistics.median(sample_frags)
    return (sample_cell_count, sample_median_frags)


# Open the QC stats file and count the total cells and
# calculate the median fragments per cell
def readQCstats(sample, inputdir):
    F = open(inputdir + sample + ".QCstats.csv")
    # Read in the header
    header = F.readline().strip().split(",")
    # Find the index of the unique nuclear fragments column
    uniqnucfrags = header.index("uniqueNuclearFrags")

    data = F.readline()
    cells = 0
    frags = []
    while data:
        # Count the cells
        cells += 1
        # Grab the fragment count for each cell
        row = data.strip().split(",")
        frags.append(int(row[uniqnucfrags]))
        data = F.readline()
    # Calculate the median
    medfrags = round(statistics.median(frags), 1)

    F.close()

    return (cells, medfrags)


# Open the QC stats file and count the total cells and calculate
# the median fragments per cell
def readTIcounts(fastqTI, sample, inputdir):
    F = open(inputdir + "fastqTIreadcounts.csv")
    # Read in the header
    header = F.readline().strip().split(",")
    # Find the index of the unique nuclear fragments column
    sampleindex = header.index("sample")
    fastqindex = header.index("fastq")
    sequenceindex = header.index("sequence")
    readcount = header.index("count")

    data = F.readline()
    reads = 0
    while data:
        # Grab the fragment count for each cell
        row = data.strip().split(",")
        mysample = row[sampleindex]
        myfastqti = row[fastqindex] + "-" + row[sequenceindex]
        if mysample == sample:
            if myfastqti == fastqTI:
                reads = int(row[readcount])
        data = F.readline()

    F.close()

    return reads


# Open the QC stats file and count the total cells and
# calculate the median fragments per cell
def readSampleTIcounts(sample, inputdir):
    F = open(inputdir + "fastqTIreadcounts.csv")
    # Read in the header
    header = F.readline().strip().split(",")
    # Find the index of the unique nuclear fragments column
    sampleindex = header.index("sample")
    readcount = header.index("count")

    data = F.readline()
    reads = 0
    while data:
        # Grab the fragment count for each cell
        row = data.strip().split(",")
        mysample = row[sampleindex]
        if mysample == sample:
            reads += int(row[readcount])
        data = F.readline()

    F.close()

    return reads


# Executing script
if __name__ == "__main__":
    # Establish inputs
    args = parser.parse_args()
    names = grabNames(args.input_directory)

    # Grab the stats for each Fastq-TI pair
    stats = []
    for i in range(len(names)):
        mystats = [names[i]]
        thresh = readDeconvolutionParams(names[i], args.input_directory)
        mystats.append(thresh)
        mystats.append(str(readbarcodeQuant(names[i], args.input_directory, thresh)))
        cells, frags = readQCstats(names[i], args.input_directory)
        reads = readTIcounts(names[i], args.sample_id, args.input_directory)
        mystats.append(str(cells))
        mystats.append(str(round(reads / cells, 1)))
        mystats.append(str(frags))
        stats.append(mystats)

    # Write the Fastq-TI summary file
    F = open(args.input_directory + args.sample_id + ".beadFiltSummary.csv", "w")
    F.write(
        "FastqTI,UniqFragThresh,MedianBeadFrags,TotalCells,TotalReadPairsPerTotalCells,MedianCellFrags\n"  # noqa E501
    )
    for i in range(len(stats)):
        F.write(",".join(stats[i]) + "\n")
    F.close()

    # Make stats for the sample as a whole

    sample_thresh = getSampleThresh(args.input_directory)
    sample_med_frag_per_bead = readSampleBarcodeQuant(
        args.input_directory, sample_thresh
    )
    sample_cell_count, sample_med_frag_per_cell = readSampleQCStats(
        args.input_directory
    )
    sample_read_count = readSampleTIcounts(args.sample_id, args.input_directory)
    sample_stats = [
        args.sample_id,
        sample_thresh,
        sample_med_frag_per_bead,
        sample_cell_count,
        round(sample_read_count / sample_cell_count, 1),
        sample_med_frag_per_cell,
    ]
    formatted_sample_stats = [str(x) for x in sample_stats]

    with open(
        args.input_directory + args.sample_id + ".sampleBeadFiltSummary.csv", "w"
    ) as outfile:
        outfile.write(
            "Sample,UniqFragThresh,MedianBeadFrags,TotalCells,TotalReadPairsPerTotalCells,MedianCellFrags\n"  # noqa #501
        )
        outfile.write(",".join(formatted_sample_stats) + "\n")
