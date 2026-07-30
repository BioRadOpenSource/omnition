#!/usr/bin/env python3

import click


# Process the filtered count matrix to calculate average antibody counts per cell
def processFilteredMatrix(filename):
    antibodies = []

    F = open(filename)
    # Read header line to get antibody names
    header = F.readline()
    header = header.strip().replace('"', "").split(",")
    antibodies = header[1:]  # Skip the first column (cell barcodes)

    filteredcounts = [0 for x in range(len(antibodies))]
    totalcells = 0

    # Process each cell (row)
    data = F.readline()
    while data:
        data = data.strip().split(",")
        if len(data) > 1:  # Make sure we have data beyond just the barcode
            totalcells += 1
            for i in range(len(data[1:])):  # Skip the first column (barcode)
                filteredcounts[i] += int(data[1 + i])
        data = F.readline()
    F.close()

    return (antibodies, filteredcounts, totalcells)


# Output the counts per cell
def writeoutfile(sample, antibodies, anticounts, totalcells):
    antistring = "Sample,Cells"
    for j in range(len(antibodies)):
        antistring += "," + antibodies[j]
    antistring += "\n"

    outfile = sample + "_antibody_counts_per_cell.csv"

    F = open(outfile, "w")
    F.write(antistring)

    fprint = sample + "," + str(totalcells)
    for j in range(len(anticounts)):
        fprint += "," + str(round(anticounts[j] / totalcells, 1))
    fprint += "\n"
    F.write(fprint)

    F.close()


# Main code
@click.command()
@click.option(
    "-c", "--count-matrix",
    required=True,
    help="CSV file of the filtered cell counts"
)
@click.option(
    "-s", "--sample-id",
    required=True,
    help="Prefix for the output files"
)
def main(count_matrix, sample_id):
    """Calculate the average antibody counts per cell from a filtered count matrix."""
    antibodies, anticounts, totalcells = processFilteredMatrix(count_matrix)
    writeoutfile(sample_id, antibodies, anticounts, totalcells)


if __name__ == "__main__":
    main()
