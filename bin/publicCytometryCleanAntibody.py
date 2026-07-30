#! /usr/bin/env python3

import argparse
import csv
import re

parser = argparse.ArgumentParser(
    description="Remove disallowed characters from antibody file."
)
parser.add_argument(
    "-f", "--file", type=str, required=True, dest="file", help="Path to antibody file."
)

args = parser.parse_args()

cleaned_data = []
reverse_cols = False
with open(args.file, "r") as file:
    data = csv.reader(file)
    # Check if the antibody file is in label,tag or tag,label order
    antibody_file_list = list(data)
    first_string = antibody_file_list[0][0]
    DNA_pattern = re.compile("[ATCGatcg]+")
    if re.fullmatch(DNA_pattern, first_string):
        # columns are tag,label
        reverse_cols = False
    else:
        reverse_cols = True

    # Replace bad chars with dashes and execute column reversal if necessary
    for row in antibody_file_list:
        if reverse_cols:
            # columns are tag,label
            clean_antibody_label = row[0].replace("/", "-")
            antibody_tag = row[1]
        else:
            # columns are label,tag
            clean_antibody_label = row[1].replace("/", "-")
            antibody_tag = row[0]

        clean_row = [antibody_tag, clean_antibody_label]
        cleaned_data.append(clean_row)

with open(args.file, "w") as file:
    outwriter = csv.writer(file, lineterminator="\n")
    for row in cleaned_data:
        outwriter.writerow(row)
