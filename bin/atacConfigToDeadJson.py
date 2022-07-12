#!/usr/bin/env python3

import argparse
import json


parser = argparse.ArgumentParser(
    description="Create a JSON config file with R1 or R2 TI sequences for DEAD."
)

# Defining argument group 'Options'
options = parser.add_argument_group("Options")
options.add_argument(
    "-i",
    "--index-list",
    required=True,
    help="Text file containing list of reference tagmentation indexes."
)
options.add_argument(
    "-m",
    "--max-dist",
    required=True,
    help="Maximum allowable distance (number of errros tolerated)."
)
options.add_argument(
    "-r",
    "--read",
    required=True,
    help="Read where TIs are present."
)
options.add_argument(
    "-t",
    "--i7asti",
    required=True,
    help="Whether or not the i7 sequence is used as TI.")


# Imports and parses the reference index list
def importTIindex(index):

    TIlist = list()
    TIlen = 0

    with open(index, "r") as f:
        data = f.readline()
        data = f.readline()
        while data:
            if len(data) > 0:
                data = data.strip().split(",")
                # Sequences are in second column
                TIlist.append(data[1])
                if len(data[1]) != TIlen:
                        if TIlen == 0:
                            TIlen = len(data[1])
                        else:
                            print(
                                "[ERROR] Reference TI sequences are of unequal lengths.  Please check pipeline config file."
                            )
                            exit(1)
            data = f.readline()

    return (TIlist,TIlen)

# Create the JSON file from the settings and indexes
def createR2JSON(indexlist,maxdist,outfile):

    # Set up dictionary for JSON
    jsondict=dict()
    jsondict["primary_fastq"]=0
    jsondict["sequences"]=list()

    # Set up sequence entry for TIs
    sequencedict=dict()
    sequencedict["type"]="Barcode"
    sequencedict["max_dist"]=int(maxdist)
    sequencedict["values"]=indexlist
    jsondict["sequences"].append(sequencedict)

    # Create JSON
    with open(outfile, 'w', encoding='utf-8') as f:
        json.dump(jsondict, f, indent=4)

# Create the JSON file from the settings and indexes
def editR1JSON(indexlist,maxdist,outfile):

    with open('atac.json') as f:
        atacjson = json.load(f)

    # Set up sequence entry for TIs
    sequencedict=dict()
    sequencedict["type"]="Barcode"
    sequencedict["max_dist"]=int(maxdist)
    sequencedict["values"]=indexlist

    atacjson['sequences'].insert(len(atacjson['sequences'])-1,sequencedict)

    # Create JSON
    with open(outfile, 'w', encoding='utf-8') as f:
        json.dump(atacjson, f, indent=4)


# Executing script
if __name__ == "__main__":

    # Establish inputs
    args = parser.parse_args()
    outfile = "TI.json"
    lenfile = "TIlen.txt"

    # Import and parse TIs
    indexlist,indexlen = importTIindex(args.index_list)

    # Write out length of TIs
    with open(lenfile, 'w') as f:
        f.write(str(indexlen))

    # Create JSON outfile
    if args.read == "r2" or args.i7asti == "true":
        createR2JSON(indexlist,args.max_dist,outfile)
    else:
        editR1JSON(indexlist,args.max_dist,outfile)

