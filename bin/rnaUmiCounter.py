#!/usr/bin/env python3

# nested dictionary for counting barcodes and umis

import pysam
import ray
import pandas as pd
from ray.exceptions import ObjectStoreFullError
from functools import reduce
import click
import sys


@click.command()
@click.option("--bam", "-b", help="A .bam file with an index.")
@click.option(
    "--barcode-tag", "-bt", default="XC", help="Tag to use for the cell barcode."
)
@click.option("--umi-tag", "-ut", default="XM", help="Tag to use for the UMI.")
@click.option(
    "--cpus", "-c", default=1, help="Number of CPUs to use in parallel processing."
)
def main(bam, barcode_tag, umi_tag, cpus):
    # get chromosomes that have reads
    idx = pysam.idxstats(bam).splitlines()
    chroms2use = []
    for chr in range(len(idx)):
        idxstat = idx[chr].split("\t")
        if int(idxstat[2]) + int(idxstat[3]) > 0:
            chroms2use.append(idxstat[0])

    # initialize parallel workers
    ray.init(num_cpus=cpus)

    # do the work
    try:
        results = ray.get(
            [
                count_umis.remote(bam, chrom, barcode_tag, umi_tag)
                for chrom in chroms2use
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
    # merge, sort
    d = (
        reduce(lambda x, y: pd.merge(x, y, on="barcode", how="outer"), results)
        .fillna(value=0)
        .astype("int")
    )

    # write out
    d.to_csv(bam.replace("bam", "sequence_counts.csv"))


@ray.remote
def count_umis(bam, chrom, barcode_tag, umi_tag):
    """
    Iterate over records and store UMI counts
    """
    bamFile = pysam.AlignmentFile(bam, "rb")
    Itr = bamFile.fetch(chrom, multiple_iterators=True, until_eof=False)
    dic = {}
    out = {}
    readCounts = {}
    for read in Itr:
        if read.has_tag(barcode_tag) & read.has_tag(umi_tag):
            barc = read.get_tag("XC")
            if readCounts.get(barc) is None:
                readCounts[barc] = 1
            else:
                readCounts[barc] += 1
            umi = read.get_tag(umi_tag)
            if dic.get(barc) is None:
                dic[barc] = {}
                dic[barc][umi] = 1
            elif dic[barc].get(umi) is None:
                dic[barc][umi] = 1
            else:
                dic[barc][umi] += 1
    for key in dic.keys():
        out[key] = len(dic[key])
    umiCount = sorted(out.items(), key=lambda kv: kv[1], reverse=True)
    readCount = sorted(readCounts.items(), key=lambda kv: kv[1], reverse=True)
    umiCount = pd.DataFrame(umiCount, columns=["barcode", "umiCount." + chrom])
    readCount = pd.DataFrame(readCount, columns=["barcode", "readCount." + chrom])
    readCount = readCount.set_index("barcode").join(umiCount.set_index("barcode"))
    return readCount


if __name__ == "__main__":
    main()
