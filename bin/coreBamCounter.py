#!/usr/bin/env python3

import pysam
import ray
import pandas as pd
from ray.exceptions import ObjectStoreFullError
from pathlib import Path
from pathlib import PurePath
from functools import reduce
import click
import sys


@click.command()
@click.option("--input", "-i", help="A .bam file with an index.")
@click.option(
    "--barcode-tag",
    "-bt",
    default="XB",
    help="Tag in the .bam file that points to the barcode.",
)
@click.option(
    "--umi-tag",
    "-ut",
    default="XU",
    help="Tag in the .bam file that points to the UMI.",
)
@click.option(
    "--gene-tag",
    "-gt",
    default="XT",
    help="Tag in the .bam file that points to the gene symbol.",
)
@click.option(
    "--cpus", "-c", default=1, help="The number of CPUs for parallel processing."
)
@click.option(
    "--mapq",
    "-m",
    default=0,
    help="Include alignments at or above this mapping quality.",
)
@click.option(
    "--out-path",
    "-o",
    default="",
    help="Output path; defaults to location of .bam file input.",
)
def main(input, barcode_tag, umi_tag, gene_tag, cpus, out_path, minmapq):
    chroms2use = getChromosomesWithAlignments(input)
    ray.init(num_cpus=cpus)
    try:
        results = ray.get(
            [
                countUniqueGenicUmi.remote(input, chrom, barcode_tag, umi_tag, gene_tag)
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
    d = reduce(lambda x, y: pd.merge(x, y, on="barcode", how="outer"), results).fillna(
        value=0
    )
    d["N"] = d.drop("barcode", axis=1).sum(axis=1).astype("int64")
    d = d[["barcode", "N"]]
    if out_path == "":
        out_dir = Path(input).parents[0]
    else:
        out_dir = Path(out_path)
    out_file = PurePath(input).name.replace("bam", "sequence_counts.csv")
    out_dir.mkdir(parents=True, exist_ok=True)
    d.to_csv(out_dir / out_file, index=False)
    ray.shutdown()


# returns list of chromosomes with any mapping alignments
def getChromosomesWithAlignments(bam):
    idx = pysam.idxstats(bam).splitlines()
    chroms2use = []
    for chr in range(len(idx)):
        idxstat = idx[chr].split("\t")
        if int(idxstat[2]) + int(idxstat[3]) > 0:
            chroms2use.append(idxstat[0])
    return chroms2use


@ray.remote
def countUniqueGenicUmi(bam, chrom, barcode_tag, umi_tag, gene_tag):
    """Counts unique genic UMI"""
    bamFile = pysam.AlignmentFile(bam, "rb")
    Itr = bamFile.fetch(chrom, multiple_iterators=True, until_eof=False)
    dic = {}
    out = {}
    for read in Itr:
        if read.has_tag(barcode_tag) & read.has_tag(umi_tag) & read.has_tag(gene_tag):
            barc = read.get_tag(barcode_tag)
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
    genicUmiCount = sorted(out.items(), key=lambda kv: kv[1], reverse=True)
    genicUmiCount = pd.DataFrame(
        genicUmiCount, columns=["barcode", "uniqueGenicUmiCount"]
    )
    return genicUmiCount


if __name__ == "__main__":
    main()
