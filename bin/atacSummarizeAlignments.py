#!/usr/bin/env python3

import pysam
import click
import ray
import re
import pandas as pd
from ray.exceptions import ObjectStoreFullError
import sys

try:
    from coreBamCounter import getChromosomesWithAlignments
except Exception:
    print("Warning! This module is inteded for use only in a Nextflow pipeline.")


@click.command()
@click.option(
    "--bam_file", "-b", type=click.Path(), required=True, help="Path to bam file."
)
@click.option(
    "--cpus", "-c", default=1, help="The number of CPUs to use in parallel processing."
)
@click.option(
    "--flagstat_file",
    "-f",
    type=click.Path(),
    required=True,
    help="Path to flagstat file for input bam file.",
)
@click.option("--sample_id", "-s", default="MySample", help="The sample ID.")
def main(bam_file, cpus, sample_id, flagstat_file):
    contigs = getChromosomesWithAlignments(bam_file)

    ray.init(num_cpus=cpus)
    try:
        results = ray.get([get_stats.remote(bam_file, contig) for contig in contigs])
    except ObjectStoreFullError:
        ray.shutdown()
        print("OOM error: ObjectStoreFullError")
        sys.exit(137)
    except Exception as err:
        ray.shutdown()
        print(f"Unexpected {err=}, {type(err)=}")
        sys.exit(1)
    ray.shutdown()
    # parse flagstat file into dictionary
    flagstat_dict = parse_flagstat(flagstat_file)

    # create AlignmentStatistics
    stats = AlignmentStatistics(flagstat_dict, results)

    # generate a data frame
    df = pd.DataFrame(data=stats.stats_dict())

    # write data frame to file
    out_file = "{}.{}".format(sample_id, "alignment_summary_qc.txt")
    with open(out_file, "w") as f:
        # All the newlines are for compatiblity with downstream parsers!
        f.write("## ALIGNMENT SUMMARY METRICS\n\n\n\n\n\n")

    df.to_csv(out_file, sep="\t", mode="a", index=False, header=True)


class AlignmentStatistics:
    def __init__(self, flagstat_dict, alignment_counts):
        self.FIRST_OF_PAIR_TOTAL_READS = flagstat_dict.get("read1")
        self.SECOND_OF_PAIR_TOTAL_READS = flagstat_dict.get("read2")
        self.PAIR_TOTAL_READS = flagstat_dict.get("total") - flagstat_dict.get(
            "secondary"
        )
        self.PF_READS_ALIGNED_FIRST_OF_PAIR = sum(
            [result[0] for result in alignment_counts]
        )
        self.PF_READS_ALIGNED_SECOND_OF_PAIR = sum(
            [result[1] for result in alignment_counts]
        )
        self.PF_READS_ALIGNED_PAIR = (
            self.PF_READS_ALIGNED_FIRST_OF_PAIR + self.PF_READS_ALIGNED_SECOND_OF_PAIR
        )
        self.READS_ALIGNED_IN_PAIRS = flagstat_dict.get("with itself and mate mapped")
        self.READS_ALIGNED_IN_PAIRS_FIRST_OF_PAIR = self.READS_ALIGNED_IN_PAIRS / 2
        self.READS_ALIGNED_IN_PAIRS_SECOND_OF_PAIR = int(
            self.READS_ALIGNED_IN_PAIRS_FIRST_OF_PAIR
        )
        self.PROPER_PAIRS = flagstat_dict.get("properly paired")
        self.PF_READS_IMPROPER_PAIRS_FIRST_OF_PAIR = int(
            self.PF_READS_ALIGNED_FIRST_OF_PAIR - (self.PROPER_PAIRS / 2)
        )
        self.PF_READS_IMPROPER_PAIRS_SECOND_OF_PAIR = int(
            self.PF_READS_ALIGNED_SECOND_OF_PAIR - (self.PROPER_PAIRS / 2)
        )
        self.PF_READS_IMPROPER_PAIRS_PAIR = int(
            self.PF_READS_IMPROPER_PAIRS_FIRST_OF_PAIR
            + self.PF_READS_IMPROPER_PAIRS_SECOND_OF_PAIR
        )

    def stats_dict(self):
        d = {
            "CATEGORY": ["FIRST_OF_PAIR", "SECOND_OF_PAIR", "PAIR"],
            "TOTAL_READS": [
                self.FIRST_OF_PAIR_TOTAL_READS,
                self.SECOND_OF_PAIR_TOTAL_READS,
                self.PAIR_TOTAL_READS,
            ],
            "PF_READS_ALIGNED": [
                self.PF_READS_ALIGNED_FIRST_OF_PAIR,
                self.PF_READS_ALIGNED_SECOND_OF_PAIR,
                self.PF_READS_ALIGNED_PAIR,
            ],
            "READS_ALIGNED_IN_PAIRS": [
                self.READS_ALIGNED_IN_PAIRS_FIRST_OF_PAIR,
                self.READS_ALIGNED_IN_PAIRS_SECOND_OF_PAIR,
                self.READS_ALIGNED_IN_PAIRS,
            ],
            "PF_READS_IMPROPER_PAIRS": [
                self.PF_READS_IMPROPER_PAIRS_FIRST_OF_PAIR,
                self.PF_READS_IMPROPER_PAIRS_SECOND_OF_PAIR,
                self.PF_READS_IMPROPER_PAIRS_PAIR,
            ],
        }
        return d


def parse_flagstat(flagstat_file):
    """
    Parse a flagstat file into a dictionary
    Remove non-qc-passed and extra text
    """

    # regexes for parsing flagstat files
    flagstat_regexes = {
        "total": r"(\d+) \+ (\d+) in total \(QC-passed reads \+ QC-failed reads\)",
        "secondary": r"(\d+) \+ (\d+) secondary",
        "supplementary": r"(\d+) \+ (\d+) supplementary",
        "duplicates": r"(\d+) \+ (\d+) duplicates",
        "mapped": r"(\d+) \+ (\d+) mapped \((.+):(.+)\)",
        "paired in sequencing": r"(\d+) \+ (\d+) paired in sequencing",
        "read1": r"(\d+) \+ (\d+) read1",
        "read2": r"(\d+) \+ (\d+) read2",
        "properly paired": r"(\d+) \+ (\d+) properly paired \((.+):(.+)\)",
        "with itself and mate mapped": r"(\d+) \+ (\d+) with itself and mate mapped",
        "singletons": r"(\d+) \+ (\d+) singletons \((.+):(.+)\)",
        "with mate mapped to a different chr": r"(\d+) \+ (\d+) with mate mapped to a different chr",  # noqa E501
        "with mate mapped to a different chr (mapQ >= 5)": r"(\d+) \+ (\d+) with mate mapped to a different chr \(mapQ>=5\)",  # noqa E501
    }
    with open(flagstat_file) as f:
        flagstat = f.read().splitlines()
    flagstat_dict = {}
    for name, regex in flagstat_regexes.items():
        r = re.compile(regex)
        result = list(filter(r.match, flagstat))
        if result:
            entry = int(result[0].split("+")[0].strip())
            flagstat_dict[name] = entry
    return flagstat_dict


@ray.remote
def get_stats(bam_file, contig):
    """
    Count individual R1 and R2 primary alignments
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    reads = bam.fetch(region=contig)
    r1_count = 0
    r2_count = 0
    for read in reads:
        if read.is_secondary or read.is_unmapped:
            continue
        if read.is_read1:
            r1_count += 1
        else:
            r2_count += 1
    return (r1_count, r2_count)


if __name__ == "__main__":
    main()
