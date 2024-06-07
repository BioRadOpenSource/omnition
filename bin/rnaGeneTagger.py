#!/usr/bin/env python3

import pysam
import click
import ray
import glob
import re
from ray.exceptions import ObjectStoreFullError
import sys

try:
    from coreBamCounter import getChromosomesWithAlignments
except Exception:
    print("Warning! This module is intended for use only in a Nextflow pipeline.")


@click.command()
@click.option("--bamfile", "-b", help="A .bam file with an index.")
@click.option("--gene-tag", "-gt", default="XT", help="Tag with gene feature.")
@click.option(
    "--feature-tag", "-ft", default="XF", help="Tag to write genic feature on."
)
@click.option(
    "--assignment-tag",
    "-at",
    default="XS",
    help="Tag in the .bam file that points to the assignment result.",
)
@click.option(
    "--unassigned-as-genes",
    default=True,
    help="In cases where GENE TAG is None, move the assignment tag to the gene tag.",
)
@click.option(
    "--tmp-dir", "-t", default="./tmp", help="Temp directory for split contig bams."
)
@click.option(
    "--out-dir", "-o", default="./", help="Output directory for final bam file."
)
@click.option(
    "--cpus", "-c", default=1, help="The number of CPUs for parallel processing."
)
@click.option(
    "--include-introns", is_flag=True, help="Allow intron alignments to count."
)
def main(
    bamfile,
    gene_tag,
    feature_tag,
    tmp_dir,
    out_dir,
    cpus,
    assignment_tag,
    unassigned_as_genes,
    include_introns,
):
    contigs = getChromosomesWithAlignments(bamfile)
    ray.init(num_cpus=cpus)
    try:
        ray.get(
            [
                tag_one_contig.remote(
                    bamfile,
                    contig,
                    tmp_dir,
                    gene_tag,
                    feature_tag,
                    assignment_tag,
                    unassigned_as_genes,
                    include_introns,
                )
                for contig in contigs
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
    cat_args = ["-o"] + [out_dir + "/tagged.bam"] + glob.glob(tmp_dir + "/*.bam")
    pysam.cat(*cat_args, catch_stout=False)
    ray.shutdown()


@ray.remote
def tag_one_contig(
    bamfile,
    contig,
    tmp_dir,
    gene_tag,
    feature_tag,
    assignment_tag,
    unassigned_as_genes,
    include_introns,
):
    """
    Iterates over a contig and moves barcodes and UMI from read name to desired tags
    Assumes the non-barcode/umi portion of the read name follows immediately
    """
    # fetch the contig
    bam = pysam.AlignmentFile(bamfile, "rb")
    Itr = bam.fetch(str(contig), multiple_iterators=True)
    # create a file to write outputs to
    splitBam = pysam.AlignmentFile(
        tmp_dir + "/" + str(contig) + ".bam", "wb", template=bam
    )
    # start processing reads
    for read in Itr:
        # check to see if the read has a gene feature annotated in its gene tag
        if read.has_tag(gene_tag):
            # get the gene tag
            tag = read.get_tag(gene_tag)
            # determine fate of read
            assign_read_tags(
                read, tag, include_introns, gene_tag, feature_tag, assignment_tag
            )
        # case where the read was unassigned and we want to move the reason
        # for it to the gene tag
        elif read.has_tag(gene_tag) is False and unassigned_as_genes:
            read.set_tag(gene_tag, read.get_tag(assignment_tag))
            # empty the assignment tag since it has now been duplicated on the gene tag
            read.set_tag(assignment_tag, None)
        # write alignment
        splitBam.write(read)
    splitBam.close()


def parse_ambiguous_alignment(tag, include_introns):
    tag = tag.split(",")
    if not include_introns:
        tag = [item for item in tag if not item.endswith("_INTRON")]
    genes = [re.sub("_EXON|_UTR|_INTRON", "", item) for item in tag]
    return (set(genes), set(tag))


def prioritize_exon(features):
    exons = set([feat for feat in features if feat.endswith("_EXON")])
    if len(exons) == 1:
        return list(exons)[0][: -len("_EXON")]
    else:
        return None


def assign_read_tags(read, tag, include_introns, gene_tag, feature_tag, assignment_tag):
    """
    Assign an alignment to a gene accounting for introns included or excluded.

    First, check for a comma in the gene_tag. If present, this indicates the read
    landed in more than one gene+feature.

    For multiple-feature alignments:
    1. If excluding introns, filter _INTRON features.
    2. If, in excluding introns, we removed all features, we know that the read
       aligned to multiple genes' introns and it is labeled as ambiguous.
    3. If there are multiple genes present, we check to see if one and only
       one of those genes was an alignment to an exon. If so, we 'prioritize'
       that exon and label the read for that gene. If there are multiple exons or
       there are multiple features considered as valid, the read is deemed ambiguous.
    4. If we align to one gene but multiple features within that gene, throw away
       any introns. If two features remain, an EXON+UTR condition exists and we
       classify the read as UTR, plus adjust the FeatureCounts XN tag to denote that
       the alignment was to a single gene. If one feature remains, the read is instead
       assigned to that feature.
    5. If the read is assigned to one feature from one gene, the appropriate tags are
       updated to reflect this.

    For single-feature alignments:
    1. Set a conditional filter for the acceptable feature types depending on whether
       we are or are not including introns.
    2. If the gene_tag ends with one of the acceptable feature types, set the gene
       tag to the gene name and the feature tag to the feature type.
    3. If the gene_tag does not end with an acceptable feature type, set it to
       INTERGENIC if the read is intergenic or, in the case of an introns-excluded
       analysis where the feature is an _INTRON, convert the gene name to INTRON.
       This matches how Omnition counts reads in intron-unaware mode.
    """
    if "," in tag:
        (genes, features) = parse_ambiguous_alignment(tag, include_introns)
        if len(genes) == 0:
            read.set_tag(gene_tag, "Unassigned_Ambiguous")
            read.set_tag(assignment_tag, None)
        elif len(genes) > 1:
            priority = prioritize_exon(features)
            if priority is not None:
                read.set_tag(gene_tag, priority)
                read.set_tag(feature_tag, "EXON")
            else:
                read.set_tag(gene_tag, "Unassigned_Ambiguous")
                read.set_tag(assignment_tag, None)
        elif len(genes) == 1 and len(features) > 1:
            features = set([feat for feat in features if not feat.endswith("_INTRON")])
            if len(features) > 1:
                feature = "UTR"
                if read.has_tag("XN"):
                    read.set_tag("XN", 1)
            else:
                feature = list(features)[0].split("_")[1]
            read.set_tag(gene_tag, list(genes)[0])
            read.set_tag(feature_tag, feature)
        elif len(genes) == 1 and len(features) == 1:
            feature = list(features)[0].split("_")[1]
            read.set_tag(gene_tag, list(genes)[0])
            read.set_tag(feature_tag, feature)
    else:
        if include_introns:
            condition = ("_EXON", "_UTR", "_INTRON")
        else:
            condition = ("_EXON", "_UTR")
        if tag.endswith(condition):
            read.set_tag(feature_tag, tag.split("_")[-1])
            read.set_tag(gene_tag, "".join(tag.split("_")[:-1]))
            read.set_tag(assignment_tag, None)
        elif tag == "INTERGENIC":
            read.set_tag(assignment_tag, None)
        elif tag.endswith("_INTRON"):
            read.set_tag(gene_tag, "INTRON")
            read.set_tag(assignment_tag, None)


if __name__ == "__main__":
    main()
