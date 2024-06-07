#!/usr/bin/env python3

import anndata as ad
import click
import os
import glob


@click.command()
@click.option("--input-dir", "-i", required=True, help="A directory of .h5ad files")
def main(input_dir):
    # Get all h5ad files
    files = gather_h5ad_files(input_dir)

    # Get all sample names
    samples = gather_sample_ids(files)

    # Build h5ad dictonary
    dic = build_h5ad_dict(files)

    # Merge h5ad dictonary
    all_samples = merge_h5ad_samples(dic, samples)

    # Write merged h5ad
    all_samples.write("all_samples.h5ad", compression="gzip")


def gather_h5ad_files(input_dir):
    """
    Gather all sample .h5ad files to merge
    """
    # Glob all .h5ad files
    files = [
        os.path.realpath(file) for file in glob.glob("{}/*.h5ad".format(input_dir))
    ]
    return files


def gather_sample_ids(files):
    """
    Use list of .h5ad files to gather sample names.
    """
    # Get list of samples
    samples = [os.path.basename(file).strip(".h5ad") for file in files]
    return samples


def build_h5ad_dict(files):
    """
    Iterate over .h5ad files and append sample and file to dict
    """
    # Initialize output dict
    out = {}

    # Build dict
    for file in files:
        sample_id = os.path.basename(file).strip(".h5ad")
        out[sample_id] = ad.read_h5ad(file)
    return out


def merge_h5ad_samples(dic, samples):
    """
    Iterate over and merge together all found .h5ad files
    """
    # Get list of all samples to iterate across
    all_samples = dic.get(samples[0])

    # Concatenate all h5ad information together
    all_samples = all_samples.concatenate(
        [dic.get(sample) for sample in samples[1:]],
        batch_key="Sample ID",
        batch_categories=samples,
    )
    return all_samples


if __name__ == "__main__":
    main()
