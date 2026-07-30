#!/usr/bin/env python3

import anndata as ad
import mudata as mu
import click
import numpy as np


@click.command()
@click.option('--rna', required=True, help='Merged RNA .h5ad file (from rnaMergeH5ad)')
@click.option('--adt', required=True, help='Merged ADT .h5ad file (from rnaMergeH5ad)')
@click.option('--output', default='multiomics.h5mu', help='Output filename for combined multiomics MuData (.h5mu)')
def main(rna, adt, output):
    # Read merged RNA and ADT files
    rna_adata = ad.read_h5ad(rna)
    adt_adata = ad.read_h5ad(adt)

    # Check that barcodes/cells match
    if not np.array_equal(rna_adata.obs.index, adt_adata.obs.index):
        # Harmonize by intersection and sorting
        common = sorted(list(set(rna_adata.obs.index) & set(adt_adata.obs.index)))
        if len(common) == 0:
            raise ValueError('No shared barcodes between RNA and ADT files.')
        print(f"Reordering and subsetting to {len(common)} shared barcodes.")
        rna_adata = rna_adata[common, :].copy()
        adt_adata = adt_adata[common, :].copy()

    # Tag each modality with assay type
    rna_adata.obs['assay'] = 'RNA'
    rna_adata.var['assay'] = 'RNA'
    adt_adata.obs['assay'] = 'ADT'
    adt_adata.var['assay'] = 'ADT'

    # Create MuData object with both modalities
    mdata = mu.MuData({'rna': rna_adata, 'adt': adt_adata})
    # Update global observations with shared metadata
    # MuData automatically creates .obs from intersection of modality .obs
    mdata.update()

    # Save combined MuData
    mdata.write(output)
    print(f'Combined multiomics MuData written to {output}')
    print(f'RNA modality: {mdata.mod["rna"].n_obs} cells, {mdata.mod["rna"].n_vars} genes')
    print(f'ADT modality: {mdata.mod["adt"].n_obs} cells, {mdata.mod["adt"].n_vars} antibodies')


if __name__ == '__main__':
    main()
