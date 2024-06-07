/*
Concatenate all h5ad files, dim reduction, clustering
*/

params.options = [:]

process PACK_BATCH_H5AD {
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    publishDir "${params.options.reportsDir}/", mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_medium'

    input:
    path h5ad
    val images_pulled

    output:
    path 'all_samples.h5ad', emit: h5ad

    script:
    """
    rnaMergeH5ad.py -i ./
    """
}
