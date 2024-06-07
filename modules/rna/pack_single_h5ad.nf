/*
Build per-sample H5AD file
*/

params.options = [:]

process PACK_SINGLE_H5AD {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/countMatrix/", mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(mtx), path(barcodes), path(features), path(summary), path(embedding)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}.h5ad"), emit: h5ad

    script:
    """
    rnaBuildH5ad.py -m ${mtx} -b ${barcodes} -f ${features} -d ${summary} -e ${embedding} -n ${sampleId}
    """
}
