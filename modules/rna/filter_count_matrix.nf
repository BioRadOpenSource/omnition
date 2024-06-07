/*
Filter count matrix to include only barcodes on the allowlist
*/

params.options = [:]

process FILTER_COUNT_MATRIX {
    tag "${sampleId}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/countMatrix/", mode: 'copy', overwrite: true
    label 'cpu_small'
    label 'memory_small'

    input:
    tuple val(sampleId), path(count), path(allowlist), path(symbols)
    val images_pulled

    output:
    tuple val(sampleId), path('*.filtered.mtx.gz'), path('*.filtered.barcodes.tsv'),
        path('*.filtered.genes.tsv'), emit: matrix

    script:
    """
    coreBuildCountMatrix.R -o ./ -p ${sampleId}.filtered -w ${allowlist} -g ${symbols} ${count}
    pigz -p ${task.cpus} ${sampleId}.filtered.mtx
    """
}
