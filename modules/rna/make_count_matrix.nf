/*
Aggregating read counts per gene per cell into a count matrix per sample
*/

params.options = [:]

process MAKE_COUNT_MATRIX {
    tag "${sampleId}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/countMatrix/", mode: 'copy', overwrite: true
    label 'cpu_small'
    label 'memory_small'

    input:
    tuple val(sampleId), path(count), path(symbols)
    val images_pulled

    output:
    tuple val(sampleId), path('*.unfiltered.mtx.gz'), path('*.unfiltered.barcodes.tsv'),
        path('*.unfiltered.genes.tsv'), emit: matrix

    script:
    """
    coreBuildCountMatrix.R -o ./ -p ${sampleId}.unfiltered -w NULL -g ${symbols} ${count}
    pigz -p ${task.cpus} ${sampleId}.unfiltered.mtx
    """
}
