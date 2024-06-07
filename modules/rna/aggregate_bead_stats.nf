/*
Aggregate bead statistics
*/

params.options = [:]

process AGGREGATE_BEAD_STATS {
    tag "${sampleId}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(sequence_counts), path(allowlist), path(cell_expression), path(bead_expression),
        path(barcode_translate)
    val images_pulled

    output:
    tuple val(sampleId), path('*.bead_summary.csv'), emit: bead_summary

    script:
    """
    rnaBeadStats.R -bc ${sequence_counts} \
        -a ${allowlist} \
        -c ${cell_expression} \
        -b ${bead_expression} \
        -bt ${barcode_translate}
    """
}
