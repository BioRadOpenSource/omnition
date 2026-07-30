/*
Calculate tss enrichment score
*/

process TSS_ENRICHMENT {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xxsmall'

    input:
    tuple val(sampleId), path(tss_matrix)
    val images_pulled

    output:
    tuple val(sampleId), path('*.tss_enrichment.csv'), emit: tss_enrichment
    path '*.tss_metric.csv', emit: tss_metric

    script:
    """
    # Computing the TSS window as the entire range, plus the central TSS base
    FULL_WINDOW=\$(( ${params.options.tssWindowSize} + 1 ))

    publicAtacTssEnrichment.py -w \${FULL_WINDOW} -m ${tss_matrix}
    """
}
