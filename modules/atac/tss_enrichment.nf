/*
Calculate tss enrichment score
*/

params.options = [:]

process TSS_ENRICHMENT {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    if (workflow.profile == 'aws') {
        label 'medium'
    } else {
        label 'cpu_xsmall'
        label 'memory_xxsmall'
    }

    input:
    tuple val(sampleId), path(tss_matrix)
    val images_pulled

    output:
    tuple val(sampleId), path('*.tss_enrichment.csv'), emit: tss_enrichment

    script:
    """
    # Computing the TSS window as the entire range, plus the central TSS base
    FULL_WINDOW=\$(( ${params.options.tssWindowSize} + 1 ))

    atacTssEnrichment.py -w \${FULL_WINDOW} -m ${tss_matrix}
    """
}
