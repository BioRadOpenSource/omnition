/*
Calculate sequence saturation for each sample
*/

params.options = [:]

process SEQUENCE_SATURATION {
    tag "${sampleId}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/", mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_medium'

    input:
        tuple val(sampleId), path(fragments), path(summary)
        val images_pulled

    output:
        tuple val(sampleId), path("*_sequence_saturation.csv"), emit: results

    script:
        """
        atacSequenceSaturation.R ${sampleId} ${fragments} ${summary}
        """
}
