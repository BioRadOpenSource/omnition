/*
Calculate sequence saturation for each sample
*/

params.options = [:]

process SEQUENCE_SATURATION {
    tag "${sample_id}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    publishDir "${params.resultsDir}/${sample_id}/", mode: 'copy', overwrite: true

    if (workflow.profile == 'aws') {
        label 'large'
    } else {
        label 'cpu_medium'
        label 'memory_medium'
    }

    input:
        tuple val(sample_id), path(fragments), path(summary)
        val images_pulled

    output:
        tuple val(sample_id), path("*_sequence_saturation.csv"), emit: results

    script:
        """
        atacSequenceSaturation.R ${sample_id} ${fragments} ${summary}
        """
}
