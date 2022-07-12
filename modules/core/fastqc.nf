/*
Assessing read quality
*/

params.options = [:]

process FASTQC {
    tag "${sample_id}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    if (workflow.profile == 'aws') {
        label 'medium'
    } else {
        label 'cpu_small'
        label 'memory_small'
    }

    input:
    tuple val(sample_id), path(read_pair)
    val images_pulled

    output:
    tuple val(sample_id), path('*_fastqc.zip'), emit: zip

    script:
    """
    fastqc ${read_pair} --threads ${task.cpus} --quiet  2>&1
    """
}
