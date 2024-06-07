/*
Split fastqs by TI
*/

params.options = [:]

process COMPRESS_SPLIT_FASTQS {
    tag "${sampleId}-${index}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_small'

    input:
    tuple val(sampleId), val(index), path(fastq)
    val images_pulled

    output:
    tuple val("${sampleId}-${index}"), path("${sampleId}-*.complete_debarcoded.split.fastq.gz"), emit: split

    script:
    """
    cat ${sampleId}-${index}_R1.complete_debarcoded.split.fastq | \
        crabz -p ${task.cpus} > ${sampleId}-${index}_R1.complete_debarcoded.split.fastq.gz
    cat ${sampleId}-${index}_R2.complete_debarcoded.split.fastq | \
        crabz -p ${task.cpus} > ${sampleId}-${index}_R2.complete_debarcoded.split.fastq.gz
    """
}
