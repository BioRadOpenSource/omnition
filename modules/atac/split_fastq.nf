/*
Split fastqs by TI
*/

params.options = [:]

process SPLIT_FASTQ {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_small'

    input:
    tuple val(sampleId), path(fastq), path(counts)
    val ti_checks_passed
    val images_pulled

    output:
    tuple val("${sampleId}"), path("${sampleId}-*.complete_debarcoded.split.fastq"), emit: split

    script:
    override_flag = params.options.tiErrorOverride == true ? "--override_errors" : ""
    """
    atacSplitFastq.py -r1 ${fastq[0]} \
        -r2 ${fastq[1]} \
        -s ${sampleId}  \
        -c ${counts} \
        -p ${task.cpus} \
        ${override_flag}
    """
}
