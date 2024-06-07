/*
Check the read counts for each TI
*/

params.options = [:]

process CHECK_TI_COUNTS {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    errorStrategy 'terminate'
    label 'cpu_xsmall'
    label 'memory_small'

    input:
    tuple val(sampleId), path(fastq)
    path asset
    path config
    val images_pulled

    output:
    path("*_error.txt"), emit: error_files
    path("*_reads.txt"), emit: split
    path("*_messages.txt"), emit: messages

    script:
    """
    touch "${sampleId}_error.txt"
    touch "${sampleId}_messages.txt"
    atacCountTiReads.py -i ${fastq[0]} -l $asset -c $config -s $sampleId -m ${sampleId}_messages.txt
    """
}
