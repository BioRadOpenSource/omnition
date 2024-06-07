/*
Check the read counts for each TI in a superloaded sample
*/

params.options = [:]

process CHECK_TI_COUNTS_SUPERLOADED {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    errorStrategy 'terminate'
    label 'cpu_xsmall'
    label 'memory_small'

    input:
    tuple val(sampleId), path(fastq)
    path asset
    val images_pulled

    output:
    path("*_reads.txt"), emit: split
    path("*_messages.txt"), emit: messages
    val true, emit:passed
    path("TIlen.txt"), emit: tilen

    script:
    """
    touch "${sampleId}_messages.txt"
    atacCountTiReadsSuperloaded.py -i ${fastq[0]} -l $asset -s $sampleId -m ${sampleId}_messages.txt
    """
}
