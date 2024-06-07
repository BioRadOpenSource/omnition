/*
Aligning reads to the reference genome(s)
*/

params.options = [:]

process STAR_ALIGN {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    label 'cpu_xlarge'
    label 'memory_large'

    input:
    tuple val(sampleId), path(fastq), path(index)
    val images_pulled

    output:
    tuple val(sampleId), path('*_Aligned.out.bam'), emit: raw_bam
    tuple val(sampleId), path('*Log.final.out'), emit: log

    script:
    """
    STAR \
        --runThreadN ${task.cpus} \
        --limitIObufferSize 60000000 \
        --limitOutSJcollapsed 2000000 \
        --genomeDir ${index} \
        --readFilesIn ${fastq[1]} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${sampleId}_ \
        --outSAMunmapped Within \
        --outSAMtype BAM Unsorted \
        --outBAMcompression 2
    """
}
