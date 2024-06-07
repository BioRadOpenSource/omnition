/*
Compile bams by sample
*/

params.options = [:]

process COMPILE_ALIGNMENTS {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/alignments",
        pattern: '*{.final.bam,.final.bam.bai,.tagged.duplicatesmarked.bam,.tagged.duplicatesmarked.bam.bai}',
        mode: 'copy', overwrite: true
    label 'cpu_large'
    label 'memory_xsmall'

    input:
    tuple val(sampleId), path(bam), path(tagged_bam), path(marked_bam)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}.final.bam"), path("${sampleId}.final.bam.bai"), emit:bam
    tuple val(sampleId), path("${sampleId}.alignments.tagged.bam"),
        path("${sampleId}.alignments.tagged.bam.bai"), emit:tagged_bam
    tuple val(sampleId), path("${sampleId}.alignments.tagged.duplicatesmarked.bam"),
        path("${sampleId}.alignments.tagged.duplicatesmarked.bam.bai"), emit:marked_bam

    script:
    """
    (samtools merge -@ ${task.cpus} ${sampleId}.final.bam ${bam}) &> ${sampleId}_finalmerge.log
    samtools index -@ ${task.cpus} ${sampleId}.final.bam
    (samtools merge -@ ${task.cpus} ${sampleId}.alignments.tagged.bam ${tagged_bam}) &> ${sampleId}_taggedmerge.log
    samtools index -@ ${task.cpus} ${sampleId}.alignments.tagged.bam
    (samtools merge -@ ${task.cpus} ${sampleId}.alignments.tagged.duplicatesmarked.bam \
    ${marked_bam}) &> ${sampleId}_markedmerge.log
    samtools index -@ ${task.cpus} ${sampleId}.alignments.tagged.duplicatesmarked.bam
    """
}
