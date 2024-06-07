/*
Tagging BAM with gene annotations
*/

params.options = [:]

process TAG_GENES {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    label 'cpu_large'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(bam), path(index), path(saf)
    val images_pulled

    output:
    tuple val(sampleId), path('*.temp.bam'), path('*.temp.bam.bai'), emit: bam
    path '*.summary', emit: summary

    script:
    """
    featureCounts \
        -T ${task.cpus} \
        --primary \
        -M \
        -s 1 \
        -Q 1 \
        -O \
        -a ${saf} \
        -o ${sampleId}.gene_counts \
        -R BAM ${bam} \
        -F SAF \
        --fracOverlap 0.80

    samtools sort -l 2 -@ ${task.cpus} ${bam}.featureCounts.bam -o ${sampleId}.temp.bam
    samtools index -@ ${task.cpus} ${sampleId}.temp.bam
    """
}
