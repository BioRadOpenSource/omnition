/*
Merge read counts from reannotate bam back together
*/

process MERGE_REANN_READ_COUNTS {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_xxsmall'

    input:
    tuple val(sampleId), path(annotate_counts), path(bam_counts), path(frag_counts)
    val images_pulled

    output:
    tuple val(sampleId), path('*_merged_reannotate_read_counts.csv'), emit: count

    script:
    """
    publicAtacMergeReannCounts.py -i ./ -s ${sampleId}
    """
}
