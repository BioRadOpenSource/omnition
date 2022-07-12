/*
Summarize alignment stats,
*/

params.options = [:]

process SUMMARIZE_ALIGNMENTS {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"

    if (workflow.profile == 'aws') {
        label 'large'
    } else {
        label 'cpu_medium'
        label 'memory_xsmall'
    }

    input:
    tuple val(sampleId), path(bam), path(index)
    path fasta
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}_flagstat.txt"), path("${sampleId}_idxstats.txt"),
      path("${sampleId}.alignment_summary_qc.txt"), emit:stats

    script:
    """
    # samtools flagstat and idxstats
    samtools flagstat ${bam} > ${sampleId}_flagstat.txt
    samtools idxstats ${bam} > ${sampleId}_idxstats.txt

    # Picard AlignmentSummaryMetrics
    picard CollectAlignmentSummaryMetrics \
      -R ${fasta} \
      -I ${bam} \
      -O ${sampleId}.alignment_summary_qc.txt \
      -Xmx${task.memory.toGiga()}g
    """
}
