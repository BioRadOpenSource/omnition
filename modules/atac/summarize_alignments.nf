/*
Summarize alignment stats
*/

params.options = [:]

process SUMMARIZE_ALIGNMENTS {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_xlarge'

    input:
    tuple val(sampleId), path(bam), path(index)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}_flagstat.txt"), path("${sampleId}_idxstats.txt"),
        path("${sampleId}.alignment_summary_qc.txt"), emit:stats

    script:
    """
    # samtools flagstat and idxstats
    samtools flagstat -@ ${task.cpus} ${bam} > ${sampleId}_flagstat.txt
    samtools idxstats -@ ${task.cpus} ${bam} > ${sampleId}_idxstats.txt

    # calculate alignment metrics
    atacSummarizeAlignments.py -b ${bam} -c ${task.cpus} -s ${sampleId} -f ${sampleId}_flagstat.txt
    """
}
