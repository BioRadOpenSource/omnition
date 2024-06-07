/*
Calculating various alignment metrics
*/

params.options = [:]

process PICARD {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(bam), path(refflat), path(interval_list)
    val images_pulled

    output:
    tuple val(sampleId), path('*.rna_seq_metrics.txt'), emit: metrics

    script:
    """
    # Reducing heap size to 80% of allocated resources account for Java overhead
    JAVA_MEMORY="\$(((${task.memory.toGiga()} * 4)/ 5))g"

    picard CollectRnaSeqMetrics \
        INPUT=${bam} \
        OUTPUT=${sampleId}.rna_seq_metrics.txt \
        REF_FLAT=${refflat} \
        STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND \
        RIBOSOMAL_INTERVALS=${interval_list} \
        -Xmx\$JAVA_MEMORY
    """
}
