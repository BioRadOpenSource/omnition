/*
Caluclate insert size metrics
*/

params.options = [:]

process CALCULATE_INSERT_SIZE_METRICS {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    if (workflow.profile == 'aws') {
        label 'xlarge'
    } else {
        label 'cpu_small'
        label 'memory_xsmall'
    }

    input:
    tuple val(sampleId), path(bam), path(index)
    val images_pulled

    output:
    tuple val(sampleId), path('*.insert_size_metrics.txt'), emit: metrics

    script:
    """
    picard CollectInsertSizeMetrics \
      -I ${bam} \
      -O ${sampleId}.insert_size_metrics.txt \
      -H ${sampleId}.insert_size_histogram.pdf \
      --VALIDATION_STRINGENCY SILENT \
      --TMP_DIR ./ \
      -Xmx${task.memory.toGiga()}g
    """
}
