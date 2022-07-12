/*
Mark duplicates
*/

params.options = [:]

process MARK_DUPLICATES {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    publishDir "${params.resultsDir}/${sampleId}/alignments", pattern:'*.{bam,bai}',
        mode: 'copy', enabled: !params.options.barcodedTn5, overwrite: true

    if (workflow.profile == 'aws') {
        label 'large'
    } else {
        label 'cpu_medium'
        label 'memory_xlarge'
    }

    input:
    tuple val(sampleId), path(bam), path(index)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}.alignments.tagged.duplicatesmarked.bam"),
        path("${sampleId}.alignments.tagged.duplicatesmarked.bam.bai"), path("${sampleId}.dedup_stats.txt"), emit: bam
    tuple val(sampleId), path("${sampleId}.dedup_stats.txt"), emit: stats

    script:
    """
    # Reducing heap size to account for Java overhead
    JAVA_MEMORY="\$((${task.memory.toGiga()} / 2))g"

    # Making TMP directory for overflow
    mkdir tmp/

    # Identifying duplicates
    picard \
      -Xms\$JAVA_MEMORY \
      -Xmx\$JAVA_MEMORY \
      MarkDuplicates \
      --TMP_DIR tmp/ \
      --BARCODE_TAG XB \
      --TAG_DUPLICATE_SET_MEMBERS true \
      --SORTING_COLLECTION_SIZE_RATIO ${params.options.sortSize} \
      -I ${bam} \
      -O ${sampleId}.alignments.tagged.duplicatesmarked.bam \
      -M ${sampleId}.dedup_stats.txt

    # Indexing output BAM
    samtools index -@ ${task.cpus} ${sampleId}.alignments.tagged.duplicatesmarked.bam
    """
}
