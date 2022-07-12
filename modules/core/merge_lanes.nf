/*
Merge lane split fastqs
*/

process MERGE_LANES {
    tag "${sampleId}"
    if (workflow.profile == 'aws') {
        label 'xlarge'
    } else {
        label 'cpu_xsmall'
        label 'memory_xxsmall'
    }

    input:
    tuple val(sampleId), path(forward), path(reverse)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}_*_001.complete.fastq.gz"), emit: complete

    script:
    """
    COUNT=\$(echo ${forward} ${reverse} | wc -w)

    if [ \$COUNT = 2 ];
    then
      mv ${forward} ${sampleId}_R1_001.complete.fastq.gz
      mv ${reverse} ${sampleId}_R2_001.complete.fastq.gz
    else
      cat ${forward} > ${sampleId}_R1_001.complete.fastq.gz
      cat ${reverse} > ${sampleId}_R2_001.complete.fastq.gz
    fi
    """
}
