/*
Mark duplicate reads
*/

params.options = [:]

process ID_DUPLICATES {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    label 'cpu_xlarge'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(bam), path(index), path(allowlist)
    val images_pulled

    output:
    tuple val(sampleId), path("*.duplicate_counts.csv"), emit:duplicate_count

    script:
    include_introns = params.options.includeIntrons == true ? "--include-introns" : ""
    """
    # Count duplicates of the Bam file
    omnition rna generate-bam-duplicates -b ${bam} \
        -a ${allowlist} \
        ${include_introns}

    mv duplicate_counts.csv ${sampleId}.duplicate_counts.csv
    """
}
