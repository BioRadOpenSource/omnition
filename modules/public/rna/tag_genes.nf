/*
Tagging BAM with gene annotations
*/

process TAG_GENES {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
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

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=7

    # Check that we have expected number of sample input files
    check_input_files "TAG_GENES" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.temp.bam" "${sampleId}.temp.bam.bai" "${sampleId}.temp.bam.bai"
    "${sampleId}.gene_counts.summary")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("task.cpus" "saf" "sampleId" "bam" "workflow.manifest.version")

    parameters=("${task.cpus}" "${saf}" "${sampleId}" "${bam}"
    "${workflow.manifest.version}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "TAG_GENES" "\$param_name" "\$parameter"
    done
    """
}
