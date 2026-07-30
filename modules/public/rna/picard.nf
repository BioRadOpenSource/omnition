/*
Calculating various alignment metrics
*/

process PICARD {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(bam), path(refflat), path(interval_list)
    val images_pulled

    output:
    tuple val(sampleId), path('*.rna_seq_metrics.txt'), emit: metrics

    script:
    // Define the strand specificity parameter
    def strandSpecificity = params.options.strandSpecificity == "first" ?
                        "FIRST_READ_TRANSCRIPTION_STRAND" : "SECOND_READ_TRANSCRIPTION_STRAND"
    """
    # Reducing heap size to 80% of allocated resources account for Java overhead
    JAVA_MEMORY="\$(((${task.memory.toGiga()} * 4)/ 5))g"

    picard CollectRnaSeqMetrics \
        INPUT=${bam} \
        OUTPUT=${sampleId}.rna_seq_metrics.txt \
        REF_FLAT=${refflat} \
        STRAND_SPECIFICITY=${strandSpecificity} \
        RIBOSOMAL_INTERVALS=${interval_list} \
        -Xmx\$JAVA_MEMORY
    """

    stub:
    """
    # Check if mixed_species is true
    if [ "${params.options.mixed}" == "true" ]; then
        # Expected number of sample input files for mixed species
        EXPECTED_INPUT_COUNT=0
    else
        # Expected number of sample input files for single species
        EXPECTED_INPUT_COUNT=1
    fi
    outputs=("${sampleId}.rna_seq_metrics.txt")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("task.memory.toGiga()" "bam" "refflat" "workflow.manifest.version"
    "interval_list" "sampleId" "strandSpecificity")

    parameters=("${task.memory.toGiga()}" "${bam}" "${refflat}"
    "${workflow.manifest.version}" "${interval_list}" "${sampleId}" "${strandSpecificity}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "PICARD" "\$param_name" "\$parameter"
    done
    """
}
