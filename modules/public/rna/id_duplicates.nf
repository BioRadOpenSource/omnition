/*
Mark duplicate reads
*/

process ID_DUPLICATES {
    tag "${sampleId}"
    container "${params.container_image.dbg_public}:${workflow.manifest.version}"
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

    stub:
    include_introns = params.options.includeIntrons == true ? "--include-introns" : ""
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=3

    # Check that we have expected number of sample input files
    check_input_files "ID_DUPLICATES" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.duplicate_counts.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("include_introns" "workflow.manifest.version" "allowlist" "sampleId" "bam")

    parameters=("${include_introns}" "${workflow.manifest.version}" "${allowlist}"
    "${sampleId}" "${bam}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "ID_DUPLICATES" "\$param_name" "\$parameter"
    done
    """
}
