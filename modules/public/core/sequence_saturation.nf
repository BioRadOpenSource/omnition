/*
Calculate sequence saturation for each sample
*/

process SEQUENCE_SATURATION {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(fragments), path(summary)
    val assay
    val images_pulled

    output:
    tuple val(sampleId), path("*_sequence_saturation.csv"), emit: results

    script:
    cDNA_counts = params.options.total_cDNA_reads?.get(sampleId) ?: false
    """
    publicCoreSequenceSaturation.R \
        --sample_id ${sampleId} \
        --feature_counts_file ${fragments} \
        --metrics_file ${summary} \
        --cDNA_counts ${cDNA_counts} \
        --assay "${assay}"
    """

    stub:
    cDNA_counts = params.options.total_cDNA_reads?.get(sampleId) ?: false
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=1

    # Check that we have expected number of sample input files
    check_input_files "SEQUENCE_SATURATION" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}_sequence_saturation.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("summary" "cDNA_counts" "assay" "sampleId" "workflow.manifest.version"
    "fragments")

    parameters=("${summary}" "${cDNA_counts}" "${assay}" "${sampleId}"
    "${workflow.manifest.version}" "${fragments}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "SEQUENCE_SATURATION" "\$param_name" "\$parameter"
    done
    """
}
