/*
Collecting metadata on edges
*/

process EDGE_METADATA {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xlarge'

    input:
    tuple val(sampleId), path(raw_edges), path(corrected_edges), path(allowlist)
    val images_pulled

    output:
    tuple val(sampleId), path("*edge_metadata.csv"), emit: metadata

    script:
    """
    publicCoreEdgeMetadata.py \
        -r ${raw_edges} \
        -c ${corrected_edges} \
        -b ${allowlist} \
        -s ${sampleId}
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=4

    # Check that we have expected number of sample input files
    check_input_files "EDGE_METADATA" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}_edge_metadata.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("workflow.manifest.version" "corrected_edges" "params.options.resultsDir"
    "raw_edges" "allowlist" "sampleId")

    parameters=("${workflow.manifest.version}" "${corrected_edges}"
    "${params.options.resultsDir}" "${raw_edges}" "${allowlist}" "${sampleId}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "EDGE_METADATA" "\$param_name" "\$parameter"
    done
    """
}
