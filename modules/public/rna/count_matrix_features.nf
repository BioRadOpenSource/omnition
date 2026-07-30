/*
Counting expression matrix features
*/

process COUNT_MATRIX_FEATURES {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_small'

    input:
    tuple val(sampleId), path(matrix), path(barcodes), path(features)
    val images_pulled

    output:
    tuple val(sampleId), path('*_matrix_features.csv'), emit: count

    script:
    """
    publicRnaCountMatrixFeatures.R ${sampleId} ${matrix}
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=3

    # Check that we have expected number of sample input files
    check_input_files "COUNT_MATRIX_FEATURES" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}_matrix_features.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("sampleId" "matrix" "workflow.manifest.version")

    parameters=("${sampleId}" "${matrix}" "${workflow.manifest.version}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "COUNT_MATRIX_FEATURES" "\$param_name" "\$parameter"
    done
    """
}
