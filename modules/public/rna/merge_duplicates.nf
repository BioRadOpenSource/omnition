/*
Convert duplicate count file from bead barcodes to drop barcodes
*/

process MERGE_DUPLICATES {
    tag "${sampleId}"
    container "${params.container_image.dbg_public}:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    tuple val(sampleId), path(duplicates), path(translate)
    val images_pulled

    output:
    tuple val(sampleId), path("*.merged_duplicate_counts.csv"), emit: duplicate_count

    script:
    """
    publicRnaMergeDuplicates.py -d ${duplicates} -b ${translate}
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=2

    # Check that we have expected number of sample input files
    check_input_files "MERGE_DUPLICATES" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.merged_duplicate_counts.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("translate" "workflow.manifest.version" "sampleId" "duplicates")

    parameters=("${translate}" "${workflow.manifest.version}" "${sampleId}" "${duplicates}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "MERGE_DUPLICATES" "\$param_name" "\$parameter"
    done
    """
}
