/*
Convert allowlist from drop barcodes to bead allowlist
*/

process UNMERGE_ALLOWLIST {
    tag "${sampleId}"
    container "${params.container_image.dbg_public}:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    tuple val(sampleId), path(allowlist), path(translate)
    val images_pulled

    output:
    tuple val(sampleId), path("*_unmerged_allowlist.csv"), emit: allowlist

    script:
    """
    publicRnaUnmergeAllowlist.py -a ${allowlist} -b ${translate}
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=2

    # Check that we have expected number of sample input files
    check_input_files "UNMERGE_ALLOWLIST" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}_barcode_unmerged_allowlist.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("workflow.manifest.version" "sampleId" "translate" "allowlist")

    parameters=("${workflow.manifest.version}" "${sampleId}" "${translate}" "${allowlist}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "UNMERGE_ALLOWLIST" "\$param_name" "\$parameter"
    done
    """
}
