/*
Combine read counts files
*/

process COMBINE_READ_COUNTS {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(counts_duplicated), path(counts_unmapped, stageAs: 'unmapped.sequence_counts.csv')
    val images_pulled

    output:
    tuple val(sampleId), path('*.counts_per_barcode.csv'), emit: sequence_counts

    script:
    """
    # Generate a per-bead per-contig count table
    publicRnaCombineReadCounts.R ${counts_duplicated} unmapped.sequence_counts.csv -p ${sampleId}
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=1

    # Check that we have expected number of sample input files
    check_input_files "COMBINE_READ_COUNTS" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.counts_per_barcode.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("params.options.resultsDir" "workflow.manifest.version" "sampleId"
    "counts_duplicated")

    parameters=("${params.options.resultsDir}" "${workflow.manifest.version}" "${sampleId}"
    "${counts_duplicated}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "COMBINE_READ_COUNTS" "\$param_name" "\$parameter"
    done
    """
}
