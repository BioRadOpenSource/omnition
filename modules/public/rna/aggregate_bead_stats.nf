/*
Aggregate bead statistics
*/

process AGGREGATE_BEAD_STATS {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(sequence_counts), path(allowlist), path(cell_expression), path(bead_expression),
        path(barcode_translate)
    val images_pulled

    output:
    tuple val(sampleId), path('*.bead_summary.csv'), emit: bead_summary
    tuple val(sampleId), path('*.bead_merge_sum_exp_metrics.csv'), emit: metrics

    script:
    """
    publicRnaBeadStats.R -bc ${sequence_counts} \
        -a ${allowlist} \
        -c ${cell_expression} \
        -b ${bead_expression} \
        -bt ${barcode_translate} \
        -s ${sampleId}
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=5

    # Check that we have expected number of sample input files
    check_input_files "AGGREGATE_BEAD_STATS" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.bead_summary.csv" "${sampleId}.bead_merge_sum_exp_metrics.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("sequence_counts" "workflow.manifest.version" "params.options.resultsDir"
    "cell_expression" "allowlist" "bead_expression" "sampleId" "barcode_translate")

    parameters=("${sequence_counts}" "${workflow.manifest.version}"
    "${params.options.resultsDir}" "${cell_expression}" "${allowlist}" "${bead_expression}"
    "${sampleId}" "${barcode_translate}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "AGGREGATE_BEAD_STATS" "\$param_name" "\$parameter"
    done
    """
}
