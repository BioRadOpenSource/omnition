/*
Calculate bead plots
*/

process CALCULATE_BEAD_PLOTS {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(bead_summary)
    val images_pulled

    output:
    tuple val(sampleId), path("*_partition_poisson.csv"), emit: distributions
    tuple val(sampleId), path("*_raw_poisson.csv"), emit: raw_distributions
    tuple val(sampleId), path("*.gel_bead_lambda_summary.csv"), emit: lambda_summary

    script:
    """
    # process above knee beads from summary to create per-drop distributions for report
    publicCoreCalcBeadPlots.R \
                    -p ${sampleId} \
                    -bs ${sampleId}.bead_summary.csv
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=6

    # Check that we have expected number of sample input files
    check_input_files "CALCULATE_BEAD_PLOTS" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.above_knee_beads_per_partition_poisson.csv"
    "${sampleId}.beads_per_partition_poisson.csv"
    "${sampleId}.above_knee_beads_per_partition_raw_poisson.csv"
    "${sampleId}.gel_bead_lambda_summary.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("sampleId"
    "workflow.manifest.version" "params.options.resultsDir"
    "params.options.estimatedDroplets")

    parameters=("${sampleId}"
    "${workflow.manifest.version}" "${params.options.resultsDir}"
    "${params.options.estimatedDroplets}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "CALCULATE_BEAD_PLOTS" "\$param_name" "\$parameter"
    done
    """
}
