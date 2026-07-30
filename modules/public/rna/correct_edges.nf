/*
Correcting DO and Bead UMI sequences
*/

process CORRECT_EDGES {
    tag "${sampleId}"
    container "${params.container_image.dbg_public}:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_large'

    input:
    tuple val(sampleId), path(edges)
    val images_pulled

    output:
    tuple val(sampleId), path("*_corrected_edges.tsv"), emit: edges
    tuple val(sampleId), path("*_corrected_edge_counts.csv"), emit: counts

    script:
    """
    omnition core correct-edge-umi -i ${edges} \
        -s ${sampleId} \
        -t ${params.options.umiHamming.get( sampleId )} \
        -c ${task.cpus} \
        -o ./
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=3

    # Check that we have expected number of sample input files
    check_input_files "CORRECT_EDGES" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}_corrected_edges.tsv" "${sampleId}_corrected_edge_counts.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("params.options.umiHamming.get( sampleId )" "workflow.manifest.version"
    "params.options.resultsDir" "edges" "task.cpus" "sampleId")

    parameters=("${params.options.umiHamming.get( sampleId )}"
    "${workflow.manifest.version}" "${params.options.resultsDir}" "${edges}" "${task.cpus}"
    "${sampleId}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "CORRECT_EDGES" "\$param_name" "\$parameter"
    done
    """
}
