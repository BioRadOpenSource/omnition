/*
Predicting number of cells loaded based on UMI knee plots
*/

process CELL_CALLING {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(count_mtx), path(barcodes_tsv), path(genes_tsv)
    val images_pulled

    output:
    tuple val(sampleId), path('*.csv'), emit: results
    tuple val(sampleId), path("${sampleId}_barcode_allowlist.csv"), emit: allowlist
    tuple val(sampleId), path("*numcell_analysis.csv"), emit: numcells
    tuple val(sampleId), path("*{allalgos,barcode_rank,numcell_analysis}*.csv"), emit: all_but_allowlist
    path("${sampleId}_CELL_CALLING_messages.txt"), optional: true, emit: messages

    script:
    // Define the override parameter conditionally
    def overrideParam = params.options.barcode.force.get(sampleId) ?
                        "--override ${params.options.barcode.force.get(sampleId)}" : ""
    """
    publicRnaCallNumCells.R \
        --count_matrix ${count_mtx} \
        --sample_name ${sampleId} \
        ${overrideParam} \
        --results_folder ./

    publicCoreSubsampleByDistance.py -i ${sampleId}_allalgos_loglog.csv -n 1000 -x 2 -y 3
    publicCoreSubsampleByDistance.py -i ${sampleId}_allalgos_cumfrac.csv -n 1000
    """

    stub:
    // Define the override parameter conditionally
    def overrideParamStub = params.options.barcode.force.get(sampleId) ?
                        "--override ${params.options.barcode.force.get(sampleId)}" : ""
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=10

    # Check that we have expected number of sample input files
    check_input_files "CELL_CALLING" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}_barcode_allowlist.csv" "${sampleId}_numcell_analysis.csv"
    "${sampleId}_allalgos_cumfrac.csv" "${sampleId}_allalgos_loglog.csv"
    "${sampleId}_barcode_rank.csv" "${sampleId}_allalgos_loglog_downsampled.csv"
    "${sampleId}_allalgos_cumfrac_downsampled.csv" "${sampleId}_barcode_allowlist.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("params.options.barcode.force.get(sampleId)"
    "params.options.cell.loaded.get(sampleId)" "sampleId" "workflow.manifest.version"
    "overrideParamStub" "count_mtx")

    parameters=("${params.options.barcode.force.get(sampleId)}"
    "${params.options.cell.loaded.get(sampleId)}" "${sampleId}"
    "${workflow.manifest.version}" "${overrideParamStub}" "${count_mtx}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "CELL_CALLING" "\$param_name" "\$parameter"
    done
    """
}
