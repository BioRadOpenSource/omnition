/*
Predicting number of cells loaded based on UMI knee plots
*/

params.options = [:]

process CELL_CALLING {
    tag "${sampleId}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(count_mtx), path(barcodes_tsv), path(genes_tsv)
    val images_pulled

    output:
    tuple val(sampleId), path('*.csv'), emit: results
    tuple val(sampleId), path("${sampleId}_barcode_allowlist.csv"), emit: allowlist
    tuple val(sampleId), path("*{allalgos,barcode_rank,numcell_analysis}*.csv"), emit: all_but_allowlist
    path("${sampleId}_CELL_CALLING_messages.txt"), optional: true, emit: messages

    script:
    // Define the override parameter conditionally
    def overrideParam = params.options.barcode.force.get(sampleId) ?
                        "--override ${params.options.barcode.force.get(sampleId)}" : ""
    """
    rnaCallNumCells.R \
        --count_matrix ${count_mtx} \
        --sample_name ${sampleId} \
        ${overrideParam} \
        --results_folder ./

    coreSubsampleByDistance.py -i ${sampleId}_allalgos_loglog.csv -n 1000 -x 2 -y 3
    coreSubsampleByDistance.py -i ${sampleId}_allalgos_cumfrac.csv -n 1000
    """
}
