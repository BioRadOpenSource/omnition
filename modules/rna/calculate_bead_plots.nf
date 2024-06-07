/*
Calculate bead plots
*/

params.options = [:]

process CALCULATE_BEAD_PLOTS {
    tag "${sampleId}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(bead_summary)
    val images_pulled

    output:
    tuple val(sampleId), path("*_partition_poisson.csv"), emit: distributions
    tuple val(sampleId), path("*_raw_poisson.csv"), emit: raw_distributions

    script:
    """
    # process above knee beads from summary to create per-drop distributions for report
    rnaCalcBeadPlots.R \
                    -p ${sampleId} \
                    -bs ${sampleId}.bead_summary.csv
    """
}
