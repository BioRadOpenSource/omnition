/*
Count the frequency of different features for mixed species
*/

params.options = [:]

process MIXED_EXPRESSION {
    tag "${sampleId}-${species}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_small'

    input:
    tuple val(sampleId), val(species), path(allowlist), path(symbols), path(count)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}.${species}_matrix_features.csv"), emit:mixed_features

    script:
    """
    coreBuildCountMatrix.R -o ./ -p ${sampleId}.${species}.filtered -w ${allowlist} -g ${symbols} ${count} -m mixed
    pigz -p ${task.cpus} ${sampleId}.${species}.filtered.mtx

    rnaCountMatrixFeatures.R ${sampleId}.${species} ${sampleId}.${species}.filtered.mtx.gz
    """
}
