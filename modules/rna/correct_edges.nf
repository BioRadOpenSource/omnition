/*
Correcting DO and Bead UMI sequences
*/

params.options = [:]

process CORRECT_EDGES {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
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
    rnaCorrectEdges.py -e ${edges} \
        -s ${sampleId} \
        -u ${params.options.umiHamming.get( sampleId )} \
        -c ${task.cpus} \
        -v
    """
}
