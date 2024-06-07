/*
Aggregating genic UMI counts per sample into a genes by sample count matrix
*/

params.options = [:]

process GENES_PER_SAMPLE {
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    publishDir "${params.options.reportsDir}/", mode:'copy', overwrite: true
    label 'cpu_small'
    label 'memory_small'

    input:
    path count_matrix
    val images_pulled

    output:
    path '*_umi_counts_per_sample.tsv', emit: count

    script:
    """
    rnaGeneReadCounts.R -o ./ ./
    """
}
