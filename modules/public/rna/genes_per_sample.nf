/*
Aggregating genic UMI counts per sample into a genes by sample count matrix
*/

process GENES_PER_SAMPLE {
    container "${params.container_image.r_public}:${workflow.manifest.version}"
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
    publicRnaGeneReadCounts.R -o ./ ./
    """

    stub:
    """
    # Create expected output files
    outputs=("gene_umi_counts_per_sample.tsv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("params.options.reportsDir" "workflow.manifest.version")

    parameters=("${params.options.reportsDir}" "${workflow.manifest.version}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "GENES_PER_SAMPLE" "\$param_name" "\$parameter"
    done
    """
}
