/*
Normalizing and clustering data with seurat
*/

process SEURAT {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/seurat/",
        pattern: '*.rds', mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_xlarge'
    beforeScript 'ulimit -u $(ulimit -Hu)'

    input:
    tuple val(sampleId), path(matrix), path(barcodes), path(features)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}_seurat.rds"), emit: rds
    tuple val(sampleId), path("${sampleId}_umap.csv"), emit: umap_csv
    tuple val(sampleId), path("${sampleId}_top_features.csv"), emit: top_features_csv
    tuple val(sampleId), path("${sampleId}_seurat_metadata.csv.gz"), emit: metadata
    path("${sampleId}_SEURAT_messages.txt"), optional: true, emit: messages

    script:
    """
    publicRnaSeurat.R ${sampleId} ./ ${params.options.mixed}
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=7

    # Check that we have expected number of sample input files
    check_input_files "SEURAT" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}_seurat.rds" "${sampleId}_umap.csv" "${sampleId}_top_features.csv"
    "${sampleId}_seurat_metadata.csv.gz")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("workflow.manifest.version" "params.options.mixed"
    "params.options.resultsDir" "sampleId")

    parameters=("${workflow.manifest.version}" "${params.options.mixed}"
    "${params.options.resultsDir}" "${sampleId}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "SEURAT" "\$param_name" "\$parameter"
    done
    """
}
