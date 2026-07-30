/*
Collecting the analysis data and generating figures for the final report
*/

process ANALYSIS {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/seurat/",
        mode: 'copy', overwrite: true, pattern: '*_seurat_metadata.csv.gz'
    label 'cpu_medium'
    label 'memory_medium'

    input:
    path counter
    tuple val(sampleId), path(files), path(rds_file), path(umi_percell), path(inital_cell_counts)
    val images_pulled

    output:
    tuple val(sampleId), path("*"), emit: analysis
    tuple val(sampleId), path("fully_processed.cells.csv"), emit: count_input // groovylint-disable-line
    tuple val(sampleId), path("${sampleId}_seurat_metadata.csv.gz"), emit: metadata

    script:
    // Define the override parameter conditionally for Barcodes
    def overrideParamBarcode = params.options.barcode.force.get(sampleId) ?
                        "--barcode_min ${params.options.barcode.force.get(sampleId)}" : ""
    """
    publicCytometrySccBasicPipeline.R --inputdir \$PWD \
        --outputdir \$PWD \
        ${overrideParamBarcode} \
        --sample ${sampleId}
    """
}
