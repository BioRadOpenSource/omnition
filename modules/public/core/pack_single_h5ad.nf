/*
Build per-sample H5AD file
*/

process PACK_SINGLE_H5AD {
    tag "${sampleId}"
    container "${params.container_image.dbg_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/countMatrix/", mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(mtx), path(barcodes), path(features), path(summary), path(embedding)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}*.h5ad"), emit: h5ad

    script:
    def summaryArg = summary.name == 'NO_FILE' ? "" : "-d ${summary}"
    """
    publicCoreBuildH5ad.py -m ${mtx} -b ${barcodes} -f ${features} ${summaryArg} -e ${embedding} -n ${sampleId}
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=5

    # Check that we have expected number of sample input files
    check_input_files "PACK_SINGLE_H5AD" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.h5ad")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("features" "summary" "sampleId" "mtx" "embedding"
    "params.options.resultsDir" "barcodes" "workflow.manifest.version")

    parameters=("${features}" "${summary}" "${sampleId}" "${mtx}" "${embedding}"
    "${params.options.resultsDir}" "${barcodes}" "${workflow.manifest.version}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "PACK_SINGLE_H5AD" "\$param_name" "\$parameter"
    done
    """
}
