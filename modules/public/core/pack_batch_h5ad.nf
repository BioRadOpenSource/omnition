/*
Concatenate all h5ad files, dim reduction, clustering
*/

process PACK_BATCH_H5AD {
    container "${params.container_image.dbg_public}:${workflow.manifest.version}"
    publishDir "${params.options.reportsDir}/", mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_medium'

    input:
    path h5ad
    val images_pulled

    output:
    path 'all_samples*.h5ad', emit: h5ad

    script:
    """
    publicCoreMergeH5ad.py -i ./
    """

    stub:
    """
    # Create expected output files
    outputs=("all_samples.h5ad")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("workflow.manifest.version" "params.options.reportsDir")

    parameters=("${workflow.manifest.version}" "${params.options.reportsDir}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "PACK_BATCH_H5AD" "\$param_name" "\$parameter"
    done
    """
}
