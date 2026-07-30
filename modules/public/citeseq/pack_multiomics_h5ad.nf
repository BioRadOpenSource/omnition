/*
Combine RNA and ADT h5ad files into a single MuData object
*/

process PACK_MULTIOMICS_H5AD {
    container "${params.container_image.dbg_public}:${workflow.manifest.version}"
    publishDir "${params.options.reportsDir}/", mode: 'copy', overwrite: true

    label 'cpu_medium'
    label 'memory_medium'

    input:
    path rna_h5ad
    path adt_h5ad
    val images_pulled

    output:
    path 'multiomics.h5mu', emit: h5mu

    script:
    """
    publicCiteseqMultiomicsH5ad.py --rna ${rna_h5ad} --adt ${adt_h5ad} --output multiomics.h5mu
    """

    stub:
    """
    # Create expected output files
    outputs=("multiomics.h5mu")

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
