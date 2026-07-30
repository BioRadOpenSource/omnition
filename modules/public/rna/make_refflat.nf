/*
Creating refFlat-formatted reference file
*/

process MAKE_REFFLAT {
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    publishDir "${params.options.reference.directory}/", mode: 'copy', overwrite: true
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    path gtf
    val images_pulled

    output:
    path 'annotation.refflat', emit: refflat

    script:
    """
    # Creating extended format refFlat file
    gtfToGenePred -genePredExt ${gtf} annotation.refflat.tmp

    # Reformatting so the output file looks like the basic refFlat output from
    # gtfToGenePred but adds the gene name as the first column
    paste <(cut -f 12 annotation.refflat.tmp) <(cut -f 1-10 annotation.refflat.tmp) > annotation.refflat
    """

    stub:
    """
    # Create expected output files
    outputs=("annotation.refflat.tmp" "annotation.refflat")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("gtf" "params.options.reference.directory" "workflow.manifest.version")

    parameters=("${gtf}" "${params.options.reference.directory}"
    "${workflow.manifest.version}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "MAKE_REFFLAT" "\$param_name" "\$parameter"
    done
    """
}
