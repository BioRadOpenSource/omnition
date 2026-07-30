/*
Creating ribosomal interval list reference file
*/

process MAKE_RIBOSOMAL_INTERVAL_LIST {
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    publishDir "${params.options.reference.directory}/", mode: 'copy', overwrite: true
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    path gtf
    path index
    val images_pulled

    output:
    path 'ribosomal.interval_list', emit: interval_list

    script:
    """
    # Adding chromosome names and lengths as header
    awk 'BEGIN{OFS="\\t";} {print "@SQ", "SN:"\$1,"LN:"\$2}' ${index}/chrNameLength.txt > ribosomal.interval_list

    # Appending ribosomal RNA feature information to file
    awk 'BEGIN{OFS="\\t";} /rRNA/ {print \$1,\$4,\$5,\$7,\$9 " " \$10}' \
        ${gtf} | tr -d ';' | uniq >> ribosomal.interval_list
    """

    stub:
    """
    # Create expected output files
    outputs=("ribosomal.interval_list")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("index" "workflow.manifest.version" "gtf"
    "params.options.reference.directory")

    parameters=("${index}" "${workflow.manifest.version}" "${gtf}"
    "${params.options.reference.directory}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "MAKE_RIBOSOMAL_INTERVAL_LIST" "\$param_name" "\$parameter"
    done
    """
}
