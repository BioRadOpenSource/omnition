/*
Decompressing files
*/

process GUNZIP {
    tag "${file}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xxsmall'

    input:
    path file
    val images_pulled

    output:
    path "${out}", emit: decompressed

    script:
    if (file.name =~ /gz$/) {
        out = file.baseName
        """
        gunzip -f ${file}
        """
    } else {
        out = file
        """
        touch ${file}
        """
    }

    stub:
    out = file
    """
    file=\${file/.gz/}
    # Create expected output files
    outputs=("$file")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("file" "workflow.manifest.version" "out")

    parameters=("${file}" "${workflow.manifest.version}" "${out}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "GUNZIP" "\$param_name" "\$parameter"
    done
    # Check that we have expected number of sample input files
    check_input_files "GUNZIP" "\$param_name" "\$EXPECTED_INPUT_COUNT"
    """
}
