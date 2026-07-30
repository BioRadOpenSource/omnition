/*
Validate FASTQ files
*/

process VALIDATE_FASTQS {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    errorStrategy "terminate"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(forward), path(reverse)
    val images_pulled

    output:
    path("*_errors.txt"), emit: error_files
    val true, emit:passed

    script:
    // Validate fastqs and capture any errors
    """
    # Reducing heap size to 80% of allocated resources account for Java overhead
    JAVA_MEMORY="\$(((${task.memory.toGiga()} * 4)/ 5))g"

    touch ${sampleId}_errors.txt
    reformat.sh \
        in=${forward} \
        in2=${reverse} \
        vpair \
        ain=t \
        trd=t \
        t=${task.cpus} \
        -Xmx\$JAVA_MEMORY \
        >> out.txt 2>&1 || \
        # Capture the error message with sample id
        printf "${sampleId}\n\$(tail -n+8 out.txt | head -n 2)\n\n" \
        >> ${sampleId}_errors.txt
    # is the file properly compressed?
    if [[ ${forward} == *".gz" ]];
    then
        gzip -l ${forward} || printf "${forward} is not properly compressed\n\n" \
            >> ${sampleId}_errors.txt
        gzip -l ${reverse} || printf "${reverse} is not properly compressed\n\n" \
            >> ${sampleId}_errors.txt
    fi
    """

    stub:
    // Validate fastqs and capture any errors
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=2

    # Check that we have expected number of sample input files
    check_input_files "VALIDATE_FASTQS" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}_errors.txt")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("task.cpus" "forward" "sampleId" "workflow.manifest.version" "reverse"
    "task.memory.toGiga()")

    parameters=("${task.cpus}" "${forward}" "${sampleId}" "${workflow.manifest.version}"
    "${reverse}" "${task.memory.toGiga()}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "VALIDATE_FASTQS" "\$param_name" "\$parameter"
    done
    """
}
