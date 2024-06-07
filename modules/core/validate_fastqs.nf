/*
Validate FASTQ files
*/

params.options = [:]

process VALIDATE_FASTQS {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
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
}
