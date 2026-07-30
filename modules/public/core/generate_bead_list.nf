/*
Generate a list of beads to be merged from the debarcoded FASTQ.
*/

process GENERATE_BEAD_LIST {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_small'

    input:
    tuple val(sampleId), path(fastq)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}_all_beads.tsv"), emit: beads

    script:
    """
    publicDoMergingGenerateBeadList.py \
        -f ${fastq[1]} \
        -o ${sampleId}_all_beads.tsv \
        -s ${sampleId}
    """
    stub:
    """
    # Create expected output files
    outputs=("${sampleId}_all_beads.tsv")
    for output in "\${outputs[@]}"
    do
        touch \$output
    done
    # Record all groovy parameters used in module
    param_names=("fastq" "workflow.manifest.version" "sampleId")

    parameters=("${fastq}" "${workflow.manifest.version}" "${sampleId}")
    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "GENERATE_BEAD_LIST" "\$param_name" "\$parameter"
    done
    """
}
