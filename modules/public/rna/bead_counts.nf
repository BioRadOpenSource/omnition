/*
Count features per barcode
*/

process BEAD_COUNTS {
    tag "${sampleId}"
    container "${params.container_image.dbg_public}:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_large'

    input:
    tuple val(sampleId), path(bam), path(index)
    val images_pulled

    output:
    tuple val(sampleId), path('*.sequence_counts.csv'), emit: sequence_counts

    script:
    """
    # Generate a per-bead per-contig count table
    publicRnaUmiCounter.py -b ${bam} -bt XC -ut XM -c ${task.cpus}
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=2

    # Check that we have expected number of sample input files
    check_input_files "BEAD_COUNTS" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.featuresTagged.sequence_counts.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("workflow.manifest.version" "bam" "sampleId" "task.cpus")

    parameters=("${workflow.manifest.version}" "${bam}" "${sampleId}" "${task.cpus}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "BEAD_COUNTS" "\$param_name" "\$parameter"
    done
    """
}
