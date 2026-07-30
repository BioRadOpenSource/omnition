/*
Tagging BAM with genomic features
*/

process TAG_FEATURES {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    tuple val(sampleId), path(bam), path(index)
    val images_pulled

    output:
    tuple val(sampleId), path('*.featuresTagged.bam'), path('*.featuresTagged.bam.bai'), emit: bam

    script:
    include_introns = params.options.includeIntrons == true ? "--include-introns" : ""
    """
    mkdir -p tmp/
    publicRnaGeneTagger.py -b ${bam} -gt XT -ft XF -t tmp/ -o ./ ${include_introns}
    rm -r tmp
    mv tagged.bam ${sampleId}.featuresTagged.bam

    samtools index -@ ${task.cpus} ${sampleId}.featuresTagged.bam
    """

    stub:
    include_introns = params.options.includeIntrons == true ? "--include-introns" : ""
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=3

    # Check that we have expected number of sample input files
    check_input_files "TAG_FEATURES" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.featuresTagged.bam" "${sampleId}.featuresTagged.bam.bai"
    "${sampleId}.featuresTagged.bam.bai")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("bam" "include_introns" "task.cpus" "workflow.manifest.version" "sampleId")

    parameters=("${bam}" "${include_introns}" "${task.cpus}" "${workflow.manifest.version}"
    "${sampleId}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "TAG_FEATURES" "\$param_name" "\$parameter"
    done
    """
}
