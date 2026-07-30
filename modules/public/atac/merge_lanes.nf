/*
Merge lane split fastqs
*/

process MERGE_LANES {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(forward, stageAs: 'forward/*'), path(reverse, stageAs: 'reverse/*')
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}_*_001.complete.fastq*"), emit: complete

    script:
    """
    FILECOUNT=\$(ls -Rp | grep -v / | grep -v ".:" | wc -w)
    GZCOUNT=\$(ls */ | grep ".gz" |wc -w)

    if [ \$GZCOUNT = 0 ];
    then
        if [ \$FILECOUNT = 2 ];
        then
            mv forward/* ${sampleId}_R1_001.complete.fastq
            mv reverse/* ${sampleId}_R2_001.complete.fastq
        else
            cat forward/* > ${sampleId}_R1_001.complete.fastq
            cat reverse/* > ${sampleId}_R2_001.complete.fastq
        fi
    else
        if [ \$GZCOUNT != \$FILECOUNT ];
        then
            echo "ERROR: All fastqs associated with ${sampleId} must be either compressed or uncompressed."
            exit 1
        elif [ \$FILECOUNT = 2 ];
        then
            mv forward/* ${sampleId}_R1_001.complete.fastq.gz
            mv reverse/* ${sampleId}_R2_001.complete.fastq.gz
        else
            cat forward/* > ${sampleId}_R1_001.complete.fastq.gz
            cat reverse/* > ${sampleId}_R2_001.complete.fastq.gz
        fi
    fi
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=2

    # Check that we have expected number of sample input files
    check_input_files "MERGE_LANES" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=(
        "${sampleId}_R1_001.complete.fastq.gz"
        "${sampleId}_R2_001.complete.fastq.gz"
    )

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("params.options.input" "sampleId" "workflow.manifest.version")

    parameters=("${params.options.input}" "${sampleId}" "${workflow.manifest.version}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "MERGE_LANES" "\$param_name" "\$parameter"
    done
    """
}
