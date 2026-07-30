/*
Counting number of reads per gene annotation per cell
*/

process UMI_TOOLS_COUNT {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(bam), path(bam_index)
    val images_pulled

    output:
    tuple val(sampleId), path('*.counts.tsv.gz'), emit: count

    script:
    include_introns = params.options.includeIntrons == true ? "^Unassigned|INTERGENIC" : "^Unassigned|INTRON|INTERGENIC"
    """
    mkdir -p tmp/

    umi_tools count \
        --method unique \
        --temp-dir "tmp/" \
        --per-cell \
        --per-gene \
        --extract-umi-method tag \
        --cell-tag XC \
        --umi-tag XM \
        --gene-tag XT \
        --skip-tags-regex "${include_introns}" \
        --random-seed 1 \
        -I ${bam} \
        -S ${sampleId}.counts.tsv.gz
    """

    stub:
    include_introns = params.options.includeIntrons == true ? "^Unassigned|INTERGENIC" : "^Unassigned|INTRON|INTERGENIC"
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=2

    # Check that we have expected number of sample input files
    check_input_files "UMI_TOOLS_COUNT" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.counts.tsv.gz")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("workflow.manifest.version" "bam" "include_introns" "sampleId")

    parameters=("${workflow.manifest.version}" "${bam}" "${include_introns}" "${sampleId}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "UMI_TOOLS_COUNT" "\$param_name" "\$parameter"
    done
    """
}
