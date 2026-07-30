/*
Count the frequency of different features
*/

process SUMMARIZE_EXPRESSION {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_xlarge'

    input:
    tuple val(sampleId), path(bam), path(bam_index), path(interval_list)
    val images_pulled

    output:
    tuple val(sampleId), path("*.scrnaseq_counts.csv"), emit: count

    script:
    """
    # Check if mixed species option is set
    if ${params.options.mixed}; then
        # get species reference prefixes off bam header
        species1=\$(samtools view -H ${bam} | grep @SQ | cut -f2 | cut -f2 -d: | \
            cut -f1 -d. | sort -T '.' | uniq | tail -n 2 | head -n 1)
        species2=\$(samtools view -H ${bam} | grep @SQ | cut -f2 | cut -f2 -d: | \
            cut -f1 -d. | sort -T '.' | uniq | tail -n 2 | tail -n 1)

        # get mitochondrial contigs as an array
        mito_contigs=( \$(samtools view -H ${bam} | grep @SQ | cut -f2 | cut -f2 -d: | grep 'chrM\\|MT' | true) )

        # check if any mitochondrial contigs were found and create a string from the array
        if (( \${#mito_contigs[@]} )); then
            mito_arg=\$(printf -- "--mito-contig %s " "\${mito_contigs[@]}")
        else
            mito_arg="--mito-contig MT"
        fi

        # read counter with additional args for counting species-specific umis
        # mito_arg is passed as an entire string
        publicRnaCountSingleCell.py -i ${bam} \
            --barcode-tag XC \
            --umi-tag XM \
            --gene-tag XT \
            --feature-tag XF \
            --species-mix True \
            --species-id-1 \${species1} \
            --species-id-2 \${species2} \
            \$mito_arg \
            -c ${task.cpus} \
            --ribosomal-interval ${interval_list} \
            -o ./
    else
        publicRnaCountSingleCell.py -i ${bam} \
            --barcode-tag XC \
            --umi-tag XM \
            --gene-tag XT \
            --feature-tag XF \
            -c ${task.cpus} \
            --mito-contig MT \
            --ribosomal-interval ${interval_list} \
            -o ./
    fi
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=2

    # Check that we have expected number of sample input files
    check_input_files "SUMMARIZE_EXPRESSION" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.final.scrnaseq_counts.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("params.options.mixed" "task.cpus" "interval_list" "sampleId"
    "workflow.manifest.version" "params.options.resultsDir" "bam")

    parameters=("${params.options.mixed}" "${task.cpus}" "${interval_list}" "${sampleId}"
    "${workflow.manifest.version}" "${params.options.resultsDir}" "${bam}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "SUMMARIZE_EXPRESSION" "\$param_name" "\$parameter"
    done
    """
}
