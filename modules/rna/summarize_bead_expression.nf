/*
Count the frequency of different features
*/

params.options = [:]

process SUMMARIZE_BEAD_EXPRESSION {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_xlarge'

    input:
    tuple val(sampleId), path(bam), path(bam_index), path(interval_list)
    val images_pulled

    output:
    tuple val(sampleId), path("*.bead_scrnaseq_counts.csv"), emit: count

    script:
    """
    # Initialize optional arguments
    optional_args=""

    # Check if mixed species option is set
    if ${params.options.mixed}; then
        header_info=\$(samtools view -H ${bam} | grep @SQ | cut -f2 | cut -f2 -d:)
        # get species reference prefixes off bam header
        species1=\$(echo "\$header_info" | cut -f1 -d. | sort -T '.' | uniq | tail -n 2 | head -n 1)
        species2=\$(echo "\$header_info" | cut -f1 -d. | sort -T '.' | uniq | tail -n 2 | tail -n 1)

        # get mitochondrial contigs as an array
        mito_contigs=( \$(echo "\$header_info" | grep 'chrM\\|MT' | true) )

        # check if any mitochondrial contigs were found and create a string from the array
        if (( \${#mito_contigs[@]} )); then
            mito_arg=\$(printf -- "--mito-contig %s " "\${mito_contigs[@]}")
        else
            mito_arg="--mito-contig MT"
        fi

        # add additional args for counting species-specific umis
        optional_args="--species-mix True --species-id-1 \${species1} --species-id-2 \${species2} \$mito_arg"
    else
        optional_args="--mito-contig MT"
    fi

    rnaCountSingleCell.py -i ${bam} --barcode-tag XB \
    --umi-tag XM --gene-tag XT \
    --feature-tag XF -c ${task.cpus} \
    --ribosomal-interval ${interval_list} \
    -o ./ \$optional_args

    mv \$(ls *.scrnaseq_counts.csv) ${sampleId}.bead_scrnaseq_counts.csv
    """
}
