/*
Combining mixed species references together
*/

process COMBINE_REFERENCES {
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    publishDir "${params.options.reference.directory}/", mode: 'copy', overwrite: true
    label 'cpu_small'
    label 'memory_small'

    input:
    path fasta
    path gtf
    val images_pulled

    output:
    path 'mixed.fa', emit: fasta
    path 'mixed.gtf', emit: gtf

    script:
    """
    # Append species names onto contigs in fasta
    for REF in \$(ls ${fasta}); do
        # Get species prefix from gtf filename
        SPECIES=\$(echo \$REF | awk -F '[.]' '{print \$1}' | sed 's;-;_;g')

        # Paste species name onto contigs in fasta
        sed "s/>/>\$SPECIES./g" \$REF >> mixed.fa
    done

    # Appending gtf headers to new file
    grep -h "#" \$(ls ${gtf}) > mixed.gtf || true

    # Append species names onto contigs in gtf
    for REF in \$(ls ${gtf}); do
        # Get species prefix from gtf filename
        SPECIES=\$(echo \$REF | awk -F '[.]' '{print \$1}' | sed 's;-;_;g')

        # Paste species name onto contigs in gtf
        grep -v "#" \$REF | sed "s/^/\$SPECIES./" >> mixed.gtf
    done
    """

    stub:
    """
    # Create expected output files
    outputs=("mixed.fa" "mixed.gtf")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("workflow.manifest.version" "fasta" "gtf"
    "params.options.reference.directory")

    parameters=("${workflow.manifest.version}" "${fasta}" "${gtf}"
    "${params.options.reference.directory}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "COMBINE_REFERENCES" "\$param_name" "\$parameter"
    done
    # Check that we have expected number of sample input files
    check_input_files "COMBINE_REFERENCES" "\$param_name" "\$EXPECTED_INPUT_COUNT"
    """
}
