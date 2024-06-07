/*
Combining mixed species references together
*/

params.options = [:]

process COMBINE_REFERENCES {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
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
}
