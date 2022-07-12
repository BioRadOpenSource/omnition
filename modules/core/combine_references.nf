/*
Combining mixed species references together
*/

params.options = [:]

process COMBINE_REFERENCES {
    publishDir "${params.options.reference.directory}/", mode: 'copy', overwrite: true
    if (workflow.profile == 'aws') {
        label 'medium'
    } else {
        label 'cpu_small'
        label 'memory_xxsmall'
    }

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
    for REF in ${fasta}; do

      # Get species prefix from gtf filename
      SPECIES=\$(echo \$REF | awk -F '[.]' '{print \$1}' | sed 's;-;_;g')

      # Paste species name onto contigs in fasta
      sed "s/>/>\$SPECIES./g" \$REF >> mixed.fa

    done

    # Appending gtf headers to new file
    grep -h "#" ${gtf} > mixed.gtf || true

    # Append species names onto contigs in gtf
    for REF in ${gtf}; do

      # Get species prefix from gtf filename
      SPECIES=\$(echo \$REF | awk -F '[.]' '{print \$1}' | sed 's;-;_;g')

      # Paste species name onto contigs in gtf
      grep -v "#" \$REF | sed "s/^/\$SPECIES./" >> mixed.gtf

    done
    """
}
