/*
Filtering Blocklist entries to only contain chromosomes/contigs in the FASTA file
*/

params.options = [:]

process FILTER_BLOCKLISTS {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    publishDir "${params.options.reference.directory}/", mode: 'copy', overwrite: true
    label 'cpu_xsmall'
    label 'memory_xxsmall'

    input:
    path fasta
    path blocklist
    val images_pulled

    output:
    path "filtered.${blocklist.baseName}.bed", emit: blocklist

    script:
    """
    # Transfering over all commented lines into the new filtered file
    awk '/^#/' ${blocklist} > filtered.${blocklist.baseName}.bed

    # Creates an array from lines in the first file that start with > and puts into it all
    # the values in the first element after the >.
    # Then for only the second file it looks to see if the first element in each line of
    # the second file is in the array previously created and only transfers those into the filtered file
    awk '/>/ && NR==FNR{ array[substr(\$1, 2, length(\$1))]; next } \$1 in array' \
        ${fasta} ${blocklist} >> filtered.${blocklist.baseName}.bed
    """
}
