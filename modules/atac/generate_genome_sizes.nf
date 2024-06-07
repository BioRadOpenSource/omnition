/*
Get contig size
*/

params.options = [:]

process GENERATE_GENOME_SIZES {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    publishDir "${params.options.reference.directory}", mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_xxsmall'

    input:
    path fasta
    val images_pulled

    output:
    path 'genome.sizes', emit: size

    script:
    """
    samtools faidx ${fasta}

    cut -f1,2 *.fa.fai > genome.sizes
    """
}
