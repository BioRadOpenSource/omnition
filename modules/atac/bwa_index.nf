/*
Indexing the reference genome for read mapping with BWA
*/

params.options = [:]

process BWA_INDEX {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    publishDir "${params.options.reference.directory}/", mode: 'copy', overwrite: true
    label 'cpu_small'
    label 'memory_large'

    input:
    path fasta
    val images_pulled

    output:
    path 'bwa-index/', emit:index

    script:
    """
    mkdir -p bwa-index/
    bwa index -p bwa-index/${fasta} ${fasta}
    """
}
