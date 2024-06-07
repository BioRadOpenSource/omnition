/*
Indexing references to prepare for alignment
*/

params.options = [:]

process STAR_INDEX {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    publishDir "${params.options.reference.directory}/", mode: 'copy', overwrite: true
    label 'cpu_large'
    label 'memory_xlarge'

    input:
    path fasta
    path gtf
    val images_pulled

    output:
    path 'star-index/', emit: index

    script:
    """
    mkdir -p star-index/

    STAR --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir star-index/ \
        --genomeFastaFiles ${fasta} \
        --sjdbGTFfile ${gtf} \
        --genomeSAsparseD 3 \
        --sjdbOverhang 100
    """
}
