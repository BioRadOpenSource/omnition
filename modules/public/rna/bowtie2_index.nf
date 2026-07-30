/*
Indexing contaminant genomes to prepare for FastQ Screen
*/

process BOWTIE2_INDEX {
    tag "${fasta}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    path fasta
    val images_pulled

    output:
    path "${fasta.simpleName}-bowtie2-index/", emit: index

    script:
    """
    mkdir ${fasta.simpleName}-bowtie2-index/
    bowtie2-build --threads ${task.cpus} ${fasta} ${fasta.simpleName}-bowtie2-index/${fasta.simpleName}
    """

    stub:
    """
    mkdir "${fasta.simpleName}-bowtie2-index"

    # Create expected output files
    outputs=(
        "${fasta.simpleName}-bowtie2-index/${fasta.simpleName}.1.bt2"
        "${fasta.simpleName}-bowtie2-index/${fasta.simpleName}.2.bt2"
        "${fasta.simpleName}-bowtie2-index/${fasta.simpleName}.3.bt2"
        "${fasta.simpleName}-bowtie2-index/${fasta.simpleName}.4.bt2"
        "${fasta.simpleName}-bowtie2-index/${fasta.simpleName}.rev.1.bt2"
        "${fasta.simpleName}-bowtie2-index/${fasta.simpleName}.rev.2.bt2"
    )

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("workflow.manifest.version" "task.cpus"
    "params.options.contaminant.directory" "fasta" "fasta.simpleName")

    parameters=("${workflow.manifest.version}" "${task.cpus}"
    "${params.options.contaminant.directory}" "${fasta}" "${fasta.simpleName}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "BOWTIE2_INDEX" "\$param_name" "\$parameter"
    done
    """
}
