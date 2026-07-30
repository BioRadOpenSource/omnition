/*
Indexing references to prepare for alignment
*/

process STAR_INDEX {
    container "${params.container_image.core_public}:${workflow.manifest.version}"
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

    stub:
    """
    mkdir "star-index"

    # Create expected output files
    outputs=(
        "star-index/chrLength.txt"
        "star-index/chrName.txt"
        "star-index/chrNameLength.txt"
        "star-index/chrStart.txt"
        "star-index/exonGeTrInfo.tab"
        "star-index/exonInfo.tab"
        "star-index/geneInfo.tab"
        "star-index/Genome"
        "star-index/genomeParameters.txt"
        "star-index/SA"
        "star-index/SAindex"
        "star-index/sjdbInfo.txt"
        "star-index/sjdbList.fromGTF.out.tab"
        "star-index/sjdbList.out.tab"
        "star-index/transcriptInfo.tab"
    )

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("fasta" "params.options.reference.directory" "task.cpus"
    "workflow.manifest.version" "gtf")

    parameters=("${fasta}" "${params.options.reference.directory}" "${task.cpus}"
    "${workflow.manifest.version}" "${gtf}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "STAR_INDEX" "\$param_name" "\$parameter"
    done
    """
}
