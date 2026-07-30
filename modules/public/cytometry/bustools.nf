/*
Generating the count matrix
*/

process BUSTOOLS {
    tag "${sampleId}"
    container "${params.container_image.cytometry_public}:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(bus), path(matrix), path(transcripts), path(t2g)
    val images_pulled

    output:
    tuple val(sampleId), path("merge.barcodes.txt"),
        path("merge.genes.txt"), path("merge.mtx"), emit: bustools_counts
    tuple val(sampleId), path("gene.barcodes.txt"), emit: beads
    path 'counter_type', emit: counter_type

    script:
    """
    bustools sort -o output_sort.bus ${bus}
    bustools count -o gene \
        -g ${t2g} \
        -e ${matrix} \
        -t ${transcripts} \
        --genecounts \
        output_sort.bus
    echo 'bustools' > counter_type

    # Modify outputs for core usage
    cp gene.barcodes.txt merge.barcodes.txt
    cp gene.mtx merge.mtx
    mv gene.genes.txt merge.genes.txt
    pigz -k -p ${task.cpus} gene.mtx

    # Get outputs for matrix transposition
    awk '{print \$0 "Null\t"}' merge.genes.txt > ${sampleId}.untransposed.genes.tsv
    cp merge.barcodes.txt ${sampleId}.untransposed.barcodes.tsv
    mv gene.mtx.gz ${sampleId}.untransposed.mtx.gz
    """
}
