/*
Mapping reads to antibody tags
*/

process KALLISTO {
    tag "${sampleId}"
    container "${params.container_image.cytometry_public}:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(cor_read_pair), path(index)
    val images_pulled

    output:
    tuple val(sampleId), path('output.bus'), path('matrix.ec'), path('transcripts.txt'), emit: kallisto
    tuple val(sampleId), path('*_run_info.json'), emit: kallisto_report

    script:
    """
    CBCSTART=0
    CBCEND=21
    UMISTART=21
    UMIEND=29

    kallisto bus -i ${index} \
        -o \$PWD \
        -x "0,\$CBCSTART,\$CBCEND:0,\$UMISTART,\$UMIEND:1,0,0" \
        -t \$(nproc) \
        ${cor_read_pair[0]} ${cor_read_pair[1]}

    # Rename run_info.json to include sample ID
    mv run_info.json ${sampleId}_run_info.json
    """
}
