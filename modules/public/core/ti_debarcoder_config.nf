/*
Create a json config file of TIs for debarcoder
*/

process TI_DEBARCODER_CONFIG {
    container "${params.container_image.dbg_public}:${workflow.manifest.version}"
    errorStrategy 'terminate'
    label 'cpu_xsmall'
    label 'memory_xxsmall'

    input:
    tuple val(sampleId), path(asset), path(fastq_valid)
    val images_pulled

    output:
    tuple val(sampleId), path("TI.json"), emit: config, optional: true

    script:
    """
    publicCoreConfigToDebarcoderJson.py \
        -i ${asset} \
        -m 1 \
        -r ${params.options.tiRead} \
        -t ${params.options.i7AsTi} \
        -c /opt/biorad/assets/dead/atac.json
    """
}
