/*
Create a json config file of TIs for debarcoder
*/

params.options = [:]

process TI_DEBARCODER_CONFIG {
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    errorStrategy 'terminate'
    label 'cpu_xsmall'
    label 'memory_xxsmall'

    input:
    tuple val(sampleId), path(asset), path(fastq_valid)
    val images_pulled

    output:
    tuple val(sampleId), path("TI.json"), emit: config

    script:
    """
    atacConfigToDebarcoderJson.py \
        -i ${asset} \
        -m 1 \
        -r ${params.options.tiRead} \
        -t ${params.options.i7AsTi} \
        -c /opt/biorad/assets/dead/atac.json
    """
}
