/*
Validate the input TI config
*/

params.options = [:]

process VALIDATE_TI_CONFIG {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    errorStrategy 'terminate'
    label 'cpu_xsmall'
    label 'memory_xxsmall'

    input:
    path config
    val images_pulled
    path asset

    output:
    path("*.validated.csv"), emit: config
    path("TIlen.txt"), emit: tilen

    script:
    """
    echo "${params.options.sampleIds}" >> samples.txt
    atacValidateTiConfig.py --config $config --sample-list samples.txt --index-list \
    $asset --override ${params.options.tiErrorOverride}
    """
}
