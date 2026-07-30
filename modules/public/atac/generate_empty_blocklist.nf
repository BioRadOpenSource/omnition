/*
Create an empty blocklist
*/

process GENERATE_EMPTY_BLOCKLIST {
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xxsmall'

    input:
        val images_pulled

    output:
        path 'blocklist.bed', emit:bed

    script:
    '''
    touch blocklist.bed
    '''
}
