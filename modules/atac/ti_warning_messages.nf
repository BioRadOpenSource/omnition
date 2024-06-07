/*
Gather TI warning messages into warning messages file
*/

process TI_WARNING_MESSAGES {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    errorStrategy 'terminate'
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    path sample_messages
    path messages

    output:
    path 'messages.txt', includeInputs: true, emit: messages

    script:
    """
    cat ${sample_messages} >> ${messages}
    """
}
