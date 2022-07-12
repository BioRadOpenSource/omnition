/*
Formatting single species blocklists
*/

params.options = [:]

process FORMAT_BLOCKLIST {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    if (workflow.profile == 'aws') {
        label 'small'
    } else {
        label 'cpu_xsmall'
        label 'memory_xxsmall'
    }

    input:
    path blocklist
    val images_pulled

    output:
    path 'formatted.*.bed', emit: bed

    script:
    """
    # Get species prefix from blocklist filename
    SPECIES=\$(echo ${blocklist} | awk -F '[.]' '{print \$1}' | sed 's;-;_;g')

    # Paste species name onto contigs in blocklist
    sed "s/^/\$SPECIES./" ${blocklist} >> formatted.${blocklist}
    """
}
