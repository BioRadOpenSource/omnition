/*
Combining mixed species blocklists together
*/

params.options = [:]

process COMBINE_BLOCKLISTS {
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
    path 'mixed.blocklist.bed', emit: bed

    script:
    """
    # Append species names onto contigs in blocklist
    for BLIST in ${blocklist}; do

        # Get species prefix from blocklist filename
        SPECIES=\$(echo \$BLIST | awk -F '[.]' '{print \$1}' | sed 's;-;_;g')

        # Paste species name onto contigs in blocklist
        sed "s/^/\$SPECIES./" \$BLIST >> mixed.blocklist.bed

    done
    """
}
