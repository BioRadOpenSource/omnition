/*
Combining mixed species blocklists together
*/

params.options = [:]

process COMBINE_BLOCKLISTS {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    errorStrategy 'terminate'
    label 'cpu_xsmall'
    label 'memory_xxsmall'

    input:
    path blocklist
    val images_pulled

    output:
    path 'mixed.blocklist.bed', emit: bed

    script:
    """
    touch reference_errors.txt

    # Append species names onto contigs in blocklist
    for BLIST in ${blocklist}; do
        # Check all second fields are integers
        # Count all second fields
        RAWCOUNT="\$(awk '/^[^#]/ { print \$2 }' \$BLIST | wc -l)"
        # Count all second fields that only contain integers
        INTCOUNT="\$(awk '\$2 ~ "^[0-9][0-9]*\$" { print \$2 }' \$BLIST | wc -l)"
        # Check that counts are equal
        if [ "\$RAWCOUNT" != "\$INTCOUNT" ]; then
            echo -e "\$BLIST second field is not numeric" >> reference_errors.txt
        fi

        # Check all third fields are integers
        RAWCOUNT="\$(awk '/^[^#]/ { print \$3 }' \$BLIST | wc -l)"
        INTCOUNT="\$(awk '\$3 ~ "^[0-9][0-9]*\$" { print \$3 }' \$BLIST | wc -l)"
        if [ "\$RAWCOUNT" != "\$INTCOUNT" ]; then
            echo -e "\$BLIST third field is not numeric" >> reference_errors.txt
        fi

        # Get species prefix from blocklist filename
        SPECIES=\$(echo \$BLIST | awk -F '[.]' '{print \$1}' | sed 's;-;_;g')

        # Paste species name onto contigs in blocklist
        sed "s/^/\$SPECIES./" \$BLIST >> mixed.blocklist.bed

    done

    ECOUNT="\$(less reference_errors.txt | wc -l)"
    if [ "\$ECOUNT" != 0 ]; then
        echo "Errors detected in reference formatting"
        cat reference_errors.txt
        exit 1
    fi
    """
}
