/*
Formatting single species blocklists
*/

params.options = [:]

process FORMAT_BLOCKLIST {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    errorStrategy 'terminate'
    label 'cpu_xsmall'
    label 'memory_xxsmall'

    input:
    path blocklist
    val images_pulled

    output:
    path 'formatted.*.bed', emit: bed

    script:
    """
    touch reference_errors.txt
    # Check all second fields are integers
    # Count all second fields
    RAWCOUNT="\$(awk '/^[^#]/ { print \$2 }' ${blocklist} | wc -l)"
    # Count all second fields that only contain integers
    INTCOUNT="\$(awk '\$2 ~ "^[0-9][0-9]*\$" { print \$2 }' ${blocklist} | wc -l)"
    # Check that counts are equal
    if [ "\$RAWCOUNT" != "\$INTCOUNT" ]; then
        echo -e "${blocklist} second field is not numeric" >> reference_errors.txt
    fi

    # Check all third fields are integers
    RAWCOUNT="\$(awk '/^[^#]/ { print \$3 }' ${blocklist} | wc -l)"
    INTCOUNT="\$(awk '\$3 ~ "^[0-9][0-9]*\$" { print \$3 }' ${blocklist} | wc -l)"
    if [ "\$RAWCOUNT" != "\$INTCOUNT" ]; then
        echo -e "${blocklist} third field is not numeric" >> reference_errors.txt
    fi

    ECOUNT="\$(less reference_errors.txt | wc -l)"
    if [ "\$ECOUNT" != 0 ]; then
        echo "Errors detected in reference formatting"
        cat reference_errors.txt
        exit 1
    fi

    # Get species prefix from blocklist filename
    SPECIES=\$(echo ${blocklist} | awk -F '[.]' '{print \$1}' | sed 's;-;_;g')

    # Paste species name onto contigs in blocklist
    sed "s/^/\$SPECIES./" ${blocklist} >> formatted.${blocklist}

    """
}
