/*
Verify reference formatting
*/

params.options = [:]

process VERIFY_REFERENCES {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    errorStrategy 'terminate'
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    path fasta
    path gtf
    val images_pulled

    script:
    """
    touch reference_errors.txt

    # Check if reference fastas have different headings for each contig
    for FILE in ${fasta}; do
        # Count all contig headers
        RAWCOUNT="\$(grep ">" \$FILE | wc -l)"
        # Count unique contig headers
        UNIQCOUNT="\$(grep ">" \$FILE | sort | uniq | wc -l)"
        # Check that counts are equal
        if [ "\$RAWCOUNT" != "\$UNIQCOUNT" ]; then
            echo -e "Contig names are not unique" >> reference_errors.txt
        fi

    done

    for FILE in ${gtf}; do
        # Check that the third field has transcript records in it
        # '/^[^#]/ - skip lines beginning with # (header)
        # { print \$3 } - check third field
        # grep -m 1 "transcript" - look for transcript, stop at first record
        # || true - if grep doesn't find anything it returns exit code 1
        FINDTEST="\$(awk '/^[^#]/ { print \$3 }' \$FILE | grep -m 1 "transcript" || true )"
        if [[ -z \$FINDTEST ]]; then
            echo -e "\$FILE is missing transcript records in the third field" >> reference_errors.txt
        fi

        # Check that the third field has exon records in it
        FINDTEST="\$(awk '/^[^#]/ { print tolower(\$3) }' \$FILE | grep -m 1 "exon" || true )"
        if [[ -z \$FINDTEST ]]; then
            echo -e "\$FILE is missing exon records in the third field" >> reference_errors.txt
        fi

        # Check all fourth fields are integers
        # Count all fourth fields
        RAWCOUNT="\$(awk '/^[^#]/ { print \$4 }' \$FILE | wc -l)"
        # Count all fourth fields that only contain integers
        INTCOUNT="\$(awk '\$4 ~ "^[0-9][0-9]*\$" { print \$4 }' \$FILE | wc -l)"
        # Check that counts are equal
        if [ "\$RAWCOUNT" != "\$INTCOUNT" ]; then
            echo -e "\$FILE fourth field is not numeric" >> reference_errors.txt
        fi

        # Check all fifth fields are integers
        RAWCOUNT="\$(awk '/^[^#]/ { print \$5 }' \$FILE | wc -l)"
        INTCOUNT="\$(awk '\$5 ~ "^[0-9][0-9]*\$" { print \$5 }' \$FILE | wc -l)"
        if [ "\$RAWCOUNT" != "\$INTCOUNT" ]; then
            echo -e "\$FILE fifth field is not numeric" >> reference_errors.txt
        fi

        # Check that the seventh field has strand information in it
        FINDTEST="\$(awk '/^[^#]/ { print \$7 }' \$FILE | grep -m 1 "+" || true )"
        if [[ -z \$FINDTEST ]]; then
            echo -e "\$FILE is missing strand information in the seventh field" >> reference_errors.txt
        fi

        # Check that records have gene id information
        FINDTEST="\$(grep -m 1 "gene_id" \$FILE || true )"
        if [[ -z \$FINDTEST ]]; then
            echo -e "\$FILE is missing gene id information" >> reference_errors.txt
        fi
    done

    ECOUNT="\$(less reference_errors.txt | wc -l)"
    if [ "\$ECOUNT" != 0 ]; then
        echo "Errors detected in reference formatting"
        cat reference_errors.txt
        exit 1
    fi
    """
}
