/*
Indexing antibody tags for mapping
*/

process KALLISTO_INDEX {
    container "${params.container_image.cytometry_public}:${workflow.manifest.version}"
    publishDir "${params.options.reference.directory}/", mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_medium'

    input:
    path antibody_file_cleaned
    path kite_fasta
    val images_pulled

    output:
    path 'adt_index.idx', emit: kallisto_index

    script:
    """
    TAGLENGTH=\$(awk -F ',' 'NR == 1 || (length(\$1) > 0 && length(\$1) < shortest) \
        {shortest = length(\$1)} END {print shortest}' ${antibody_file_cleaned})
    if ! (( \$TAGLENGTH % 2 )); then
        echo 'The min tags length was even, subtracting one to make it odd...';
        TAGLENGTH=\$((\$TAGLENGTH - 1));
    fi
    kallisto index -i adt_index.idx -k \$TAGLENGTH ${kite_fasta}
    """
}
