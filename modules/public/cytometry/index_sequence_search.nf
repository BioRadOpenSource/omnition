/*
Indexing the R2 antibody tag regions for mapping
*/

process INDEX_SEQUENCE_SEARCH {
    container "${params.container_image.cytometry_public}:${workflow.manifest.version}"
    publishDir "${params.options.reference.directory}/", mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_medium'

    input:
    path cleaned
    path fasta
    val images_pulled

    output:
    path 'seqsearch_index.idx', emit: seqsearch_index

    script:
    """
    TAGLENGTH=\$(awk -F ',' 'NR == 1 || (length(\$1) > 0 && length(\$1) < shortest) \
        {shortest = length(\$1)} END {print shortest}' ${cleaned})
    if ! (( \$TAGLENGTH % 2 )); then
        echo 'Minimum tag length was even, subtracting one to make it odd...';
        TAGLENGTH=\$((\$TAGLENGTH - 1));
    fi

    cp /opt/biorad/assets/references/* ./

    cat ${fasta} mitochondrial_transcripts.fasta ribosomal_genome.fasta > seqsearch.fasta
    kallisto index -i seqsearch_index.idx -k \$TAGLENGTH seqsearch.fasta
    """
}
