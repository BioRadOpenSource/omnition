/*
Creating map of antibody tag reference sequences and potential mutations
*/

process KITE {
    container "${params.container_image.cytometry_public}:${workflow.manifest.version}"
    publishDir "${params.options.reference.directory}/", mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_medium'

    input:
    path antibody_file_cleaned
    val validated
    val images_pulled

    output:
    path 'features_processed.fasta', emit: fasta
    path 'features_processed.t2g', emit: t2g

    script:
    """
    # Oligo file, flip it and reverse it
    awk -F, '{ print \$2 "," \$1 }' ${antibody_file_cleaned} > antibodies_reverse.csv
    publicCytometryFeatureMap.py antibodies_reverse.csv --fa features_processed.fasta --t2g features_processed.t2g
    """
}
