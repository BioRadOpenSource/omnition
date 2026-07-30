/*
Cleaning up the antibody tag file by converting to unix format and removing special characters
*/

process PREPARE_ANTIBODY_FILE {
    container "${params.container_image.cytometry_public}:${workflow.manifest.version}"
    publishDir "${params.options.reference.directory}/", mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_medium'

    input:
    path antibody_file
    val images_pulled

    output:
    path 'antibodies.csv', emit: cleaned
    path 'antibodies.json', emit: adt_debarcoder_ref

    script:
    """
    dos2unix -R -n ${antibody_file} antibodies.csv
    publicCytometryCleanAntibody.py -f antibodies.csv
    publicCytometryAdtDebarcoderRef.py antibodies.csv
    """
}
