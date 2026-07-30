/*
Checking that antibody tags are distinguishable enough from each other
*/

process VALIDATE_ANTIBODY_FILE {
    container "${params.container_image.cytometry_public}:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    path antibody_file_cleaned
    val images_pulled

    output:
    val true, emit: validated

    script:
    """
    MINLEVDIST=\$(publicCytometryMinimumLevDistanceCheck.py ${antibody_file_cleaned})
    if [ \$MINLEVDIST \\< 4 ]; then
        echo "Warning: Minimum Levenshtein distance between antibodies is only \$MINLEVDIST";
    fi;
    """
}
