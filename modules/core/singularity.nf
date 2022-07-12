/*
Check ATAC containers
*/

process SINGULARITY_PULL {
    input:
    val module_files

    output:
    val true, emit: completed

    script:
    """
    modules=\$(echo ${module_files} | sed 's/[][]//g' | tr -d \\,)

    mkdir -p ${projectDir}/singularity/

    # Parallel Singularity image pulls create race conditions, do not parallelize this process
    for FILE in \$modules; do

        echo \$FILE

        DOCKER=\$(grep -o -m 1 'container.*' \$FILE | cut -f2- -d\\" | cut -f1 -d:)
        SINGULARITY="\$(echo \$DOCKER | sed 's;[/:];-;g')"

        # Some modules dont use containers, this if passes for modules that do
        if [[ ! -z "\$DOCKER" ]]; then
            if [[ ! -f ${projectDir}/singularity/\$SINGULARITY-${workflow.manifest.version}.img ]]; then
                singularity pull -F ${projectDir}/singularity/\$SINGULARITY-${workflow.manifest.version}.img \
                    docker://\$DOCKER:${workflow.manifest.version}
            fi
        fi

    done
    """
}
