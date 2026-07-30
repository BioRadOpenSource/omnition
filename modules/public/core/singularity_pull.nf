/*
Pull Singularity images from Dockerhub
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
    IMAGE_VARS=\$(echo ${params.container_image} | sed 's/[][]//g')
    echo \$IMAGE_VARS

    # Parallel Singularity image pulls create race conditions, do not parallelize this process
    for FILE in \$modules; do

        echo \$FILE

        IMAGE_VAR_NEEDED=\$(grep 'params.container_image' \$FILE | sed 's/.*params.container_image.//g' | sed 's/\\}.*//g')
        echo \$IMAGE_VAR_NEEDED
        if [[ -z "\$IMAGE_VAR_NEEDED" ]]; then
            # Use the old parsing method
            DOCKER=\$(grep -o -m 1 'container.*' \$FILE | cut -f2- -d\\" | cut -f1 -d:)
            echo \$DOCKER
            SINGULARITY="\$(echo \$DOCKER | sed 's;[/:];-;g')"
            echo \$SINGULARITY
        else
            # Use the new parsing method
            DOCKER=\$(echo \$IMAGE_VARS | sed 's/, /\\n/g' | grep "^\$IMAGE_VAR_NEEDED:" | cut -f2 -d:)
            echo \$DOCKER
            SINGULARITY="\$(echo \$DOCKER | sed 's;[/:];-;g')"
            echo \$SINGULARITY
        fi
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
