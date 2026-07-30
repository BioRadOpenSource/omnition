/*
Check Singularity images
*/

include { SINGULARITY_PULL } from '../../../modules/public/core/singularity_pull.nf'

workflow PULL_SINGULARITY_IMAGES {
    take:

    ch_workflow_modules // Modules used in a workflow

    main:

    SINGULARITY_PULL(
        ch_workflow_modules.filter { !(it =~ /singularity_pull.nf/) }.collect()
    )

    emit:
        images_pulled = SINGULARITY_PULL.out.completed
}
