/*
Check Singularity images
*/

include { SINGULARITY_PULL        } from '../../modules/core/singularity.nf'

workflow PULL_IMAGES {
    take:

    ch_workflow_modules // Modules used in a workflow

    main:

    SINGULARITY_PULL(
        ch_workflow_modules.filter { !(it =~ /singularity.nf/) }.collect()
    )

    emit:
        images_pulled       = SINGULARITY_PULL.out.completed
}
