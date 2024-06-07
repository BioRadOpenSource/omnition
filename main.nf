#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

log.info(Core.bioradLogo(params.monochrome_logs))
messages = []

// Check if help flag was given
if (params.help) {
    message = Core.helpMessage(params.monochrome_logs)
    log.info(message)
    System.exit(0)
}

// Import workflows
include { PULL_IMAGES } from './workflows/core/singularity.nf'
include { RNA } from './workflows/rna/rna.nf'
include { ATAC_NORMAL } from './workflows/atac_normal/atac_normal.nf'
include { ATAC_COMBINATORIAL } from './workflows/atac_combinatorial/atac_combinatorial.nf'

workflow {
    // Analyze 3' RNA Data
    if (params.rna) {
        RNA(
            messages
        )
    }

    // Analyze ATAC data
    if (params.atac) {
        ATAC_NORMAL(
            messages
        )
    }

    // Analyze ATAC data
    if (params.catac) {
        ATAC_COMBINATORIAL(
            messages
        )
    }
}
