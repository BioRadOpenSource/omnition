/*
Validate parameters and execute core workflows
*/

// Import workflows
include { DO_MERGING_ANALYSIS           } from "${params.omnition.workflowsDir}/public/do_merging/do_merging_analysis.nf"

workflow DO_MERGING {
    take:
    ch_images_pulled // boolean: images pulled indicator
    ch_val_pass // tuple: [ sampleId, empty file for validation ]
    ch_debarcoder_fastq // tuple: [ sampleId, debarcoded R1 FASTQ files, debarcoded R2 FASTQ files ]
    ch_debarcoder_edges // tuple: [ sampleId, debarcoded edges file ]
    ch_bead_list // tuple: [ sampleId, bead list file ]
    messages // Channel for messages

    main:
    // Set global assay params
    params.do_merging.assay      = "DO Merging"
    params.do_merging.sampleIds  = params.core.sampleIds
    params.do_merging.resultsDir = "${params.core.outputDir}/Sample_Files/"
    params.do_merging.reportsDir = "${params.core.outputDir}/report/"

    PublicCore.loadSampleLevelParams(params.do_merging)

    // Execute workflow
    DO_MERGING_ANALYSIS(
        ch_val_pass,
        ch_debarcoder_fastq,
        ch_debarcoder_edges,
        ch_bead_list,
        ch_images_pulled
    )

    emit:
    edges_input = DO_MERGING_ANALYSIS.out.edges_input
    corrected_edges = DO_MERGING_ANALYSIS.out.corrected_edges
    barcode_translate = DO_MERGING_ANALYSIS.out.barcode_translate
}
