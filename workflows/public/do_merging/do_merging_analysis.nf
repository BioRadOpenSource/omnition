/*
Analyze single-cell data with do_merging modules
*/

include { CORRECT_EDGES                                   } from "${params.omnition.modulesDir}/public/rna/correct_edges.nf"
include { MERGE_BEADS                                     } from "${params.omnition.modulesDir}/public/rna/merge_beads.nf"
include { CORRECT_CELL_BARCODES                           } from "${params.omnition.modulesDir}/public/core/correct_cell_barcodes.nf"

workflow DO_MERGING_ANALYSIS {
    take:
    ch_val_pass // tuple: [ sampleId, empty file for validation ]
    ch_debarcoded_fastq // tuple: [ sampleId, debarcoded R1 FASTQ files, debarcoded R2 FASTQ files ]
    ch_debarcoded_edges // tuple: [ sampleId, debarcoded edges file ]
    ch_bead_list // tuple: [ sampleId, barcodes file ]
    ch_images_pulled // boolean: true if Singularity pull completed

    main:
    ch_edges_input = ch_debarcoded_edges

    // Correct DOs and UMIs
    CORRECT_EDGES(
        ch_edges_input,
        ch_images_pulled
    )

    // Merge beads
    MERGE_BEADS(
        CORRECT_EDGES.out.edges
            .join(ch_bead_list, by: 0),
        ch_images_pulled
    )

    // Correct barcode formatting
    CORRECT_CELL_BARCODES(
        MERGE_BEADS.out.barcode_translate
            .join(MERGE_BEADS.out.filtered_beads),
        ch_images_pulled
    )

    emit:
    edges_input = ch_edges_input
    corrected_edges = CORRECT_EDGES.out.edges
    barcode_translate = CORRECT_CELL_BARCODES.out.barcode_translate
}
