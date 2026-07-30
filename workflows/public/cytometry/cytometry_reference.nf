/*
Prepare single-cell cytometry Droplet analysis reference files
*/

include { VALIDATE_ANTIBODY_FILE } from "${params.omnition.modulesDir}/public/cytometry/validate_antibody_file.nf"
include { KITE                   } from "${params.omnition.modulesDir}/public/cytometry/kite.nf"
include { KALLISTO_INDEX         } from "${params.omnition.modulesDir}/public/cytometry/kallisto_index.nf"
include { INDEX_SEQUENCE_SEARCH  } from "${params.omnition.modulesDir}/public/cytometry/index_sequence_search.nf"

workflow CYTOMETRY_REFERENCE {
    take:
    ch_prepared_antibody_file
    ch_images_pulled

    main:
    // Check that antibody tags are distinguishable
    VALIDATE_ANTIBODY_FILE(
        ch_prepared_antibody_file,
        ch_images_pulled
   )

    // Create map of antibody tag sequences and mutations
    KITE(
        ch_prepared_antibody_file,
        VALIDATE_ANTIBODY_FILE.out.validated,
        ch_images_pulled
   )

    // Indexing antibody tags for mapping
    KALLISTO_INDEX(
        ch_prepared_antibody_file,
        KITE.out.fasta,
        ch_images_pulled
   )

    // Indexing the R2 antibody tag regions for mapping
    INDEX_SEQUENCE_SEARCH(
        ch_prepared_antibody_file,
        KITE.out.fasta,
        ch_images_pulled
   )

    emit:
        reference_kite_t2g              = KITE.out.t2g // path: features_processed.t2g
        reference_kallisto_index        = KALLISTO_INDEX.out.kallisto_index // path: adt_index.idx
        reference_seqsearch_index       = INDEX_SEQUENCE_SEARCH.out.seqsearch_index // path: seqsearch_index.idx
}
