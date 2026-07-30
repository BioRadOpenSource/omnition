/*
Prepare single-cell 3' RNA Droplet analysis reference files
*/

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF               } from "${params.omnition.modulesDir}/public/core/gunzip.nf"
include { VERIFY_REFERENCES            } from "${params.omnition.modulesDir}/public/rna/verify_references.nf"
include { COMBINE_REFERENCES           } from "${params.omnition.modulesDir}/public/core/combine_references.nf"
include { FILTER_REFERENCES            } from "${params.omnition.modulesDir}/public/core/filter_references.nf"
include { VALIDATE_REFERENCES          } from "${params.omnition.modulesDir}/public/rna/validate_references.nf"
include { MIXED_GENE_SYMBOLS           } from "${params.omnition.modulesDir}/public/rna/mixed_gene_symbols.nf"
include { STAR_INDEX                   } from "${params.omnition.modulesDir}/public/rna/star_index.nf"
include { MAKE_SAF                     } from "${params.omnition.modulesDir}/public/rna/make_saf.nf"
include { MAKE_REFFLAT                 } from "${params.omnition.modulesDir}/public/rna/make_refflat.nf"
include { MAKE_RIBOSOMAL_INTERVAL_LIST } from "${params.omnition.modulesDir}/public/rna/make_ribosomal_interval_list.nf"

workflow RNA_REFERENCE {
    take:
    ch_reference_fasta // path: reference fasta file(s)
    ch_reference_gtf // path: reference gtf file(s)
    ch_images_pulled // boolean: true if Singularity pull completed

    main:
    // Decompress FASTA files if needed
    GUNZIP_FASTA(
            ch_reference_fasta,
            ch_images_pulled
    )

    // Decompress GTF files if needed
    GUNZIP_GTF(
            ch_reference_gtf,
            ch_images_pulled
    )

    // Verify reference formatting
    VERIFY_REFERENCES(
            GUNZIP_FASTA.out.decompressed.collect(),
            GUNZIP_GTF.out.decompressed.collect(),
            ch_images_pulled
    )

    // Prepare mixed species references
    if (params.rna.mixed) {
        COMBINE_REFERENCES(
            GUNZIP_FASTA.out.decompressed.collect().sort(),
            GUNZIP_GTF.out.decompressed.collect().sort(),
            ch_images_pulled
       )
        ch_reference_fasta_prep = COMBINE_REFERENCES.out.fasta
        ch_reference_gtf_prep = COMBINE_REFERENCES.out.gtf
    } else {
        ch_reference_fasta_prep = GUNZIP_FASTA.out.decompressed.collect().first()
        ch_reference_gtf_prep = GUNZIP_GTF.out.decompressed.first()
    }

    // Filter GTF records to match provided FASTA
    FILTER_REFERENCES(
        ch_reference_fasta_prep,
        ch_reference_gtf_prep,
        ch_images_pulled
   )

    // Verify that FASTA and GTF match each other
    VALIDATE_REFERENCES(
        ch_reference_fasta_prep,
        FILTER_REFERENCES.out.gtf,
        ch_images_pulled
   )
    ch_final_reference_fasta = VALIDATE_REFERENCES.out.fasta
    ch_final_reference_gtf = VALIDATE_REFERENCES.out.gtf

    // Prepare mixed species gene symbols files
    if (params.rna.mixed) {
        MIXED_GENE_SYMBOLS(
            ch_final_reference_gtf,
            ch_images_pulled
       )
        ch_mixed_symbols = MIXED_GENE_SYMBOLS.out.symbols
    } else {
        ch_mixed_symbols = Channel.empty()
    }

    // Index references with STAR
    STAR_INDEX(
        ch_final_reference_fasta,
        ch_final_reference_gtf,
        ch_images_pulled
   )

    // Create refFlat-formatted reference file
    MAKE_REFFLAT(
        ch_final_reference_gtf,
        ch_images_pulled
   )

    // Create SAF-formatted reference file
    MAKE_SAF(
        ch_final_reference_gtf,
        STAR_INDEX.out.index,
        ch_images_pulled
   )

    // Create list of ribosomal feature locations
    MAKE_RIBOSOMAL_INTERVAL_LIST(
        ch_final_reference_gtf,
        STAR_INDEX.out.index,
        ch_images_pulled
   )

    emit:
    reference_fasta         = ch_final_reference_fasta // path: params.fasta
    reference_gtf           = ch_final_reference_gtf // path: *.filtered.gtf
    reference_mixed_symbols = ch_mixed_symbols // path {species}_gene_symbols.txt
    reference_index         = STAR_INDEX.out.index // path: star-index/
    reference_saf           = MAKE_SAF.out.saf // path: annotation.saf
    reference_symbols       = MAKE_SAF.out.symbols // path: gene_symbols.txt
    reference_refflat       = MAKE_REFFLAT.out.refflat // path: annotation.refflat
    reference_interval_list = MAKE_RIBOSOMAL_INTERVAL_LIST.out.interval_list // path: ribosomal.interval_list
}
