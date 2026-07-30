/*
Prepare atac analysis reference files
*/

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF           } from "${params.omnition.modulesDir}/public/core/gunzip.nf"
include { VERIFY_REFERENCES        } from "${params.omnition.modulesDir}/public/atac/verify_references.nf"
include { COMBINE_REFERENCES       } from "${params.omnition.modulesDir}/public/core/combine_references.nf"
include { GENERATE_EMPTY_BLOCKLIST } from "${params.omnition.modulesDir}/public/atac/generate_empty_blocklist.nf"
include { FORMAT_BLOCKLIST         } from "${params.omnition.modulesDir}/public/atac/format_blocklist.nf"
include { COMBINE_BLOCKLISTS       } from "${params.omnition.modulesDir}/public/atac/combine_blocklists.nf"
include { FILTER_REFERENCES        } from "${params.omnition.modulesDir}/public/core/filter_references.nf"
include { FILTER_BLOCKLISTS        } from "${params.omnition.modulesDir}/public/atac/filter_blocklists.nf"
include { BWA_INDEX                } from "${params.omnition.modulesDir}/public/atac/bwa_index.nf"
include { ARCHR_REFERENCE          } from "${params.omnition.modulesDir}/public/atac/archr_reference.nf"
include { GENERATE_GENOME_SIZES    } from "${params.omnition.modulesDir}/public/atac/generate_genome_sizes.nf"
include { GENERATE_TSS_WINDOWS     } from "${params.omnition.modulesDir}/public/atac/generate_tss_windows.nf"

workflow ATAC_NORMAL_REFERENCE {
    take:
    ch_fasta // path: reference fasta file(s)
    ch_gtf // path: reference gtf file(s)
    ch_images_pulled // boolean: true if Singularity pull completed

    main:
    // Decompress FASTA files if needed
    GUNZIP_FASTA(
            ch_fasta,
            ch_images_pulled
    )

    // Decompress GTF files if needed
    GUNZIP_GTF(
            ch_gtf,
            ch_images_pulled
    )

    // Verify reference formatting
    VERIFY_REFERENCES(
            GUNZIP_FASTA.out.decompressed.collect(),
            GUNZIP_GTF.out.decompressed.collect(),
            ch_images_pulled
    )

    // Make mixed references if needed
    if (params.atac.mixed) {
        COMBINE_REFERENCES(
            GUNZIP_FASTA.out.decompressed.collect().sort(),
            GUNZIP_GTF.out.decompressed.collect().sort(),
            ch_images_pulled
       )
        ch_fasta_prep = COMBINE_REFERENCES.out.fasta
        ch_gtf_prep = COMBINE_REFERENCES.out.gtf
    } else {
        ch_fasta_prep = GUNZIP_FASTA.out.decompressed.first()
        ch_gtf_prep = GUNZIP_GTF.out.decompressed.first()
    }

    // Create the blocklist channel
    if (!params.atac.reference.blocklist) {
        // If no blocklist is supplied make an empty one
        GENERATE_EMPTY_BLOCKLIST(
            ch_images_pulled
       )
        ch_blocklist = GENERATE_EMPTY_BLOCKLIST.out.bed
    } else {
        block_count = Core.countElement(params.atac.reference.blocklist)
        if (block_count == 1) {
            // If one blocklist is supplied format it
            ch_blocklist_pre = Channel.fromPath(params.atac.reference.blocklist).first()
            FORMAT_BLOCKLIST(
                ch_blocklist_pre,
                ch_images_pulled
           )
            ch_blocklist = FORMAT_BLOCKLIST.out.bed
        } else if (block_count == 2) {
            // If two blocklists are supplied format and combine them
            ch_blocklist_pre = Channel.fromPath(params.atac.reference.blocklist)
            COMBINE_BLOCKLISTS(
                ch_blocklist_pre.collect().sort(),
                ch_images_pulled
           )
            ch_blocklist = COMBINE_BLOCKLISTS.out.bed
        }
    }

    // Filter GTF records to match provided FASTA
    FILTER_REFERENCES(
        ch_fasta_prep,
        ch_gtf_prep,
        ch_images_pulled
   )

    // Filter blocklist to match provided FASTA
    FILTER_BLOCKLISTS(
        ch_fasta_prep,
        ch_blocklist,
        ch_images_pulled
   )

    ch_final_fasta = ch_fasta_prep
    ch_final_gtf = FILTER_REFERENCES.out.gtf

    // Indexing the reference genome for read mapping with STA
    BWA_INDEX(
        ch_final_fasta,
        ch_images_pulled
   )

    // Generate BSgenome package for ArchR reference
    ARCHR_REFERENCE(
        ch_final_fasta,
        ch_final_gtf,
        ch_images_pulled
   )

    // Get contig sizes
    GENERATE_GENOME_SIZES(
        ch_final_fasta,
        ch_images_pulled
   )

    // Generate TSS window
    GENERATE_TSS_WINDOWS(
        ch_final_gtf,
        ch_final_fasta,
        GENERATE_GENOME_SIZES.out.size,
        ch_images_pulled
   )

    emit:
        fasta         = ch_final_fasta // path: params.fasta
        gtf           = ch_final_gtf // path: *.filtered.gtf
        tss           = GENERATE_TSS_WINDOWS.out.tss
        bwaIndex      = BWA_INDEX.out.index
        genomeSize    = GENERATE_GENOME_SIZES.out.size
        archrRef      = ARCHR_REFERENCE.out.pkg
        blocklist     = FILTER_BLOCKLISTS.out.blocklist
}
