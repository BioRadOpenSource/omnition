/*
Analyze single-cell Cytometry Droplet data for CITE-seq
*/

include { KALLISTO                } from "${params.omnition.modulesDir}/public/cytometry/kallisto.nf"
include { BUSTOOLS                } from "${params.omnition.modulesDir}/public/cytometry/bustools.nf"
include { CELL_FILTERING_MERGING  } from "${params.omnition.modulesDir}/public/citeseq/cell_filtering_merging.nf"
include { MAS_MATRIX_CITE         } from "${params.omnition.modulesDir}/public/citeseq/mas_matrix_cite.nf"
include { SEQUENCE_SEARCH         } from "${params.omnition.modulesDir}/public/cytometry/sequence_search.nf"
include { ANALYSIS                } from "${params.omnition.modulesDir}/public/cytometry/analysis.nf"
include { ANTIBODY_COUNTS         } from "${params.omnition.modulesDir}/public/cytometry/antibody_counts.nf"
include { PACK_SINGLE_H5AD        } from "${params.omnition.modulesDir}/public/core/pack_single_h5ad.nf"
include { PACK_BATCH_H5AD         } from "${params.omnition.modulesDir}/public/core/pack_batch_h5ad.nf"
include { PACK_MULTIOMICS_H5AD    } from "${params.omnition.modulesDir}/public/citeseq/pack_multiomics_h5ad.nf"
include { AGGREGATE_METRICS       } from "${params.omnition.modulesDir}/public/core/aggregate_metrics.nf"

workflow CITESEQ_ANALYSIS {
    take:
    ch_debarcoded_fastq // tuple: [ sampleId, debarcoded R1 FASTQ files, raw R2 FASTQ files ]
    ch_params_file // path: input params file
    ch_reference_adt_file // path: prepared reference antibody file
    ch_reference_kite_t2g // path: prepared reference kite t2g file
    ch_reference_kallisto_index // path: prepared reference kallisto index
    ch_reference_seqsearch_index // path: prepared reference seqsearch index
    ch_cell_calling_results // tuple: [ sampleId, csvs with cell calling results ]
    ch_reference_allowlist // tuple: [ sampleId, allowlist file ]
    ch_reference_translate // tuple: [ sampleId, translate file ]
    ch_rna_filtered_mtx // tuple: [ sampleId, mtx, barcodes, genes files ]
    ch_debarcoder_counts // tuple: [ sampleId, debarcoder counts csv ]
    ch_fastqc_zip // tuple: [ sampleId, zip of fastqc results ]
    ch_fastqc_counts // tuple: [ sampleId, fastqc read counts csv ]
    ch_rna_batch_h5ad // path: rna batch h5ad file
    assay // string: the name of the assay
    ch_sample_sheet_map // map: sample sheet metadata
    ch_messages // any warning or info messages created by the pipeline
    ch_images_pulled // boolean: true if Singularity pull completed

    main:
    // Organizing input FASTQ files as RNA does and using the same channel name for consistency
    ch_complete_fastq = ch_debarcoded_fastq
        .map { out -> out.flatten() }
        .map { out -> tuple(out[0], [out[1], out[2]]) }

    // Create file with sample IDs
    ch_samples = ch_debarcoded_fastq
        .map { i -> i[0] }
        .collectFile(name: 'samples.txt', newLine: true, cache: 'lenient')

    // Mapping reads to antibody tags
    KALLISTO(
        ch_debarcoded_fastq
            .combine(ch_reference_kallisto_index),
        ch_images_pulled
    )

    // Generating the count matrix
    BUSTOOLS(
        KALLISTO.out.kallisto
            .combine(ch_reference_kite_t2g),
        ch_images_pulled
    )

    // Cell filtering and merging using RNA results for ADT reads
    CELL_FILTERING_MERGING(
        BUSTOOLS.out.bustools_counts
            .join(ch_reference_allowlist)
            .join(ch_reference_translate),
        ch_images_pulled
    )

    // Create output matrix for use with MAS
    MAS_MATRIX_CITE(
        CELL_FILTERING_MERGING.out.matrix
            .join(ch_rna_filtered_mtx, by:0),
        ch_images_pulled
    )

    // Checking R2 for known contaminants: mito, ribosomal, unused adts
    SEQUENCE_SEARCH(
        ch_debarcoded_fastq.map { i -> tuple(i[0], i[1][1]) }
            .combine(ch_reference_adt_file)
            .combine(ch_reference_seqsearch_index),
        ch_images_pulled
    )

    // Collecting the analysis data and generating figures for the final report
    ANALYSIS(
        BUSTOOLS.out.counter_type,
        ch_cell_calling_results
            .join(CELL_FILTERING_MERGING.out.cell_calling_filtering_results, by:0)
            .join(CELL_FILTERING_MERGING.out.initial_cell_count, by:0),
        ch_images_pulled
    )

    // Counting the antibodies per cell in the fully processed cells
    ANTIBODY_COUNTS(
        CELL_FILTERING_MERGING.out.filtered_cells,
        ch_images_pulled
    )

    // Adding antibody counts in a future release - create placeholder channel
    ch_summarize_expression_count = CELL_FILTERING_MERGING.out.matrix
        .map { tuple(it[0], file('NO_FILE')) }

    // Aggregate data into an H5AD format
    PACK_SINGLE_H5AD(
        CELL_FILTERING_MERGING.out.matrix
            .join(ch_summarize_expression_count)
            .join(ANALYSIS.out.metadata),
        ch_images_pulled
    )

    // Combine all H5AD files together into a single file
    PACK_BATCH_H5AD(
        PACK_SINGLE_H5AD.out.h5ad
            .map { i -> i[1] }
            .collect()
            .ifEmpty { [] },
        ch_images_pulled
    )

    // Combine all H5AD files together into a single file
    PACK_MULTIOMICS_H5AD(
        ch_rna_batch_h5ad,
        PACK_BATCH_H5AD.out.h5ad,
        ch_images_pulled
    )

    // Combine useful metrics into a single file
    AGGREGATE_METRICS(
        Channel.empty()
            .mix(ch_debarcoder_counts,
                ch_fastqc_zip,
                ch_fastqc_counts,
                SEQUENCE_SEARCH.out.metrics,
                CELL_FILTERING_MERGING.out.filtering_metrics)
            .flatMap { i -> i[1] }
            .collect(),
        assay,
        ch_images_pulled
    )
}
