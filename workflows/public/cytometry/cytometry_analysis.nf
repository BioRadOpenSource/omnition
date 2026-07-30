/*
Analyze single-cell Cytometry Droplet data
*/

include { COMPILE_QC_REPORT       } from "${params.omnition.modulesDir}/public/cytometry/compile_qc_report.nf"
include { KALLISTO                } from "${params.omnition.modulesDir}/public/cytometry/kallisto.nf"
include { BUSTOOLS                } from "${params.omnition.modulesDir}/public/cytometry/bustools.nf"
include { MERGE_BARCODES          } from "${params.omnition.modulesDir}/public/cytometry/merge_barcodes.nf"
include { CELL_CALLING            } from "${params.omnition.modulesDir}/public/rna/cell_calling.nf"
include { SEQUENCE_SEARCH         } from "${params.omnition.modulesDir}/public/cytometry/sequence_search.nf"
include { CELL_CALLING_FILTERING  } from "${params.omnition.modulesDir}/public/cytometry/cell_calling_filtering.nf"
include { MAS_MATRIX_CYTO         } from "${params.omnition.modulesDir}/public/cytometry/mas_matrix_cyto.nf"
include { ANALYSIS                } from "${params.omnition.modulesDir}/public/cytometry/analysis.nf"
include { ANTIBODY_COUNTS         } from "${params.omnition.modulesDir}/public/cytometry/antibody_counts.nf"
include { PACK_SINGLE_H5AD        } from "${params.omnition.modulesDir}/public/core/pack_single_h5ad.nf"
include { PACK_BATCH_H5AD         } from "${params.omnition.modulesDir}/public/core/pack_batch_h5ad.nf"
include { AGGREGATE_METRICS       } from "${params.omnition.modulesDir}/public/core/aggregate_metrics.nf"

workflow CYTOMETRY_ANALYSIS {
    take:
    ch_input_params_yaml // path: input params yaml file
    ch_params_yaml // path: params yaml file
    ch_command_txt // path: command txt file
    ch_fastqc_zip // path: fastqc zip file
    ch_fastqc_count // path: fastqc count file
    ch_complete_fastq // tuple: [ sampleId, raw R1 FASTQ files, raw R2 FASTQ files ]
    ch_debarcoded_fastq // tuple: [ sampleId, debarcoded R1 FASTQ files, debarcoded R2 FASTQ files ]
    ch_debarcoded_count // tuple: [ sampleId, debarcoded counts ]
    ch_barcode_translate // tuple: [ sampleId, barcode translate file ]
    ch_prepared_antibody_file // path: prepared reference antibody file
    ch_reference_kite_t2g // path: prepared reference kite t2g file
    ch_reference_kallisto_index // path: prepared reference kallisto index
    ch_reference_seqsearch_index // path: prepared reference seqsearch index
    assay // string: the name of the assay
    ch_sample_sheet_map // map: sample sheet metadata
    ch_messages // any warning or info messages created by the pipeline
    ch_images_pulled // boolean: true if Singularity pull completed

    main:
    // Create file with sample IDs
    ch_samples = ch_complete_fastq
        .map { i -> i[0] }
        .collectFile(name: 'samples.txt', newLine: true, cache: 'lenient')

    // Combining library QC reports into a single report
    COMPILE_QC_REPORT(
        ch_fastqc_zip.flatMap { i -> i[1] } .collect(),
        ch_images_pulled
    )

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

    // Use bead merge list to merge count matrix
    MERGE_BARCODES(
        BUSTOOLS.out.bustools_counts
            .join(ch_barcode_translate),
        ch_images_pulled
    )

    // Calculating the proposed number of cells present
    CELL_CALLING(
        MERGE_BARCODES.out.merged_counts,
        ch_images_pulled
    )

    // Checking R2 for known contaminants: mito, ribosomal, unused adts
    SEQUENCE_SEARCH(
        ch_debarcoded_fastq.map { i -> tuple(i[0], i[1][1]) }
            .combine(ch_prepared_antibody_file)
            .combine(ch_reference_seqsearch_index),
        ch_images_pulled
    )

    // Cell Calling Filtering Process
    CELL_CALLING_FILTERING(
        CELL_CALLING.out.allowlist
            .join(MERGE_BARCODES.out.old_merged_counts, by:0),
        ch_images_pulled
    )

    // Create output matrix for use with MAS
    MAS_MATRIX_CYTO(
        CELL_CALLING_FILTERING.out.matrix,
        ch_images_pulled
    )

    // Collecting the analysis data and generating figures for the final report
    ANALYSIS(
        BUSTOOLS.out.counter_type,
        CELL_CALLING.out.results
            .join(CELL_CALLING_FILTERING.out.cell_calling_filtering_results, by:0)
            .join(CELL_CALLING_FILTERING.out.initial_cell_count, by:0),
        ch_images_pulled
    )

    // Counting the antibodies per cell in the fully processed cells
    ANTIBODY_COUNTS(
        CELL_CALLING_FILTERING.out.filtered_cells,
        ch_images_pulled
    )

    // Adding antibody counts in a future release - create placeholder channel
    ch_summarize_expression_count = CELL_CALLING_FILTERING.out.matrix
        .map { tuple(it[0], file('NO_FILE')) }

    // Aggregate data into an H5AD format
    PACK_SINGLE_H5AD(
        CELL_CALLING_FILTERING.out.matrix
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

    // Combine useful metrics into a single file
    AGGREGATE_METRICS(
        Channel.empty()
            .mix(ch_fastqc_zip,
                ch_fastqc_count,
                SEQUENCE_SEARCH.out.metrics,
                CELL_CALLING_FILTERING.out.filtering_metrics,
                CELL_CALLING.out.all_but_allowlist,
                CELL_CALLING.out.allowlist,
                ch_debarcoded_count)
            .flatMap { i -> i[1] }
            .collect(),
        assay,
        ch_images_pulled
    )
}
