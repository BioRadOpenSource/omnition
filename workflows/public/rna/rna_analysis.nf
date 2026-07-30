/*
Analyze single-cell 3' RNA Droplet data
*/

include { EDGE_METADATA                         } from "${params.omnition.modulesDir}/public/core/edge_metadata.nf"
include { STAR_ALIGN                            } from "${params.omnition.modulesDir}/public/rna/star_align.nf"
include { POST_ALIGN_PROCESSING                 } from "${params.omnition.modulesDir}/public/rna/post_align_processing.nf"
include {
    TAG_BARCODES
    TAG_BARCODES as TAG_UNMAPPED_BARCODES       } from "${params.omnition.modulesDir}/public/rna/tag_barcodes.nf"
include { TAG_MERGED_BARCODES                   } from "${params.omnition.modulesDir}/public/rna/tag_merged_barcodes.nf"
include { TAG_GENES                             } from "${params.omnition.modulesDir}/public/rna/tag_genes.nf"
include { TAG_FEATURES                          } from "${params.omnition.modulesDir}/public/rna/tag_features.nf"
include {
    BEAD_COUNTS
    BEAD_COUNTS as BEAD_COUNTS_UNMAPPED         } from "${params.omnition.modulesDir}/public/rna/bead_counts.nf"
include { COMBINE_READ_COUNTS                   } from "${params.omnition.modulesDir}/public/rna/combine_read_counts.nf"
include { UMI_TOOLS_DEDUP                       } from "${params.omnition.modulesDir}/public/rna/umi_tools_dedup.nf"
include { UMI_TOOLS_COUNT                       } from "${params.omnition.modulesDir}/public/rna/umi_tools_count.nf"
include { UNMERGE_ALLOWLIST                     } from "${params.omnition.modulesDir}/public/rna/unmerge_allowlist.nf"
include { ID_DUPLICATES                         } from "${params.omnition.modulesDir}/public/rna/id_duplicates.nf"
include { MERGE_DUPLICATES                      } from "${params.omnition.modulesDir}/public/rna/merge_duplicates.nf"
include { MAKE_COUNT_MATRIX                     } from "${params.omnition.modulesDir}/public/rna/make_count_matrix.nf"
include { CELL_CALLING                          } from "${params.omnition.modulesDir}/public/rna/cell_calling.nf"
include { FILTER_COUNT_MATRIX                   } from "${params.omnition.modulesDir}/public/rna/filter_count_matrix.nf"
include { SUMMARIZE_EXPRESSION                  } from "${params.omnition.modulesDir}/public/rna/summarize_expression.nf"
include { SPLIT_MIXED_BAM                       } from "${params.omnition.modulesDir}/public/rna/split_mixed_bam.nf"
include { SUMMARIZE_BEAD_EXPRESSION             } from "${params.omnition.modulesDir}/public/rna/summarize_bead_expression.nf"
include { AGGREGATE_BEAD_STATS                  } from "${params.omnition.modulesDir}/public/rna/aggregate_bead_stats.nf"
include { SUMMARIZE_MIXED_SPECIES               } from "${params.omnition.modulesDir}/public/rna/summarize_mixed_species.nf"
include { MIXED_EXPRESSION                      } from "${params.omnition.modulesDir}/public/rna/mixed_expression.nf"
include { SEURAT                                } from "${params.omnition.modulesDir}/public/rna/seurat.nf"
include { PACK_SINGLE_H5AD                      } from "${params.omnition.modulesDir}/public/core/pack_single_h5ad.nf"
include { PACK_BATCH_H5AD                       } from "${params.omnition.modulesDir}/public/core/pack_batch_h5ad.nf"
include {
    PICARD
    PICARD as PICARD_MIXED_SPECIES              } from "${params.omnition.modulesDir}/public/rna/picard.nf"
include { GENES_PER_SAMPLE                      } from "${params.omnition.modulesDir}/public/rna/genes_per_sample.nf"
include { COUNT_MATRIX_FEATURES                 } from "${params.omnition.modulesDir}/public/rna/count_matrix_features.nf"
include { SEQUENCE_SATURATION                   } from "${params.omnition.modulesDir}/public/core/sequence_saturation.nf"
include { SEQUENCE_SATURATION_MIXED             } from "${params.omnition.modulesDir}/public/rna/sequence_saturation_mixed.nf"
include { AGGREGATE_METRICS                     } from "${params.omnition.modulesDir}/public/core/aggregate_metrics.nf"
include { RENDER_REPORT                         } from "${params.omnition.modulesDir}/public/rna/render_report.nf"
include { CALCULATE_BEAD_PLOTS                  } from "${params.omnition.modulesDir}/public/core/calculate_bead_plots.nf"

workflow RNA_ANALYSIS {
    take:
    ch_val_pass // tuple: [ sampleId, empty file for validation ]
    ch_input_params_yaml // path: input params yaml file
    ch_params_yaml // path: params yaml file
    ch_command_txt // path: command txt file
    ch_fastqc_zip // path: fastqc zip file
    ch_fastqc_count // path: fastqc count file
    ch_fastqc_sequence_traces // path: fastqc sequence traces file
    ch_fastqc_quality_scores // path: fastqc quality scores file
    ch_complete_fastq // tuple: [ sampleId, raw R1 FASTQ files, raw R2 FASTQ files ]
    ch_cutadapt_metrics // tuple: [ sampleId, cutadapt metrics ]
    ch_debarcoded_fastq // tuple: [ sampleId, debarcoded R1 FASTQ files, debarcoded R2 FASTQ files ]
    ch_debarcoded_count // tuple: [ sampleId, debarcoded counts ]
    ch_edges_input // tuple: [ sampleId, edges file ]
    ch_corrected_edges // tuple: [ sampleId, corrected edges file ]
    ch_barcode_translate // tuple: [ sampleId, barcode translate file ]
    ch_reference_index // path: prepared reference index directory
    ch_reference_saf // path: prepared reference SAF file
    ch_reference_symbols // path: prepared reference gene symbol file
    ch_reference_refflat // path: prepared reference refFlat file
    ch_reference_interval_list // path: prepared ribosomal interval list file
    ch_mixed_symbols // paths: prepared reference gene symbol files for each species
    assay // string: the name of the assay
    ch_sample_sheet_map // map: sample sheet metadata
    ch_messages // any warning or info messages created by the pipeline
    ch_images_pulled // boolean: true if Singularity pull completed

    main:
    // Create file with sample IDs
    ch_samples = ch_complete_fastq
        .map { i -> i[0] }
        .collectFile(name: 'samples.txt', newLine: true, cache: 'lenient')

    // Create file with sample settings
    ch_settings = Channel.from(PublicCore.getSampleSettings(params.rna.bead),
        PublicCore.getSampleSettings(params.rna.barcode))
        .flatMap()
        .flatMap()
        .map { set -> String.join(',', set) }
        .collectFile(name: 'settings.txt', newLine: true, cache: 'lenient')

    ch_parsed_fastq = ch_debarcoded_fastq

    // Align reads to reference
    STAR_ALIGN(
        ch_parsed_fastq
            .combine(ch_reference_index),
        ch_images_pulled
    )

    POST_ALIGN_PROCESSING(
        ch_parsed_fastq
            .join(STAR_ALIGN.out.raw_bam)
            .join(STAR_ALIGN.out.log),
        ch_images_pulled
    )

    // Annotate reads with parse capture oligos
    TAG_BARCODES(
        POST_ALIGN_PROCESSING.out.bam,
        ch_images_pulled
    )

    // Annotate reads with parse capture oligos
    TAG_UNMAPPED_BARCODES(
        POST_ALIGN_PROCESSING.out.unmapped,
        ch_images_pulled
   )

    // Annotate reads with genomic features
    TAG_GENES(
        TAG_BARCODES.out.bam
            .combine(ch_reference_saf),
        ch_images_pulled
    )

    // Organize feature tags
    TAG_FEATURES(
        TAG_GENES.out.bam,
        ch_images_pulled
    )

    // Remove duplicate reads
    UMI_TOOLS_DEDUP(
        TAG_FEATURES.out.bam,
        ch_images_pulled
    )
    ch_dedup_bam = UMI_TOOLS_DEDUP.out.bam
    ch_dedup_count = UMI_TOOLS_DEDUP.out.count

    // Tag cell barcodes
    TAG_MERGED_BARCODES(
        ch_dedup_bam
            .join(ch_barcode_translate),
        ch_images_pulled
    )

    ch_final_bam = TAG_MERGED_BARCODES.out.bam

    // Generate bead counts pre-dedup
    BEAD_COUNTS(
        TAG_FEATURES.out.bam,
        ch_images_pulled
    )

    // Generate bead counts for unmapped reads
    BEAD_COUNTS_UNMAPPED(
        TAG_UNMAPPED_BARCODES.out.bam,
        ch_images_pulled
    )

    // Combine bead counts files
    COMBINE_READ_COUNTS(
        BEAD_COUNTS.out.sequence_counts
            .join(BEAD_COUNTS_UNMAPPED.out.sequence_counts),
        ch_images_pulled
    )

    // Summarize read counts per feature per cell
    UMI_TOOLS_COUNT(
        ch_final_bam,
        ch_images_pulled
    )

    // Output feature count matrix
    MAKE_COUNT_MATRIX(
        UMI_TOOLS_COUNT.out.count
            .combine(ch_reference_symbols),
        ch_images_pulled
    )

    // Identify putative cells
    CELL_CALLING(
        MAKE_COUNT_MATRIX.out.matrix,
        ch_images_pulled
    )
    ch_umi_counts = UMI_TOOLS_COUNT.out.count
    ch_allowlist = CELL_CALLING.out.allowlist

    // Convert allowlist to bead barcodes
    UNMERGE_ALLOWLIST(
        ch_allowlist
            .join(ch_barcode_translate),
        ch_images_pulled
    )

    // Edge metadata collection
    EDGE_METADATA(
        ch_edges_input
            .join(ch_corrected_edges)
            .join(UNMERGE_ALLOWLIST.out.allowlist),
        ch_images_pulled
    )

    // Mark duplicate reads
    ID_DUPLICATES(
        TAG_FEATURES.out.bam
            .join(UNMERGE_ALLOWLIST.out.allowlist),
        ch_images_pulled
    )

    // Convert duplicate counts to drop barcodes
    MERGE_DUPLICATES(
        ID_DUPLICATES.out.duplicate_count
            .join(ch_barcode_translate),
        ch_images_pulled
    )

    ch_id_dup_out = MERGE_DUPLICATES.out.duplicate_count

    // Output feature count matrix with only putative cells
    FILTER_COUNT_MATRIX(
        ch_umi_counts
            .join(ch_allowlist)
            .combine(ch_reference_symbols),
        ch_images_pulled
    )

    // Calculate summary stats
    SUMMARIZE_EXPRESSION(
        ch_final_bam
            .combine(ch_reference_interval_list),
        ch_images_pulled
    )

    // Get counts per bead if it's a bead merging experiment
    SUMMARIZE_BEAD_EXPRESSION(
        ch_final_bam
        .combine(ch_reference_interval_list),
        ch_images_pulled
    )

    // Aggregate bead stats into a table
    AGGREGATE_BEAD_STATS(
        COMBINE_READ_COUNTS.out.sequence_counts
            .join(ch_allowlist)
            .join(SUMMARIZE_EXPRESSION.out.count)
            .join(SUMMARIZE_BEAD_EXPRESSION.out.count)
            .join(ch_barcode_translate),
        ch_images_pulled
    )

    // Calculate summary stats for mixed species experiments
    ch_mixed_metrics = Channel.empty()
    ch_mixed_report = Channel.empty()
    ch_mixed_features = Channel.empty()
    ch_crosstalk_density = Channel.empty()
    if (params.rna.mixed) {
        ch_species = ch_reference_saf
            .splitText()
            .map { it.split('\\t')[1] }
            .map { it.split('\\.')[0] }
            .unique()
            .collect()

        SUMMARIZE_MIXED_SPECIES(
            ch_species,
            SUMMARIZE_EXPRESSION.out.count
                .join(ch_allowlist),
            ch_images_pulled
        )

        ch_mixed_expression = SUMMARIZE_MIXED_SPECIES.out.s1_allowlist
            .mix(SUMMARIZE_MIXED_SPECIES.out.s2_allowlist)
            .combine(ch_mixed_symbols.flatten().map { item -> tuple(item.name - '_gene_symbols.txt', item) }, by:0)
            .map { i -> [ i[1], i[0], i[2], i[3] ] }
            .combine(ch_umi_counts, by:0)

        MIXED_EXPRESSION(
            ch_mixed_expression,
            ch_images_pulled
        )

        ch_mixed_metrics = SUMMARIZE_MIXED_SPECIES.out.metrics
        ch_mixed_report = SUMMARIZE_MIXED_SPECIES.out.count
        ch_mixed_features = MIXED_EXPRESSION.out.mixed_features
        ch_crosstalk_density = SUMMARIZE_MIXED_SPECIES.out.crosstalk_density

        // Split the BAM files on a per species basis
        SPLIT_MIXED_BAM(
            POST_ALIGN_PROCESSING.out.bam
                .join(POST_ALIGN_PROCESSING.out.count),
            ch_species,
            ch_images_pulled
        )

        ch_new_bam = SPLIT_MIXED_BAM.out.bam1
            .mix(SPLIT_MIXED_BAM.out.bam2)

        // Calculate alignment metrics
        PICARD_MIXED_SPECIES(
            ch_new_bam
            .combine(ch_reference_refflat)
            .combine(ch_reference_interval_list),
            ch_images_pulled
        )

        ch_split_mixed_bam_stats = SPLIT_MIXED_BAM.out.count
        ch_picard_mixed_species = PICARD_MIXED_SPECIES.out.metrics
    } else {
        ch_split_mixed_bam_stats = Channel.empty()
        ch_picard_mixed_species = Channel.empty()
    }

    // Normalize and cluster data
    SEURAT(
        FILTER_COUNT_MATRIX.out.matrix,
        ch_images_pulled
    )

    // Aggregate data into an H5AD format
    PACK_SINGLE_H5AD(
        FILTER_COUNT_MATRIX.out.matrix
            .join(SUMMARIZE_EXPRESSION.out.count)
            .join(SEURAT.out.metadata),
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

    // Calculate alignment metrics
    PICARD(
        POST_ALIGN_PROCESSING.out.bam.map { i -> i.take(2) }
            .combine(ch_reference_refflat)
            .combine(ch_reference_interval_list),
        ch_images_pulled
    )

    // Summarize feature counts per sample
    GENES_PER_SAMPLE(
        MAKE_COUNT_MATRIX.out.matrix.map { i -> i.drop(1) } .collect(),
        ch_images_pulled
    )

    // Calculate summary matrix characteristics
    COUNT_MATRIX_FEATURES(
        FILTER_COUNT_MATRIX.out.matrix,
        ch_images_pulled
    )

    // Calculate bead distribution plots
    CALCULATE_BEAD_PLOTS(
        AGGREGATE_BEAD_STATS.out.bead_summary,
        ch_images_pulled
    )

    // Combine useful metrics into a single file
    AGGREGATE_METRICS(
        Channel.empty()
            .mix(ch_fastqc_zip,
                ch_fastqc_count,
                ch_cutadapt_metrics,
                ch_mixed_metrics,
                STAR_ALIGN.out.log,
                PICARD.out.metrics,
                ch_picard_mixed_species,
                SUMMARIZE_EXPRESSION.out.count,
                AGGREGATE_BEAD_STATS.out.metrics,
                CELL_CALLING.out.all_but_allowlist,
                ch_allowlist,
                COUNT_MATRIX_FEATURES.out.count,
                ch_debarcoded_count,
                POST_ALIGN_PROCESSING.out.count,
                ch_split_mixed_bam_stats,
                ch_dedup_count,
                EDGE_METADATA.out.metadata,
                ch_mixed_features,
                CALCULATE_BEAD_PLOTS.out.lambda_summary)
            .flatMap { i -> i[1] }
            .collect(),
        assay,
        ch_images_pulled
    )

    ch_seq_sat = Channel.empty()
    if (params.rna.mixed) {
        ch_mixed_seq_sat = ch_mixed_expression
            .map { i -> [ i[0], i[1], i[2], i[3] ] }
            .combine(ch_id_dup_out, by:0)
            .combine(AGGREGATE_METRICS.out.summary)

        // Calculate sequence saturation curves for mixed species
        SEQUENCE_SATURATION_MIXED(
            ch_mixed_seq_sat,
            assay,
            ch_images_pulled
        )

        ch_seq_sat = SEQUENCE_SATURATION_MIXED.out.results
    } else {
        // Calculate sequence saturation curves
        SEQUENCE_SATURATION(
            ch_id_dup_out
            .combine(AGGREGATE_METRICS.out.summary),
            assay,
            ch_images_pulled
        )

        ch_seq_sat = SEQUENCE_SATURATION.out.results
    }

    // Render customer facing analysis report
    RENDER_REPORT(
        ch_samples,
        ch_settings,
        Channel.empty()
            .mix(ch_input_params_yaml,
                ch_params_yaml,
                ch_command_txt,
                ch_fastqc_sequence_traces,
                ch_fastqc_quality_scores,
                SEURAT.out.umap_csv,
                SEURAT.out.top_features_csv,
                ch_seq_sat,
                CELL_CALLING.out.results,
                AGGREGATE_METRICS.out.summary,
                CALCULATE_BEAD_PLOTS.out.distributions,
                ch_crosstalk_density,
                ch_messages,
                CELL_CALLING.out.messages,
                SEURAT.out.messages)
            .flatMap { it.getClass() == ArrayList ? it[1] : it }
            .collect(),
        ch_images_pulled
    )

    emit:
    unfiltered_matrix    = MAKE_COUNT_MATRIX.out.matrix
    // tuple: [ sampleId, [ *.unfiltered.mtx.gz, *.unfiltered.barcodes.tsv, *.unfiltered.genes.tsv ] ]
    filtered_matrix      = FILTER_COUNT_MATRIX.out.matrix
    // tuple: [ sampleId, [ *.filtered.mtx.gz, *.filtered.barcodes.tsv, *.filtered.genes.tsv ] ]
    batch_h5ad           = PACK_BATCH_H5AD.out.h5ad // path: all_samples.h5ad
    cell_calling_results = CELL_CALLING.out.results
    // tuple: [ sampleId, *.csv ]
    allowlist            = ch_allowlist
    // tuple: [ sampleId, ${sampleId}_barcode_allowlist.csv ]
    translate            = ch_barcode_translate
    // tuple: [ sampleId, *_barcodeTranslate.tsv ]
}
