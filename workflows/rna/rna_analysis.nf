/*
Analyze single-cell 3' RNA Droplet data
*/

include { PUBLISH_PARAMETERS      } from '../../modules/core/publish_parameters.nf' \
          addParams(options: params.rna)
include { VALIDATE_FASTQS         } from '../../modules/core/validate_fastqs.nf' \
          addParams(options: params.rna)
include { FASTQC                  } from '../../modules/core/fastqc.nf'           addParams(options: params.rna)
include { CUTADAPT_TRIM           } from '../../modules/rna/cutadapt_trim.nf'     addParams(options: params.rna)
include { DEBARCODER              } from '../../modules/rna/debarcoder.nf'  addParams(options: params.rna)
include { CORRECT_EDGES           } from '../../modules/rna/correct_edges.nf'     addParams(options: params.rna)
include { MERGE_BEADS             } from '../../modules/rna/merge_beads.nf'       addParams(options: params.rna)
include { CORRECT_CELL_BARCODES   } from '../../modules/rna/correct_cell_barcodes.nf' \
            addParams(options: params.rna)
include { EDGE_METADATA           } from '../../modules/rna/edge_metadata.nf'     addParams(options: params.rna)
include { STAR_ALIGN              } from '../../modules/rna/star_align.nf'        addParams(options: params.rna)
include { POST_ALIGN_PROCESSING   } from '../../modules/rna/post_align_processing.nf' \
            addParams(options: params.rna)
include { TAG_BARCODES            } from '../../modules/rna/tag_barcodes.nf'      addParams(options: params.rna)
include { TAG_BARCODES as TAG_UNMAPPED_BARCODES } from '../../modules/rna/tag_barcodes.nf' \
          addParams(options: params.rna)
include { BEAD_COUNTS as BEAD_COUNTS_UNMAPPED } from '../../modules/rna/bead_counts.nf' \
          addParams(options: params.rna)
include { TAG_MERGED_BARCODES     } from '../../modules/rna/tag_merged_barcodes.nf' \
          addParams(options: params.rna)
include { TAG_GENES               } from '../../modules/rna/tag_genes.nf'         addParams(options: params.rna)
include { TAG_FEATURES            } from '../../modules/rna/tag_features.nf'      addParams(options: params.rna)
include { BEAD_COUNTS             } from '../../modules/rna/bead_counts.nf'       addParams(options: params.rna)
include { COMBINE_BEAD_COUNTS     } from '../../modules/rna/combine_read_counts.nf' \
          addParams(options: params.rna)
include { UMI_TOOLS_DEDUP         } from '../../modules/rna/umi_tools_dedup.nf'   addParams(options: params.rna)
include { UMI_TOOLS_COUNT         } from '../../modules/rna/umi_tools_count.nf'   addParams(options: params.rna)
include { UNMERGE_ALLOWLIST       } from '../../modules/rna/unmerge_allowlist.nf' addParams(options: params.rna)
include { ID_DUPLICATES           } from '../../modules/rna/id_duplicates.nf'     addParams(options: params.rna)
include { MERGE_DUPLICATES        } from '../../modules/rna/merge_duplicates.nf'  addParams(options: params.rna)
include { MAKE_COUNT_MATRIX       } from '../../modules/rna/make_count_matrix.nf' \
          addParams(options: params.rna)
include { CELL_CALLING            } from '../../modules/rna/cell_calling.nf'      addParams(options: params.rna)
include { FILTER_COUNT_MATRIX     } from '../../modules/rna/filter_count_matrix.nf' \
          addParams(options: params.rna)
include { SUMMARIZE_EXPRESSION    } from '../../modules/rna/summarize_expression.nf' \
          addParams(options: params.rna)
include { SPLIT_MIXED_BAM         } from '../../modules/rna/split_mixed_bam.nf'   addParams(options: params.rna)
include { SUMMARIZE_BEAD_EXPRESSION } from '../../modules/rna/summarize_bead_expression.nf' \
          addParams(options: params.rna)
include { AGGREGATE_BEAD_STATS    } from '../../modules/rna/aggregate_bead_stats.nf' \
          addParams(options: params.rna)
include { SUMMARIZE_MIXED_SPECIES } from '../../modules/rna/summarize_mixed_species.nf' \
          addParams(options: params.rna)
include { MIXED_EXPRESSION        } from '../../modules/rna/mixed_expression.nf'  addParams(options: params.rna)
include { SEURAT                  } from '../../modules/rna/seurat.nf'            addParams(options: params.rna)
include { PACK_SINGLE_H5AD        } from '../../modules/rna/pack_single_h5ad.nf'  addParams(options: params.rna)
include { PACK_BATCH_H5AD         } from '../../modules/rna/pack_batch_h5ad.nf'   addParams(options: params.rna)
include {
    PICARD as PICARD
    PICARD as PICARD_MIXED_SPECIES } from '../../modules/rna/picard.nf'           addParams(options: params.rna)
include { GENES_PER_SAMPLE        } from '../../modules/rna/genes_per_sample.nf'  addParams(options: params.rna)
include { COUNT_MATRIX_FEATURES   } from '../../modules/rna/count_matrix_features.nf' \
          addParams(options: params.rna)
include { SEQUENCE_SATURATION     } from '../../modules/core/sequence_saturation.nf' \
          addParams(options: params.rna)
include { SEQUENCE_SATURATION_MIXED } from '../../modules/rna/sequence_saturation_mixed.nf' \
          addParams(options: params.rna)
include { AGGREGATE_METRICS       } from '../../modules/rna/aggregate_metrics.nf' addParams(options: params.rna)
include { RENDER_REPORT  } from '../../modules/rna/render_report.nf' \
          addParams(options: params.rna)
include { CALCULATE_BEAD_PLOTS    } from '../../modules/rna/calculate_bead_plots.nf' \
          addParams(options: params.rna)

workflow RNA_ANALYSIS {
    take:
    ch_fastq // tuple: [ sampleId, raw R1 FASTQ files, raw R2 FASTQ files ]
    ch_params_file // path: input params file
    ch_reference_index // path: prepared reference index directory
    ch_reference_saf // path: prepared reference SAF file
    ch_reference_symbols // path: prepared reference gene symbol file
    ch_reference_refflat // path: prepared reference refFlat file
    ch_reference_interval_list // path: prepared ribosomal interval list file
    ch_mixed_symbols // paths: prepared reference gene symbol files for each species
    assay // string: the name of the assay
    ch_messages // any warning or info messages created by the pipeline
    ch_images_pulled // boolean: true if Singularity pull completed

    main:
    ch_fastq_complete = ch_fastq
        .map { out -> out.flatten() }
        .map { out -> tuple(out[0], [out[1], out[2]]) }

    // Validate the fastq files
    VALIDATE_FASTQS(
        ch_fastq_complete
            .map { out -> out.flatten() }
            .map { out -> tuple(out[0], [out[1]], [out[2]]) },
        ch_images_pulled
    )

    // Sort error files into passed and failed
    compile_fastq_validations = VALIDATE_FASTQS.out.error_files
        .collect()
        .flatten()
        .branch {
            passed: it.size() == 0
            failed: it.size() != 0
        }

    // Halt the pipeline if any samples failed
    compile_fastq_validations.failed
        .flatMap { it }
        .collectFile(name:'fastq_validation_errors.txt', sort: { it.name },
            newLine: false, storeDir: "${params.rna.reportsDir}")
        .subscribe {
            log.error("ERROR: [$params.rna.assay] FASTQ read files cannot be validated. " \
                    + "Please see ${params.rna.reportsDir}fastq_validation_errors.txt " \
                    + "for more details.")
            exit(1)
        }

    // Pass the empty files if validations passed
    ch_val_pass = compile_fastq_validations.passed
        .map { item -> tuple(item.name - "_errors.txt", item) }

    // Publish parameters
    PUBLISH_PARAMETERS(
        ch_params_file.ifEmpty { [] },
        ch_images_pulled
    )

    // Create file with sample IDs
    ch_samples = ch_fastq_complete
        .map { i -> i[0] }
        .collectFile(name: 'samples.txt', newLine: true, cache: 'lenient')

    // Create file with sample settings
    ch_settings = Channel.from(Core.getSampleSettings(params.rna.bead),
        Core.getSampleSettings(params.rna.barcode))
        .flatMap()
        .flatMap()
        .map { set -> String.join(',', set) }
        .collectFile(name: 'settings.txt', newLine: true, cache: 'lenient')

    ch_complete_fastq = ch_fastq_complete

    // Run fastqc analysis
    FASTQC(
        ch_complete_fastq
            .join(ch_val_pass),
        ch_images_pulled
    )

    // Trim and remove low-quality sequences
    CUTADAPT_TRIM(
        ch_complete_fastq
            .join(ch_val_pass),
        ch_images_pulled
    )

    ch_cutadapt_metrics = CUTADAPT_TRIM.out.readlength
        .mix(CUTADAPT_TRIM.out.count)
    ch_cutadapt_report = CUTADAPT_TRIM.out.log

    // Parse and correct capture and deconvolution oligos
    DEBARCODER(
        CUTADAPT_TRIM.out.fastq,
        ch_images_pulled
    )

    ch_parsed_fastq = DEBARCODER.out.fastq

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

    ch_edges_input = DEBARCODER.out.edges

    // Correct DOs and UMIs
    CORRECT_EDGES(
        ch_edges_input,
        ch_images_pulled
    )

    // Merge beads
    MERGE_BEADS(
        CORRECT_EDGES.out.edges
            .join(TAG_BARCODES.out.beads),
        ch_images_pulled
    )

    // Correct barcode formatting
    CORRECT_CELL_BARCODES(
        MERGE_BEADS.out.barcode_translate
            .join(MERGE_BEADS.out.filtered_beads),
        ch_images_pulled
    )

    // Tag cell barcodes
    TAG_MERGED_BARCODES(
        ch_dedup_bam
            .join(CORRECT_CELL_BARCODES.out.barcode_translate),
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
    COMBINE_BEAD_COUNTS(
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
            .join(CORRECT_CELL_BARCODES.out.barcode_translate),
        ch_images_pulled
    )

    // Edge metadata collection
    EDGE_METADATA(
        ch_edges_input
            .join(CORRECT_EDGES.out.edges)
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
            .join(CORRECT_CELL_BARCODES.out.barcode_translate),
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
        COMBINE_BEAD_COUNTS.out.sequence_counts
            .join(ch_allowlist)
            .join(SUMMARIZE_EXPRESSION.out.count)
            .join(SUMMARIZE_BEAD_EXPRESSION.out.count)
            .join(CORRECT_CELL_BARCODES.out.barcode_translate),
        ch_images_pulled
    )

    // Calculate summary stats for mixed species experiments
    ch_mixed_stats = Channel.empty()
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

        ch_mixed_stats = SUMMARIZE_MIXED_SPECIES.out.count
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
            .mix(FASTQC.out.zip,
                FASTQC.out.counts,
                ch_cutadapt_metrics,
                ch_mixed_stats,
                STAR_ALIGN.out.log,
                PICARD.out.metrics,
                ch_picard_mixed_species,
                SUMMARIZE_EXPRESSION.out.count,
                CELL_CALLING.out.all_but_allowlist,
                ch_allowlist,
                COUNT_MATRIX_FEATURES.out.count,
                DEBARCODER.out.count,
                POST_ALIGN_PROCESSING.out.count,
                ch_split_mixed_bam_stats,
                ch_dedup_count,
                EDGE_METADATA.out.metadata,
                ch_mixed_features)
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

    // Render QC report
    ch_merge_beads_report = CALCULATE_BEAD_PLOTS.out.distributions

    // Render customer facing analysis report
    RENDER_REPORT(
        ch_samples,
        ch_settings,
        Channel.empty()
            .mix(PUBLISH_PARAMETERS.out.input_params_yaml,
                PUBLISH_PARAMETERS.out.params_yaml,
                PUBLISH_PARAMETERS.out.commmand_txt,
                FASTQC.out.sequence_traces,
                FASTQC.out.quality_scores,
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
    unfiltered_matrix = MAKE_COUNT_MATRIX.out.matrix
    // tuple: [ sampleId, [ *.unfiltered.mtx.gz, *.unfiltered.barcodes.tsv, *.unfiltered.genes.tsv ] ]
    filtered_matrix   = FILTER_COUNT_MATRIX.out.matrix
    // tuple: [ sampleId, [ *.filtered.mtx.gz, *.filtered.barcodes.tsv, *.filtered.genes.tsv ] ]
    batch_h5ad        = PACK_BATCH_H5AD.out.h5ad // path: all_samples.h5ad
}
