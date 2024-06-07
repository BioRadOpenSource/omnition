/*
Analyze single-cell ATAC data
*/

include { PUBLISH_PARAMETERS             } from '../../modules/core/publish_parameters.nf'            \
    addParams(options: params.atac)
include { VALIDATE_FASTQS                } from '../../modules/core/validate_fastqs.nf'               \
    addParams(options: params.atac)
include { MERGE_LANES                    } from '../../modules/atac/merge_lanes.nf'                   \
    addParams(options: params.atac)
include { CHECK_MITO_CONTIG              } from '../../modules/atac/check_mito_contig'                \
    addParams(options: params.atac)
include { FASTQC                         } from '../../modules/core/fastqc.nf'                        \
    addParams(options: params.atac)
include { DEBARCODER                     } from '../../modules/atac/debarcoder.nf'                    \
    addParams(options: params.atac)
include { CUTADAPT_HEADCROP              } from '../../modules/atac/cutadapt_headcrop.nf'             \
    addParams(options: params.atac)
include { BWA_ALIGNMENT                  } from '../../modules/atac/bwa_alignment.nf'                 \
    addParams(options: params.atac)
include { MARK_DUPLICATES                } from '../../modules/atac/mark_duplicates.nf'               \
    addParams(options: params.atac)
include { SPLIT_BAM                      } from '../../modules/atac/split_bam.nf'                     \
    addParams(options: params.atac)
include { ASSEMBLE_FRAGMENTS             } from '../../modules/atac/assemble_fragments.nf'            \
    addParams(options: params.atac)
include { ANNOTATE_FRAGMENTS             } from '../../modules/atac/annotate_fragments.nf'            \
    addParams(options: params.atac)
include { DETERMINE_BARCODE_ALLOWLIST    } from '../../modules/atac/determine_barcode_allowlist.nf'   \
    addParams(options: params.atac)
include { COMPUTE_DECONVOLUTION_STAT_CHR } from '../../modules/atac/compute_deconvolution_stat_chr.nf'\
    addParams(options: params.atac)
include { DETERMINE_BARCODE_MERGES       } from '../../modules/atac/determine_barcode_merges.nf'      \
    addParams(options: params.atac)
include { REANNOTATE_FRAGMENTS           } from '../../modules/atac/reannotate_fragments.nf'          \
    addParams(options: params.atac)
include { REANNOTATE_BAM                 } from '../../modules/atac/reannotate_bam.nf'                \
    addParams(options: params.atac)
include { FINAL_BAM_MERGE                } from '../../modules/atac/final_bam_merge.nf'               \
    addParams(options: params.atac)
include { MERGE_REANN_READ_COUNTS        } from '../../modules/atac/merge_reann_read_counts.nf'       \
    addParams(options: params.atac)
include { FINAL_FRAG_MERGE               } from '../../modules/atac/final_frag_merge.nf'              \
    addParams(options: params.atac)
include { ASSEMBLE_BASIC_QC              } from '../../modules/atac/assemble_basic_qc.nf'             \
    addParams(options: params.atac)
include { FINAL_QC_SE                    } from '../../modules/atac/final_qc_se.nf'                   \
    addParams(options: params.atac)
include { CALCULATE_BEADS_PER_DROP       } from '../../modules/atac/calc_beads_per_drop.nf'           \
    addParams(options: params.atac)
include { SUMMARIZE_ALIGNMENTS           } from '../../modules/atac/summarize_alignments.nf'          \
    addParams(options: params.atac)
include { COMPUTE_TSS_MATRIX             } from '../../modules/atac/compute_tss_matrix.nf'            \
    addParams(options: params.atac)
include { CALL_PEAKS                     } from '../../modules/atac/call_peaks.nf'                    \
    addParams(options: params.atac)
include { CLEAN_PEAKS                    } from '../../modules/atac/clean_peaks.nf'                   \
    addParams(options: params.atac)
include { MAKE_COUNT_MATRIX              } from '../../modules/atac/make_count_matrix.nf'         \
    addParams(options: params.atac)
include { FRACTION_OF_READS_IN_PEAKS     } from '../../modules/atac/fraction_of_reads_in_peaks.nf'    \
    addParams(options: params.atac)
include { FRACTION_OF_READS_IN_TSS       } from '../../modules/atac/fraction_of_reads_in_tss.nf'      \
    addParams(options: params.atac)
include { CALCULATE_INSERT_SIZE_METRICS  } from '../../modules/atac/calculate_insert_size_metrics.nf' \
    addParams(options: params.atac)
include { SUMMARIZE_MIXED_SPECIES        } from '../../modules/atac/summarize_mixed_species.nf'       \
    addParams(options: params.atac)
include { SEQUENCE_SATURATION            } from '../../modules/core/sequence_saturation.nf'           \
    addParams(options: params.atac)
include { ARCHR                          } from '../../modules/atac/archr.nf'                         \
    addParams(options: params.atac)
include { TSS_ENRICHMENT                 } from '../../modules/atac/tss_enrichment.nf'                \
    addParams(options: params.atac)
include { AGGREGATE_METRICS              } from '../../modules/atac/aggregate_metrics.nf'             \
    addParams(options: params.atac)
include { BUILD_REPORT_CONTENTS          } from '../../modules/atac/build_report_contents.nf'         \
    addParams(options: params.atac)
include { RENDER_REPORT                  } from '../../modules/atac/render_report.nf'                 \
    addParams(options: params.atac)

workflow ATAC_NORMAL_ANALYSIS {
    take:
    ch_fastq // tuple: [ sampleId, raw R1 FASTQ files, raw R2 FASTQ files ]
    ch_params_file // path: input params file
    ch_reference_fasta // Prepared reference FASTA file
    ch_gtf // Prepared reference GTF file
    ch_index // Prepared indexed references
    ch_sizes // List of chromosome sizes
    ch_archr_ref // tar.gz of BSgenome package for reference fasta
    ch_blocklist // Prepared blocklist
    ch_tss // BED file of TSS
    ch_messages // any warning or info messages created by the pipeline
    assay // string: name of the assay, used in core functions
    ch_images_pulled // boolean: true if Singularity pull completed

    main:
    // Publish parameters
    PUBLISH_PARAMETERS(
        ch_params_file.ifEmpty { [] },
        ch_images_pulled
    )

    // Merge lanesplit fastqs
    MERGE_LANES(
        ch_fastq,
        ch_images_pulled
    )

    // Validate the fastq files
    VALIDATE_FASTQS(
        MERGE_LANES.out.complete
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
            newLine: false, storeDir: "${params.atac.reportsDir}")
        .subscribe {
            log.error("ERROR: [$params.atac.assay] FASTQ read files cannot be validated. " \
                    + "Please see ${params.atac.reportsDir}fastq_validation_errors.txt " \
                    + "for more details.")
            exit(1)
        }

    // Pass the empty files if validations passed
    ch_val_pass = compile_fastq_validations.passed
        .map { item -> tuple(item.name - "_errors.txt", item) }

    // Check mitochondrial contig
    CHECK_MITO_CONTIG(
        ch_sizes,
        ch_images_pulled
    )

    // Run fastqc analysis
    FASTQC(
        MERGE_LANES.out.complete
            .join(ch_val_pass),
        ch_images_pulled
    )

    // Parse and correct capture oligos
    DEBARCODER(
        MERGE_LANES.out.complete
            .map { out -> tuple(out[0], out[1], []) }
            .join(ch_val_pass),
        ch_images_pulled
    )

    ch_zero = 0

    ch_sample_map = DEBARCODER.out.fastq
        .map { item  -> "${item[0]},${item[0]}" }
        .collectFile(name:'sample_map.csv', seed: 'sample,fastq', newLine: true, // groovylint-disable-line
            sort: true, cache: 'lenient', storeDir: "${params.atac.reportsDir}")
        .first()

    ch_complete_fastq = DEBARCODER.out.fastq

    // Remove bases from beginning of R2 reads
    CUTADAPT_HEADCROP(
        ch_complete_fastq,
        ch_images_pulled
    )

    // Align reads to refernce genome then tag with cell barcode
    BWA_ALIGNMENT(
        CUTADAPT_HEADCROP.out.fastq
            .combine(ch_index),
        ch_images_pulled
    )

    // Mark duplicated reads
    MARK_DUPLICATES(
        BWA_ALIGNMENT.out.bam,
        ch_images_pulled
    )

    // Split bam files by chromosome
    SPLIT_BAM(
        MARK_DUPLICATES.out.bam,
        ch_sizes,
        ch_images_pulled
    )

    // Format chromosome-specific bam channels for downstream processes
    ch_split_bam_formatted = SPLIT_BAM.out.split_bam
        .transpose()
        .map { item -> tuple(item[0], item[1].name - "${item[0]}." - '.raw.bam', item[1], item[2], item[3]) }

    // Create chromosome-specific fragment files
    ASSEMBLE_FRAGMENTS(
        ch_split_bam_formatted,
        ch_images_pulled
    )

    // Annotate chromosome-specific fragment files
    ANNOTATE_FRAGMENTS(
        ch_split_bam_formatted
            .map { tuple(it[0], it[1], it[4]) }
            .join(ASSEMBLE_FRAGMENTS.out.assemble_fragments, by: [0, 1])
            .combine(ch_blocklist),
        ch_images_pulled
    )

    // Bead filtration with knee-calling
    DETERMINE_BARCODE_ALLOWLIST(
        ANNOTATE_FRAGMENTS.out.bead_counts
            .groupTuple(),
        ch_images_pulled
    )

    // Calculate chromosome-specific stats
    COMPUTE_DECONVOLUTION_STAT_CHR(
        ANNOTATE_FRAGMENTS.out.annotate_fragments
            .combine(DETERMINE_BARCODE_ALLOWLIST.out.allowlist, by: 0),
        ch_images_pulled
    )

    ch_in_determine_merges = DETERMINE_BARCODE_ALLOWLIST.out.quant
        .join(DETERMINE_BARCODE_ALLOWLIST.out.allowlist, by: 0)
        .join(DETERMINE_BARCODE_ALLOWLIST.out.params)
        .join(COMPUTE_DECONVOLUTION_STAT_CHR.out.overlap_count
            .filter { i -> !(i[1] =~ /${params.atac.mitoContig}/) }
            .map { i  -> [ i[0], i[2] ] }
            .groupTuple(), by: 0)

    // Determining which barcodes to merge
    DETERMINE_BARCODE_MERGES(
        ch_in_determine_merges,
        ch_zero,
        ch_images_pulled
    )

    // parse the paths to get out the samplId, will be fastq for no TI and fastqti for with TI
    ch_determine_barcode_merges_implicated_barcodes = DETERMINE_BARCODE_MERGES.out.implicated_barcodes
        .flatten().map {
            i  -> [ i.toString().substring(i.toString().lastIndexOf('/') + 1,
            i.toString().indexOf('.implicatedBarcodes.csv.gz')), i ] }

    ch_determine_barcode_merges_barcode_translate = DETERMINE_BARCODE_MERGES.out.barcode_translate
        .flatten().map {
            i  -> [ i.toString().substring(i.toString().lastIndexOf('/') + 1, // groovylint-disable-line
            i.toString().indexOf('.barcodeTranslate.tsv')), i ] }

    ch_determine_barcode_merges_params = DETERMINE_BARCODE_MERGES.out.params
        .flatten().map {
            i  -> [ i.toString().substring(i.toString().lastIndexOf('/') + 1, // groovylint-disable-line
            i.toString().indexOf('.deconvolutionParams.csv')), i ] }

    // Reannotate fragments files based on merging scheme
    REANNOTATE_FRAGMENTS(
        ANNOTATE_FRAGMENTS.out.annotate_fragments
            .combine(ch_determine_barcode_merges_barcode_translate, by: 0),
        ch_images_pulled
    )

    // Reannotate bam files based on merging scheme
    REANNOTATE_BAM(
        ch_split_bam_formatted
            .map { tuple(it[0], it[1], it[2], it[3]) }
            .join(REANNOTATE_FRAGMENTS.out.chr_reannotate_fragments, by: [0, 1])
            .combine(ch_determine_barcode_merges_barcode_translate, by: 0),
        ch_images_pulled
    )

    // Merge chromosome bam files back together
    FINAL_BAM_MERGE(
        REANNOTATE_BAM.out.reannotate_bam
            .filter { i -> !(i[1] =~ /${params.atac.mitoContig}/) }
            .map { tuple(it[0], it[2], it[3]) }
            .groupTuple(),
        ch_images_pulled
    )

    MERGE_REANN_READ_COUNTS(
        REANNOTATE_BAM.out.count
            .map { tuple(it[0], it[1]) }
            .groupTuple(),
        ch_images_pulled
    )

    // Merge chromosome fragment files back together
    FINAL_FRAG_MERGE(
        REANNOTATE_FRAGMENTS.out.reannotate_fragments
            .filter { i -> !(i[1].name =~ /${params.atac.mitoContig}/) }
            .groupTuple(),
        ch_images_pulled
    )

    // Calculate basic QC stats
    ASSEMBLE_BASIC_QC(
        REANNOTATE_FRAGMENTS.out.frag_sumstats
        .groupTuple(),
        ch_images_pulled
    )

    // Calculate final qc Stats
    FINAL_QC_SE(
        FINAL_FRAG_MERGE.out.final_frag_merge
            .join(ASSEMBLE_BASIC_QC.out.assemble_basic_qc)
            .join(ch_determine_barcode_merges_barcode_translate)
            .combine(ch_tss),
        ch_images_pulled
    )

    // Creates csv for cells and marks as true or false
    CALCULATE_BEADS_PER_DROP(
        FINAL_QC_SE.out.qc_stats,
        ch_images_pulled
    )

    // Collect alignment stats
    SUMMARIZE_ALIGNMENTS(
        MARK_DUPLICATES.out.bam.map { tuple(it[0], it[1], it[2]) },
        ch_images_pulled
    )

    // Compute TSS matrix
    COMPUTE_TSS_MATRIX(
        FINAL_BAM_MERGE.out.final_bam
            .combine(ch_tss),
        ch_images_pulled
    )

    // Call peaks of fragment coverage
    CALL_PEAKS(
        FINAL_BAM_MERGE.out.final_bam,
        ch_images_pulled
    )

    // Clean peaks to set width and remove low probability peaks and peaks in blocklist
    CLEAN_PEAKS(
        CALL_PEAKS.out.peaks
            .combine(ch_blocklist)
            .combine(ch_sizes),
        ch_images_pulled
    )

    // Generate count matrix for called peaks
    MAKE_COUNT_MATRIX(
        FINAL_QC_SE.out.qc_stats
            .join(FINAL_BAM_MERGE.out.final_bam, by:0)
            .combine(CLEAN_PEAKS.out.peaks, by:0),
        ch_images_pulled
    )

    // Calculate fraction of reads in called peaks
    FRACTION_OF_READS_IN_PEAKS(
        FINAL_BAM_MERGE.out.final_bam
            .join(CLEAN_PEAKS.out.peaks, by:0),
        ch_images_pulled
    )

    // Calculate fraction of reads in TSS
    FRACTION_OF_READS_IN_TSS(
        FINAL_BAM_MERGE.out.final_bam
            .combine(ch_tss),
        ch_images_pulled
    )

    // Calculate insert size metrics
    CALCULATE_INSERT_SIZE_METRICS(
        FINAL_BAM_MERGE.out.final_bam,
        ch_images_pulled
    )

    ch_report_crosstalk = Channel.empty()

    if (params.atac.mixed) {
        ch_species = ch_gtf
            .splitText()
            .filter { !(it =~ /#/) }
            .map { it.split('\\.')[0] } // groovylint-disable-line
            .unique()
            .collect()

        // Calculate mixed species stats
        SUMMARIZE_MIXED_SPECIES(
            FINAL_QC_SE.out.qc_stats,
            ch_species,
            ch_images_pulled
        )

        ch_report_crosstalk = SUMMARIZE_MIXED_SPECIES.out.stats
    }

    ch_seq_sat_out = Channel.empty()
    // Calculate fragment saturation curve
    SEQUENCE_SATURATION(
        FINAL_FRAG_MERGE.out.final_frag_merge
            .map { i -> [i[0], i[1]] }
            .join(SUMMARIZE_ALIGNMENTS.out.stats
                .map { i -> [i[0], i[3]] }, by: 0),
        assay,
        ch_images_pulled
    )
    ch_seq_sat_out = SEQUENCE_SATURATION.out.results

    // Use ArchR
    ch_archr_output = Channel.empty()
    ARCHR(
        FINAL_BAM_MERGE.out.final_bam.join(CLEAN_PEAKS.out.peaks, by:0)
            .combine(ch_gtf)
            .combine(ch_archr_ref)
            .combine(ch_blocklist),
        ch_images_pulled
    )
    ch_archr_output = ARCHR.out.clusterinfo

    // Calculate TSS enrichment score for use in aggregate metrics
    TSS_ENRICHMENT(
        COMPUTE_TSS_MATRIX.out.tss_matrix,
        ch_images_pulled
    )

    // Combine useful metrics into a single file
    AGGREGATE_METRICS(
        Channel.empty()
            .mix(FASTQC.out.zip.map { i  -> i[1] },
                FASTQC.out.counts.map { i -> i[1] },
                CUTADAPT_HEADCROP.out.readlength.map { i  -> i[1] },
                SUMMARIZE_ALIGNMENTS.out.stats.map { i  -> i[1, 2, 3] },
                MARK_DUPLICATES.out.stats.map { i  -> i[1] },
                DETERMINE_BARCODE_ALLOWLIST.out.allowlist.map { i  -> i[1] },
                DETERMINE_BARCODE_ALLOWLIST.out.quant.map { i  -> i[1] },
                ch_determine_barcode_merges_params.map { i  -> i[1] },
                ch_determine_barcode_merges_implicated_barcodes.map { i  -> i[1] },
                TSS_ENRICHMENT.out.tss_enrichment.map { i  -> i[1] },
                FRACTION_OF_READS_IN_PEAKS.out.frip.map { i  -> i[1] },
                FRACTION_OF_READS_IN_TSS.out.frit.map { i  -> i[1] },
                FINAL_QC_SE.out.qc_stats.map { i  -> i[1] },
                ch_report_crosstalk.map { i  -> i[1, 2] },
                CLEAN_PEAKS.out.peaks.map { i  -> i[1] },
                FINAL_BAM_MERGE.out.count.map { i  -> i[1] },
                DEBARCODER.out.count.map { i  -> i[1] },
                BWA_ALIGNMENT.out.count.map { i  -> i[1] },
                REANNOTATE_BAM.out.count.map { i  -> i[1] },
                ANNOTATE_FRAGMENTS.out.stats.map { i  -> i[2] },
                ASSEMBLE_FRAGMENTS.out.stats.map { i  -> i[2] },
                REANNOTATE_FRAGMENTS.out.stats.map { i  -> i[2] },
                MERGE_REANN_READ_COUNTS.out.count.map { i  -> i[1] })
            .collect(),
        assay,
        ch_images_pulled
    )

    // Build contents for export to reporting step
    BUILD_REPORT_CONTENTS(
        Channel.empty()
            .mix(CALCULATE_INSERT_SIZE_METRICS.out.metrics.map { i -> i[1] }.collect() .ifEmpty { [] },
                AGGREGATE_METRICS.out.summary .collect() .ifEmpty { [] },
                )
            .collect(),
        ch_images_pulled
    )

    // Generate the report
    RENDER_REPORT(
        Channel.empty()
            .mix(PUBLISH_PARAMETERS.out.input_params_yaml,
                PUBLISH_PARAMETERS.out.params_yaml,
                PUBLISH_PARAMETERS.out.commmand_txt,
                BUILD_REPORT_CONTENTS.out.metric_summary,
                BUILD_REPORT_CONTENTS.out.pipeline_summary,
                FASTQC.out.zip.map { i  -> i[1] },
                DEBARCODER.out.count.map { i  -> i[1] },
                SUMMARIZE_ALIGNMENTS.out.stats.map { i  -> i[1, 2, 3] },
                FINAL_QC_SE.out.qc_stats.map { i  -> i[1] },
                ch_determine_barcode_merges_params.map { i  -> i[1] },
                ch_determine_barcode_merges_implicated_barcodes.map { i  -> i[1] },
                DETERMINE_BARCODE_ALLOWLIST.out.quant.map { i  -> i[1] },
                DETERMINE_BARCODE_ALLOWLIST.out.allowlist.map { i  -> i[1] },
                CALCULATE_BEADS_PER_DROP.out.cells.map { i  -> i[1] },
                TSS_ENRICHMENT.out.tss_enrichment.map { i  -> i[1] },
                CALCULATE_INSERT_SIZE_METRICS.out.metrics.map { i  -> i[1] },
                ch_seq_sat_out.map { i  -> i[1] },
                ch_report_crosstalk.map { i  -> i[1, 2] },
                ch_archr_output,
                AGGREGATE_METRICS.out.summary,
                ch_messages,
                ch_sample_map)
            .collect(),
        ch_images_pulled
    )

    emit:
    MAKE_COUNT_MATRIX.out.matrix // tuple: [ sampleId, column_names.txt, row_names.txt, read_counts.txt]
}
