/*
Analyze single-cell ATAC data
*/

include { MERGE_LANES                    } from '../../modules/core/merge_lanes.nf'                   \
    addParams(options: params.atac)
include { CHECK_MITO_CONTIG              } from '../../modules/atac/check_mito_contig'                \
    addParams(options: params.atac)
include { PUBLISH_PARAMETERS             } from '../../modules/core/publish_parameters.nf'                   \
    addParams(options: params.atac)
include { FASTQC                         } from '../../modules/core/fastqc.nf'                        \
    addParams(options: params.atac)
include { TI_DEAD_CONFIG                 } from '../../modules/atac/ti_dead_config.nf'                \
    addParams(options: params.atac)
include { DEAD                           } from '../../modules/atac/dead.nf'                          \
    addParams(options: params.atac)
include { CUTADAPT_HEADCROP              } from '../../modules/core/cutadapt_headcrop.nf'             \
    addParams(options: params.atac)
include { VALIDATE_TI_CONFIG             } from '../../modules/atac/validate_TI_config.nf'            \
    addParams(options: params.atac)
include { CHECK_TI_COUNTS                } from '../../modules/atac/check_ti_counts.nf'               \
    addParams(options: params.atac)
include { CHECK_TI_COUNTS_SUPERLOADED    } from '../../modules/atac/check_ti_counts_superloaded.nf'   \
    addParams(options: params.atac)
include { TI_WARNING_MESSAGES            } from '../../modules/atac/ti_warning_messages.nf'           \
    addParams(options: params.atac)
include { SPLIT_FASTQ                    } from '../../modules/atac/split_fastq.nf'                   \
    addParams(options: params.atac)
include { TI_ERROR_CHECK                 } from '../../modules/atac/ti_error_check.nf'                \
    addParams(options: params.atac)
include { BWA_ALIGNMENT                  } from '../../modules/atac/bwa_alignment.nf'                 \
    addParams(options: params.atac)
include { TAG_BARCODES                   } from '../../modules/atac/tag_barcodes.nf'                  \
    addParams(options: params.atac)
include { MARK_DUPLICATES                } from '../../modules/atac/mark_duplicates.nf'               \
    addParams(options: params.atac)
include { SUMMARIZE_ALIGNMENTS           } from '../../modules/atac/summarize_alignments.nf'          \
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
include { SEQUENCE_SATURATION            } from '../../modules/core/sequence_saturation.nf'           \
    addParams(options: params.atac)
include { COMPILE_ALIGNMENTS             } from '../../modules/atac/compile_alignments.nf'            \
    addParams(options: params.atac)
include { COMPILE_FRAGMENTS              } from '../../modules/atac/compile_fragments.nf'             \
    addParams(options: params.atac)
include { COMPILE_QC_STATS               } from '../../modules/atac/compile_qc_stats.nf'              \
    addParams(options: params.atac)
include { BEAD_FILT_SUMMARY              } from '../../modules/atac/bead_filt_summary.nf'             \
    addParams(options: params.atac)
include { COMPUTE_TSS_MATRIX             } from '../../modules/atac/compute_tss_matrix.nf'            \
    addParams(options: params.atac)
include { CALL_PEAKS                     } from '../../modules/atac/call_peaks.nf'                    \
    addParams(options: params.atac)
include { CLEAN_PEAKS                    } from '../../modules/atac/clean_peaks.nf'                   \
    addParams(options: params.atac)
include { FRACTION_OF_READS_IN_PEAKS     } from '../../modules/atac/fraction_of_reads_in_peaks.nf'    \
    addParams(options: params.atac)
include { FRACTION_OF_READS_IN_TSS       } from '../../modules/atac/fraction_of_reads_in_tss.nf'      \
    addParams(options: params.atac)
include { CALCULATE_INSERT_SIZE_METRICS  } from '../../modules/atac/calculate_insert_size_metrics.nf' \
    addParams(options: params.atac)
include { CALCULATE_BEADS_PER_DROP       } from '../../modules/atac/calc_beads_per_drop.nf'           \
    addParams(options: params.atac)
include { SUMMARIZE_MIXED_SPECIES        } from '../../modules/atac/summarize_mixed_species.nf'       \
    addParams(options: params.atac)
include { MAKE_COUNT_MATRIX              } from '../../modules/atac/make_count_matrix_tmp.nf'         \
    addParams(options: params.atac)
include { ARCHR                          } from '../../modules/atac/archr.nf'                         \
    addParams(options: params.atac)
include { TSS_ENRICHMENT                 } from '../../modules/atac/tss_enrichment.nf'                \
    addParams(options: params.atac)
include { AGGREGATE_METRICS              } from '../../modules/atac/aggregate_metrics.nf'             \
    addParams(options: params.atac)
include { BUILD_REPORT_CONTENTS          } from '../../modules/atac/build_report_contents.nf'         \
    addParams(options: params.atac)
include { GENERATE_REPORT                } from '../../modules/atac/generate_report.nf'               \
    addParams(options: params.atac)

workflow ATAC_ANALYSIS {
    take:
    ch_fastq // tuple: [ sample_id, raw R1 FASTQ files, raw R2 FASTQ files ]
    ch_reference_fasta // Prepared reference FASTA file
    ch_gtf // Prepared reference GTF file
    ch_index // Prepared indexed references
    ch_sizes // List of chromosome sizes
    ch_archr_ref // tar.gz of BSgenome package for reference fasta
    ch_blocklist // Prepared blocklist
    ch_tss // BED file of TSS
    ch_TI_config // CSV config file of fastq and TIs
    ch_messages // any warning or info messages created by the pipeline
    assay // string: name of the assay, used in core functions
    ch_images_pulled // boolean: true if Singularity pull completed

    main:
    // Publish parameters
    PUBLISH_PARAMETERS(
        ch_images_pulled
    )

    // Merge lanesplit fastqs
    MERGE_LANES(
        ch_fastq,
        ch_images_pulled
    )

    // Check mitochondrial contig
    CHECK_MITO_CONTIG(
        ch_sizes,
        ch_images_pulled
    )

    // Run fastqc analysis
    FASTQC(
        MERGE_LANES.out.complete,
        ch_images_pulled
    )

    ch_optional_config = Channel.empty().first()
    ch_ti_len = 0

    if (params.atac.barcodedTn5) {
        // Create CSV file of TIs from Nextflow parameter
        ti_csv_ch = Channel
            .from("${ params.atac.ti }")
            .map { it.replaceAll('\\s', '') }
            .map { it.replaceAll(',', '\n') }
            .map { it.replaceAll(':', ',') } // groovylint-disable-line
            .map { it.replaceAll('\\[', '') } // groovylint-disable-line
            .map { it.replaceAll('\\]', '') } // groovylint-disable-line
            .collectFile(name:'TI.csv', seed: 'name,sequence\n', cache: 'lenient')
            .first()

        // Create a TI config for dead
        TI_DEAD_CONFIG(
            ti_csv_ch,
            Channel.fromPath("${projectDir}/assets/atac.json").first(),
            ch_images_pulled
        )

        ch_optional_config = TI_DEAD_CONFIG.out.config
        ch_ti_len = TI_DEAD_CONFIG.out.tilen.text
    }

    // Parse and correct capture oligos
    DEAD(
        MERGE_LANES.out.complete,
        ch_optional_config.ifEmpty([]),
        ch_images_pulled
    )

    ch_fastq_ti_pairs_read_counts = Channel.empty()
    ch_validated_ti_config = Channel.empty()

    if (params.atac.barcodedTn5) {
        if (params.atac.barcodedTn5Config) {
            // Validate TI config
            VALIDATE_TI_CONFIG(
                ch_TI_config,
                ch_images_pulled,
                ti_csv_ch
            )

            ch_validated_ti_config = VALIDATE_TI_CONFIG.out.config

            // Check counts for each TI
            CHECK_TI_COUNTS(
                DEAD.out.fastq,
                ti_csv_ch,
                VALIDATE_TI_CONFIG.out.config.first(),
                ch_images_pulled
            )

            // Combine Warning messages into one output file
            TI_WARNING_MESSAGES(
                CHECK_TI_COUNTS.out.messages.collect(),
                ch_messages
            )

            // Compile all errors for each fastq file
            TI_ERROR_CHECK(
                CHECK_TI_COUNTS.out.error_files.collect(),
                ch_images_pulled
            )

            ch_fastq_ti_pairs_read_counts = CHECK_TI_COUNTS.out.split
                .flatMap { it }
                .collectFile(name:'fastqTIreadcounts.csv', seed: 'sample,fastq,ti,sequence,count',
                    newLine: true, cache: 'lenient', storeDir: "${params.reportsDir}")

            ch_fastq_ti_pairs = ch_fastq_ti_pairs_read_counts
                .splitCsv(header:true)
                .map { row -> tuple(row.fastq, row.sequence) }
                .combine(DEAD.out.fastq, by:0)

            // Split each fastq by TI
            SPLIT_FASTQ(
                ch_fastq_ti_pairs,
                TI_ERROR_CHECK.out.passed,
                ch_images_pulled
            )

            ch_prepped_fastq = SPLIT_FASTQ.out.split
        } else {
            // Check counts for each TI
            CHECK_TI_COUNTS_SUPERLOADED(
                DEAD.out.fastq,
                ti_csv_ch,
                ch_images_pulled
            )

            // Combine Warning messages into one output file
            TI_WARNING_MESSAGES(
                CHECK_TI_COUNTS_SUPERLOADED.out.messages.collect(),
                ch_messages
            )

            ch_fastq_ti_pairs_read_counts = CHECK_TI_COUNTS_SUPERLOADED.out.split
                .flatMap { it }
                .collectFile(name:'fastqTIreadcounts_superloaded.csv', seed: 'sample,fastq,ti,sequence,count', // groovylint-disable-line
                    newLine: true, cache: 'lenient', storeDir: "${params.reportsDir}")

            ch_fastq_ti_pairs = ch_fastq_ti_pairs_read_counts
                .splitCsv(header:true)
                .map { row -> tuple(row.fastq, row.sequence) }
                .combine(DEAD.out.fastq, by:0)

            // Split each fastq by TI
            SPLIT_FASTQ(
                ch_fastq_ti_pairs,
                CHECK_TI_COUNTS_SUPERLOADED.out.passed.first(),
                ch_images_pulled
            )

            ch_prepped_fastq = SPLIT_FASTQ.out.split
        }

        ch_sample_map = ch_fastq_ti_pairs_read_counts
            .splitCsv(header:true)
            .map { row -> "${row.sample},${row.fastq},${row.ti}" }
            .collectFile(name:'sample_map.csv', seed: 'sample,fastq,ti', newLine: true,
                sort: true, cache: 'lenient', storeDir: "${params.reportsDir}")
            .first()
    } else {
        ch_prepped_fastq = DEAD.out.fastq

        ch_sample_map = ch_prepped_fastq
            .map { item  -> "${item[0]},${item[0]}" }
            .collectFile(name:'sample_map.csv', seed: 'sample,fastq', newLine: true, // groovylint-disable-line
                sort: true, cache: 'lenient', storeDir: "${params.reportsDir}")
            .first()
    }

    // Remove bases from beginning of R2 reads
    CUTADAPT_HEADCROP(
        ch_prepped_fastq,
        ch_images_pulled
    )

    // Align reads to refernce genome then tag with cell barcode
    BWA_ALIGNMENT(
        CUTADAPT_HEADCROP.out.fastq,
        ch_index,
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
                      .flatMap { item  ->
                          fileTupleList = []
                          for (i = 0; i < item[1].size(); i++) {
                chr = item[1][i].name - "${item[0]}." - '.raw.bam'
                fileTuple = tuple(item[0], chr, item[1][i], item[2][i], item[3][i])
                fileTupleList.add(fileTuple)
                          }
                          return(fileTupleList)
                      }

    // Create chromosome-specific fragment files
    ASSEMBLE_FRAGMENTS(
        ch_split_bam_formatted,
        ch_images_pulled
    )

    // Annotate chromosome-specific fragment files
    ANNOTATE_FRAGMENTS(
        ch_split_bam_formatted
            .map { tuple(it[0], it[1], it[4]) }
            .join(ASSEMBLE_FRAGMENTS.out.assemble_fragments, by: [0, 1]),
        ch_blocklist,
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

    if (params.atac.barcodedTn5) {
        // if there are tis, substring removes ti and merge these different channels by fastq
        ch_in_determine_merges = DETERMINE_BARCODE_ALLOWLIST.out.quant
            .join(DETERMINE_BARCODE_ALLOWLIST.out.allowlist, by: 0)
            .join(DETERMINE_BARCODE_ALLOWLIST.out.params, by: 0)
            .map { i  -> [ i[0].substring(0, i[0].lastIndexOf('-')), i[1], i[2], i[3] ] }
            .groupTuple()
            .join(COMPUTE_DECONVOLUTION_STAT_CHR.out.overlap_count
                .filter { !(it =~ /${params.atac.mitoContig}/) }
                .map { i  -> [ i[0].substring(0, i[0].lastIndexOf('-')), i[2] ] }
                .groupTuple())
    } else {
        // if there are not tis, merge these different channels by fastq
        ch_in_determine_merges = DETERMINE_BARCODE_ALLOWLIST.out.quant
            .join(DETERMINE_BARCODE_ALLOWLIST.out.allowlist, by: 0)
            .join(DETERMINE_BARCODE_ALLOWLIST.out.params)
            .join(COMPUTE_DECONVOLUTION_STAT_CHR.out.overlap_count
                .filter { !(it =~ /${params.atac.mitoContig}/) }
                .map { i  -> [ i[0], i[2] ] }
                .groupTuple(), by: 0)
    }

    // Determining which barcodes to merge
    DETERMINE_BARCODE_MERGES(
        ch_in_determine_merges,
        ch_ti_len,
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
            .map { tuple(it[0], it[2], it[3]) }
            .filter { !(it =~ /${params.atac.mitoContig}/) }
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
            .filter { !(it =~ /${params.atac.mitoContig}/) }
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
            .join(ch_determine_barcode_merges_barcode_translate),
        ch_tss,
        ch_images_pulled
    )

    // Creates csv for cells and marks as true or false
    CALCULATE_BEADS_PER_DROP(
        FINAL_QC_SE.out.qc_stats,
        ch_images_pulled
    )

    ch_indexes = Channel.empty()
    ch_final_bam = Channel.empty()

    if (params.atac.barcodedTn5) {
        if (params.atac.barcodedTn5Config) {
            // Use the config file to reorganize the fastq-TI pairs into samples
            ch_ci_parse = VALIDATE_TI_CONFIG.out.config
                .splitCsv(header:true)
                .map { row -> tuple(row.fastq + '-' + row.sequence, row.sample) }

            ch_bam_parse = ch_ci_parse
                .combine(FINAL_BAM_MERGE.out.final_bam
                    .map { i -> [i[0], i[1]] }, by:0)
                .combine(BWA_ALIGNMENT.out.bam
                    .map { i -> [i[0], i[1]] }, by:0)
                .combine(MARK_DUPLICATES.out.bam
                    .map { i -> [i[0], i[1]] }, by:0)
                .map { i -> i.drop(1) }
                .groupTuple()

            ch_fragments_parse = ch_ci_parse
                .combine(FINAL_FRAG_MERGE.out.final_frag_merge, by:0)
                .map { i -> i.drop(1) }
                .groupTuple()

            ch_stats_parse = VALIDATE_TI_CONFIG.out.config
                .splitCsv(header:true)
                .map { row -> tuple(row.fastq + '-' + row.sequence, row.sample) }
                .combine(FINAL_QC_SE.out.qc_compile, by:0)
                .combine(MARK_DUPLICATES.out.stats, by:0)
                .combine(CALCULATE_BEADS_PER_DROP.out.cells, by:0)
                .map { tuple(it[1], it[2], it[3], it[4], it[5]) }
                .groupTuple()

            ch_bead_summary = VALIDATE_TI_CONFIG.out.config
                .splitCsv(header:true)
                .map { row -> tuple(row.fastq + '-' + row.sequence, row.sample) }
                .combine(FINAL_QC_SE.out.qc_stats, by:0)
                .combine(DETERMINE_BARCODE_ALLOWLIST.out.quant, by:0)
                .combine(ch_determine_barcode_merges_params, by:0)
                .map { tuple(it[1], it[2], it[3], it[4]) }
                .groupTuple()
        } else {
            // Gather all of the fastq-TI pairs into the SuperloadedSample
            ch_bam_parse = FINAL_BAM_MERGE.out.final_bam
                .map { i -> [i[0], i[1]] }
                .combine(BWA_ALIGNMENT.out.bam
                    .map { i -> [i[0], i[1]] }, by:0)
                .combine(MARK_DUPLICATES.out.bam
                    .map { i -> [i[0], i[1]] }, by:0)
                .combine(FINAL_FRAG_MERGE.out.final_frag_merge, by:0)
                .map { i -> ['SuperloadedSample', i.drop(1)].flatten() }
                .map { i -> [i[0], i[1], i[2], i[3]] }
                .groupTuple()

            ch_fragments_parse = FINAL_FRAG_MERGE.out.final_frag_merge
                .map { i -> ['SuperloadedSample', i.drop(1)].flatten() }
                .groupTuple()

            ch_stats_parse = FINAL_QC_SE.out.qc_compile
                .combine(MARK_DUPLICATES.out.stats, by:0)
                .combine(CALCULATE_BEADS_PER_DROP.out.cells, by:0)
                .map { i -> ['SuperloadedSample', i.drop(1)].flatten() }
                .groupTuple()

            ch_bead_summary = FINAL_QC_SE.out.qc_stats
                .combine(DETERMINE_BARCODE_ALLOWLIST.out.quant, by:0)
                .combine(ch_determine_barcode_merges_params, by:0)
                .map { i -> ['SuperloadedSample', i.drop(1)].flatten() }
                .groupTuple()
        }

        // Compile bam files from fastq-TI pairs
        COMPILE_ALIGNMENTS(
            ch_bam_parse,
            ch_images_pulled
        )

        // Compile fragments files from fastq-TI pairs
        COMPILE_FRAGMENTS(
            ch_fragments_parse,
            ch_images_pulled
        )

        // Compile stats files from fastq-TI pairs
        COMPILE_QC_STATS(
            ch_stats_parse,
            ch_images_pulled
        )

        // Create bead filtration summary table
        BEAD_FILT_SUMMARY(
            ch_bead_summary,
            ch_fastq_ti_pairs_read_counts.first(),
            ch_images_pulled
        )

        ch_marked_bam = COMPILE_ALIGNMENTS.out.marked_bam
        ch_final_bam = COMPILE_ALIGNMENTS.out.bam
        ch_final_fragments = COMPILE_FRAGMENTS.out.fragments
        ch_final_stats = COMPILE_QC_STATS.out.stats
        ch_report_stats = COMPILE_QC_STATS.out.stats.mix(FINAL_QC_SE.out.qc_stats)
        ch_final_dedup = COMPILE_QC_STATS.out.dedup
        ch_final_cells = COMPILE_QC_STATS.out.cells
        ch_report_bead_index_summary = BEAD_FILT_SUMMARY.out.index_summary
        ch_report_bead_sample_summary = BEAD_FILT_SUMMARY.out.sample_summary
    } else {
        ch_marked_bam = MARK_DUPLICATES.out.bam.map { tuple(it[0], it[1], it[2]) }
        ch_final_bam = FINAL_BAM_MERGE.out.final_bam
        ch_final_fragments = FINAL_FRAG_MERGE.out.final_frag_merge
        ch_final_stats = FINAL_QC_SE.out.qc_stats
        ch_report_stats = FINAL_QC_SE.out.qc_stats
        ch_final_dedup = MARK_DUPLICATES.out.stats
        ch_final_cells = CALCULATE_BEADS_PER_DROP.out.cells
        ch_report_bead_index_summary = Channel.empty()
        ch_report_bead_sample_summary = Channel.empty()
    }

    // Collect alignment stats
    SUMMARIZE_ALIGNMENTS(
        ch_marked_bam,
        ch_reference_fasta,
        ch_images_pulled
    )

    // Compute TSS matrix
    COMPUTE_TSS_MATRIX(
        ch_final_bam,
        ch_tss,
        ch_images_pulled
    )

    // Call peaks of fragment coverage
    CALL_PEAKS(
        ch_final_bam,
        ch_images_pulled
    )

    // Clean peaks to set width and remove low probability peaks and peaks in blocklist
    CLEAN_PEAKS(
        CALL_PEAKS.out.peaks,
        ch_blocklist,
        ch_sizes,
        ch_images_pulled
    )

    // Calculate fraction of reads in called peaks
    FRACTION_OF_READS_IN_PEAKS(
        ch_final_bam
            .join(CLEAN_PEAKS.out.peaks, by:0),
        ch_images_pulled
    )

    // Calculate fraction of reads in TSS
    FRACTION_OF_READS_IN_TSS(
        ch_final_bam,
        ch_tss,
        ch_images_pulled
    )

    // Calculate insert size metrics
    CALCULATE_INSERT_SIZE_METRICS(
        ch_final_bam,
        ch_images_pulled
    )

    ch_report_crosstalk = Channel.empty()
    ch_html_crosstalk = Channel.empty()
    ch_metrics_crosstalk = Channel.empty()

    if (params.atac.mixed) {
        ch_species = ch_gtf
            .splitText()
            .filter { !(it =~ /#/) }
            .map { it.split('\\.')[0] }
            .unique()
            .collect()

        // Calculate mixed species stats
        SUMMARIZE_MIXED_SPECIES(
            ch_final_stats,
            ch_species,
            ch_images_pulled
        )

        ch_report_crosstalk = SUMMARIZE_MIXED_SPECIES.out.stats
        ch_html_crosstalk = SUMMARIZE_MIXED_SPECIES.out.stats
        ch_metrics_crosstalk = SUMMARIZE_MIXED_SPECIES.out.stats
    }

    // Generate count matrix for called peaks
    MAKE_COUNT_MATRIX(
        ch_final_stats
            .join(ch_final_bam, by:0)
            .combine(CLEAN_PEAKS.out.peaks, by:0),
        ch_images_pulled
    )

    // Calculate fragment saturation curve
    SEQUENCE_SATURATION(
        ch_final_fragments
            .map { i -> [i[0], i[1]] }
            .join(SUMMARIZE_ALIGNMENTS.out.stats
                .map { i -> [i[0], i[3]] }, by: 0),
        assay,
        ch_images_pulled
    )

    // Use ArchR
    ARCHR(
        ch_gtf,
        ch_archr_ref,
        ch_blocklist,
        ch_final_bam.join(CLEAN_PEAKS.out.peaks, by:0),
        ch_images_pulled
    )

    // Calculate TSS enrichment score for use in aggregate metrics
    TSS_ENRICHMENT(
        COMPUTE_TSS_MATRIX.out.tss_matrix,
        ch_images_pulled
    )

    // Combine useful metrics into a single file
    AGGREGATE_METRICS(
        Channel.empty()
            .mix(FASTQC.out.zip.map { i  -> i[1] } .collect() .ifEmpty { [] },
                CUTADAPT_HEADCROP.out.readlength.map { i  -> i[1] } .collect() .ifEmpty { [] },
                SUMMARIZE_ALIGNMENTS.out.stats.map { i  -> i[1, 2, 3] } .collect() .ifEmpty { [] },
                ch_final_dedup.map { i  -> i[1] } .collect() .ifEmpty { [] },
                DETERMINE_BARCODE_ALLOWLIST.out.allowlist.map { i  -> i[1] } .collect() .ifEmpty { [] },
                DETERMINE_BARCODE_ALLOWLIST.out.quant.map { i  -> i[1] } .collect() .ifEmpty { [] },
                ch_determine_barcode_merges_params.map { i  -> i[1] } .collect() .ifEmpty { [] },
                ch_determine_barcode_merges_implicated_barcodes.map { i  -> i[1] } .collect() .ifEmpty { [] },
                TSS_ENRICHMENT.out.tss_enrichment.map { i  -> i[1] } .collect() .ifEmpty { [] },
                FRACTION_OF_READS_IN_PEAKS.out.frip.map { i  -> i[1] } .collect() .ifEmpty { [] },
                FRACTION_OF_READS_IN_TSS.out.frit.map { i  -> i[1] } .collect() .ifEmpty { [] },
                FINAL_QC_SE.out.qc_stats.map { i  -> i[1] } .collect() .ifEmpty { [] },
                ch_metrics_crosstalk.map { i  -> i[1, 2] } .collect() .ifEmpty { [] },
                CLEAN_PEAKS.out.peaks.map { i  -> i[1] } .collect() .ifEmpty { [] },
                FINAL_BAM_MERGE.out.count.map { i  -> i[1] } .collect() .ifEmpty { [] },
                DEAD.out.count.map { i  -> i[1] } .collect() .ifEmpty { [] },
                BWA_ALIGNMENT.out.count.map { i  -> i[1] } .collect() .ifEmpty { [] },
                REANNOTATE_BAM.out.count.map { i  -> i[1] } .collect() .ifEmpty { [] },
                ANNOTATE_FRAGMENTS.out.stats.map { i  -> i[2] } .collect() .ifEmpty { [] },
                ASSEMBLE_FRAGMENTS.out.stats.map { i  -> i[2] } .collect() .ifEmpty { [] },
                REANNOTATE_FRAGMENTS.out.stats.map { i  -> i[2] } .collect() .ifEmpty { [] },
                MERGE_REANN_READ_COUNTS.out.count.map { i  -> i[1] } .collect() .ifEmpty { [] },
                ch_fastq_ti_pairs_read_counts .ifEmpty { [] },
                ch_report_bead_sample_summary.map { i  -> i[1] } .collect() .ifEmpty { [] })
            .collect(),
        assay,
        ch_images_pulled
    )

    // Build contents for export to reporting step
    BUILD_REPORT_CONTENTS(
        Channel.empty()
            .mix(CALCULATE_INSERT_SIZE_METRICS.out.metrics.map { i -> i[1] }.collect() .ifEmpty { [] },
                AGGREGATE_METRICS.out.summary .collect() .ifEmpty { [] },
                ch_validated_ti_config .ifEmpty { [] },
                ch_fastq_ti_pairs_read_counts .ifEmpty { [] }
                )
            .collect(),
        ch_images_pulled
    )

    if (params.atac.barcodedTn5) {
        ch_index_counts = BUILD_REPORT_CONTENTS.out.index_counts
        ch_messages_file = TI_WARNING_MESSAGES.out.messages
    } else {
        ch_index_counts = Channel.empty()
        ch_messages_file = ch_messages
    }

    // Generate the report
    GENERATE_REPORT(
        BUILD_REPORT_CONTENTS.out.metric_summary .collect() .ifEmpty { [] },
        BUILD_REPORT_CONTENTS.out.pipeline_summary .collect() .ifEmpty { [] },
        FASTQC.out.zip.map { i  -> i[1] } .collect() .ifEmpty { [] },
        DEAD.out.count.map { i  -> i[1] } .collect() .ifEmpty { [] },
        SUMMARIZE_ALIGNMENTS.out.stats.map { i  -> i[1, 2, 3] } .collect() .ifEmpty { [] },
        ch_report_stats.map { i  -> i[1] } .collect() .ifEmpty { [] },
        ch_determine_barcode_merges_params.map { i  -> i[1] } .collect() .ifEmpty { [] },
        ch_determine_barcode_merges_implicated_barcodes.map { i  -> i[1] } .collect() .ifEmpty { [] },
        DETERMINE_BARCODE_ALLOWLIST.out.quant.map { i  -> i[1] } .collect() .ifEmpty { [] },
        DETERMINE_BARCODE_ALLOWLIST.out.allowlist.map { i  -> i[1] } .collect() .ifEmpty { [] },
        ch_final_cells.map { i  -> i[1] } .collect() .ifEmpty { [] },
        TSS_ENRICHMENT.out.tss_enrichment.map { i  -> i[1] } .collect() .ifEmpty { [] },
        CALCULATE_INSERT_SIZE_METRICS.out.metrics.map { i  -> i[1] }.collect() .ifEmpty { [] },
        SEQUENCE_SATURATION.out.results.map { i  -> i[1] } .collect() .ifEmpty { [] },
        ch_report_crosstalk.map { i  -> i[1, 2] } .collect() .ifEmpty { [] },
        ch_report_bead_index_summary.map { i  -> i[1] } .collect() .ifEmpty { [] },
        ch_report_bead_sample_summary.map { i  -> i[1] } .collect() .ifEmpty { [] },
        ARCHR.out.clusterinfo .collect() .ifEmpty { [] },
        AGGREGATE_METRICS.out.summary .collect() .ifEmpty { [] },
        ch_fastq_ti_pairs_read_counts .collect() .ifEmpty { [] },
        ch_validated_ti_config .ifEmpty { [] },
        ch_messages_file .ifEmpty { [] },
        ch_sample_map .ifEmpty { [] },
        ch_images_pulled
    )

    emit:
    MAKE_COUNT_MATRIX.out.matrix // tuple: [ sample_id, column_names.txt, row_names.txt, read_counts.txt]
}
