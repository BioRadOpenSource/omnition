/*
Analyze single-cell ATAC data
*/

include { CHECK_MITO_CONTIG                  } from "${params.omnition.modulesDir}/public/atac/check_mito_contig.nf"
include { CUTADAPT_HEADCROP                  } from "${params.omnition.modulesDir}/public/atac/cutadapt_headcrop.nf"
include { VALIDATE_TI_CONFIG                 } from "${params.omnition.modulesDir}/public/atac/validate_ti_config.nf"
include { CHECK_TI_COUNTS                    } from "${params.omnition.modulesDir}/public/atac/check_ti_counts.nf"
include { CHECK_TI_COUNTS_SUPERLOADED        } from "${params.omnition.modulesDir}/public/atac/check_ti_counts_superloaded.nf"
include { TI_WARNING_MESSAGES                } from "${params.omnition.modulesDir}/public/atac/ti_warning_messages.nf"
include { SPLIT_FASTQ                        } from "${params.omnition.modulesDir}/public/atac/split_fastq.nf"
include { COMPRESS_SPLIT_FASTQS              } from "${params.omnition.modulesDir}/public/atac/compress_split_fastqs.nf"
include { BWA_ALIGNMENT                      } from "${params.omnition.modulesDir}/public/atac/bwa_alignment.nf"
include { MARK_DUPLICATES                    } from "${params.omnition.modulesDir}/public/atac/mark_duplicates.nf"
include { SUMMARIZE_ALIGNMENTS               } from "${params.omnition.modulesDir}/public/atac/summarize_alignments.nf"
include { SPLIT_BAM                          } from "${params.omnition.modulesDir}/public/atac/split_bam.nf"
include { ASSEMBLE_FRAGMENTS                 } from "${params.omnition.modulesDir}/public/atac/assemble_fragments.nf"
include { DETERMINE_BARCODE_ALLOWLIST        } from "${params.omnition.modulesDir}/public/atac/determine_barcode_allowlist.nf"
include { COMPUTE_DECONVOLUTION_STAT_CHR     } from "${params.omnition.modulesDir}/public/atac/compute_deconvolution_stat_chr.nf"
include { DETERMINE_BARCODE_MERGES           } from "${params.omnition.modulesDir}/public/atac/determine_barcode_merges.nf"
include { REANNOTATE_FRAGMENTS               } from "${params.omnition.modulesDir}/public/atac/reannotate_fragments.nf"
include { REANNOTATE_BAM                     } from "${params.omnition.modulesDir}/public/atac/reannotate_bam.nf"
include { FINAL_BAM_MERGE                    } from "${params.omnition.modulesDir}/public/atac/final_bam_merge.nf"
include { MERGE_REANN_READ_COUNTS            } from "${params.omnition.modulesDir}/public/atac/merge_reann_read_counts.nf"
include { FINAL_FRAG_MERGE                   } from "${params.omnition.modulesDir}/public/atac/final_frag_merge.nf"
include { ASSEMBLE_BASIC_QC                  } from "${params.omnition.modulesDir}/public/atac/assemble_basic_qc.nf"
include { FINAL_QC_SE                        } from "${params.omnition.modulesDir}/public/atac/final_qc_se.nf"
include { SEQUENCE_SATURATION                } from "${params.omnition.modulesDir}/public/core/sequence_saturation.nf"
include { COMPILE_ALIGNMENTS                 } from "${params.omnition.modulesDir}/public/atac/compile_alignments.nf"
include { COMPILE_FRAGMENTS                  } from "${params.omnition.modulesDir}/public/atac/compile_fragments.nf"
include { COMPILE_QC_STATS                   } from "${params.omnition.modulesDir}/public/atac/compile_qc_stats.nf"
include { BEAD_FILT_SUMMARY                  } from "${params.omnition.modulesDir}/public/atac/bead_filt_summary.nf"
include { COMPUTE_TSS_MATRIX                 } from "${params.omnition.modulesDir}/public/atac/compute_tss_matrix.nf"
include { CALL_PEAKS                         } from "${params.omnition.modulesDir}/public/atac/call_peaks.nf"
include { CLEAN_PEAKS                        } from "${params.omnition.modulesDir}/public/atac/clean_peaks.nf"
include { FRACTION_OF_READS_IN_PEAKS         } from "${params.omnition.modulesDir}/public/atac/fraction_of_reads_in_peaks.nf"
include { FRACTION_OF_READS_IN_TSS           } from "${params.omnition.modulesDir}/public/atac/fraction_of_reads_in_tss.nf"
include { CALCULATE_INSERT_SIZE_METRICS      } from "${params.omnition.modulesDir}/public/atac/calculate_insert_size_metrics.nf"
include { CALCULATE_BEADS_PER_DROP           } from "${params.omnition.modulesDir}/public/atac/calc_beads_per_drop.nf"
include { SUMMARIZE_MIXED_SPECIES            } from "${params.omnition.modulesDir}/public/atac/summarize_mixed_species.nf"
include { MAKE_COUNT_MATRIX                  } from "${params.omnition.modulesDir}/public/atac/make_count_matrix.nf"
include { ARCHR                              } from "${params.omnition.modulesDir}/public/atac/archr.nf"
include { TSS_ENRICHMENT                     } from "${params.omnition.modulesDir}/public/atac/tss_enrichment.nf"
include { AGGREGATE_METRICS                  } from "${params.omnition.modulesDir}/public/core/aggregate_metrics.nf"
include { BUILD_REPORT_CONTENTS              } from "${params.omnition.modulesDir}/public/atac/build_report_contents.nf"
include { RENDER_REPORT                      } from "${params.omnition.modulesDir}/public/atac/render_report.nf"

workflow ATAC_COMBINATORIAL_ANALYSIS {
    take:
    ch_input_params_yaml // path: input params yaml file
    ch_params_yaml // path: params yaml file
    ch_command_txt // path: command txt file
    ch_fastqc_zip // tuple: [ sampleId, fastqc zip files ]
    ch_fastqc_count // tuple: [ sampleId, fastqc per base sequence content ]
    ch_fastqc_sequence_traces // tuple: [ sampleId, fastqc sequence traces ]
    ch_fastqc_quality_scores // tuple: [ sampleId, fastqc quality scores ]
    ch_debarcoded_fastq // tuple: [ sampleId, debarcoded R1
    ch_debarcoder_count // tuple: [ sampleId, debarcoder counts ]
    ti_csv_ch // CSV config file of fastq and TIs
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
    ch_sample_sheet_map // map: sample sheet metadata
    ch_ti_map // map: TI metadata
    ch_images_pulled // boolean: true if Singularity pull completed

    main:

    // Check mitochondrial contig
    CHECK_MITO_CONTIG(
        ch_sizes,
        ch_images_pulled
    )

    ch_validated_ti_config = Channel.empty()

    ch_ti_len = 0
    if (params.catac.barcodedTn5Config) {
        // Validate TI config
        VALIDATE_TI_CONFIG(
            ch_TI_config,
            ch_images_pulled,
            ti_csv_ch
        )

        ch_ti_len = VALIDATE_TI_CONFIG.out.tilen.first().text
        ch_validated_ti_config = VALIDATE_TI_CONFIG.out.config

        // Check counts for each TI
        CHECK_TI_COUNTS(
            ch_debarcoded_fastq,
            ti_csv_ch,
            VALIDATE_TI_CONFIG.out.config.first(),
            ch_images_pulled
        )

        // Combine Warning messages into one output file
        TI_WARNING_MESSAGES(
            CHECK_TI_COUNTS.out.messages.collect(),
            ch_messages
        )

        // Sort error files into passed and failed
        compile_ti_errors = CHECK_TI_COUNTS.out.error_files
            .collect()
            .flatten()
            .branch {
                passed: it.size() == 0
                failed: it.size() != 0
            }

        // Halt the pipeline if any tis failed
        compile_ti_errors.failed
            .flatMap { it }
            .collectFile(name:'TI_run_errors.txt', sort: { it.name },
                newLine: false, storeDir: "${params.catac.reportsDir}")
            .subscribe {
                if (params.catac.tiErrorOverride) {
                    log.info("INFO: [$params.catac.assay] TI errors detected, but overridden. " \
                            + "Please see ${params.catac.reportsDir}TI_run_errors.txt for more details.")
                } else {
                    log.error("ERROR: [$params.catac.assay] TI error(s) detected. " \
                            + "Please see ${params.catac.reportsDir}TI_run_errors.txt for more details.") // groovylint-disable-line
                    exit(1)
                }
            }

        ch_fastq_ti_pairs_read_counts = CHECK_TI_COUNTS.out.split
            .flatMap { it }
            .collectFile(name:'fastqTIreadcounts.csv', seed: 'sample,fastq,ti,sequence,count',
                newLine: true, cache: 'lenient', storeDir: "${params.catac.reportsDir}")

        if (params.catac.tiErrorOverride) {
            ch_error_pass = Channel.value(true)
        } else {
            ch_error_pass = compile_ti_errors.passed.first()
        }
    } else {
        // Check counts for each TI
        CHECK_TI_COUNTS_SUPERLOADED(
            ch_debarcoded_fastq,
            ti_csv_ch,
            ch_images_pulled
        )

        ch_ti_len = CHECK_TI_COUNTS_SUPERLOADED.out.tilen.first().text

        // Combine Warning messages into one output file
        TI_WARNING_MESSAGES(
            CHECK_TI_COUNTS_SUPERLOADED.out.messages.collect(),
            ch_messages
        )

        ch_fastq_ti_pairs_read_counts = CHECK_TI_COUNTS_SUPERLOADED.out.split
            .flatMap { it }
            .collectFile(name:'fastqTIreadcounts_superloaded.csv', seed: 'sample,fastq,ti,sequence,count', // groovylint-disable-line
                newLine: true, cache: 'lenient', storeDir: "${params.catac.reportsDir}")

        ch_error_pass = CHECK_TI_COUNTS_SUPERLOADED.out.passed.first()
    }

    // Split each fastq by TI
    SPLIT_FASTQ(
        ch_debarcoded_fastq
            .combine(ch_fastq_ti_pairs_read_counts),
        ch_error_pass,
        ch_images_pulled
    )

    // Split each fastq by TI
    COMPRESS_SPLIT_FASTQS(
        SPLIT_FASTQ.out.split
            .transpose()
            .map { out -> tuple(out[0], out[1].name.split('\\.')[0] - "${out[0]}-" - '_R1' - '_R2', out[1]) }
            .groupTuple(by: [0, 1]),
        ch_images_pulled
    )

    ch_sample_map = ch_fastq_ti_pairs_read_counts
        .splitCsv(header:true)
        .map { row -> "${row.sample},${row.fastq},${row.ti}" }
        .collectFile(name:'sample_map.csv', seed: 'sample,fastq,ti', newLine: true,
            sort: true, cache: 'lenient', storeDir: "${params.catac.reportsDir}")
        .first()

    ch_complete_fastq = COMPRESS_SPLIT_FASTQS.out.split

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
        ch_split_bam_formatted
            .combine(ch_blocklist),
        ch_images_pulled
    )

    // Bead filtration with knee-calling
    DETERMINE_BARCODE_ALLOWLIST(
        ASSEMBLE_FRAGMENTS.out.bead_counts
            .groupTuple(),
        ch_images_pulled
    )

    // Calculate chromosome-specific stats
    COMPUTE_DECONVOLUTION_STAT_CHR(
        ASSEMBLE_FRAGMENTS.out.assemble_fragments
            .combine(DETERMINE_BARCODE_ALLOWLIST.out.allowlist, by: 0),
        ch_images_pulled
    )

    ch_in_determine_merges = DETERMINE_BARCODE_ALLOWLIST.out.quant
        .join(DETERMINE_BARCODE_ALLOWLIST.out.allowlist, by: 0)
        .join(DETERMINE_BARCODE_ALLOWLIST.out.params, by: 0)
        .map { i  -> [ i[0].substring(0, i[0].lastIndexOf('-')), i[1], i[2], i[3] ] }
        .groupTuple()
        .join(COMPUTE_DECONVOLUTION_STAT_CHR.out.overlap_count
            .filter { i -> !(i[1] =~ /${params.catac.mitoContig}/) }
            .map { i  -> [ i[0].substring(0, i[0].lastIndexOf('-')), i[2] ] }
            .groupTuple())

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
        ASSEMBLE_FRAGMENTS.out.assemble_fragments
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
            .filter { i -> !(i[1] =~ /${params.catac.mitoContig}/) }
            .map { tuple(it[0], it[2], it[3]) }
            .groupTuple(),
        ch_images_pulled
    )

    MERGE_REANN_READ_COUNTS(
        ASSEMBLE_FRAGMENTS.out.count
            .groupTuple()
            .join(
                REANNOTATE_BAM.out.count
                    .groupTuple(),
                by:0
            )
            .join(
                REANNOTATE_FRAGMENTS.out.count
                    .groupTuple(),
                by:0
            ),
        ch_images_pulled
    )

    // Merge chromosome fragment files back together
    FINAL_FRAG_MERGE(
        REANNOTATE_FRAGMENTS.out.reannotate_fragments
            .filter { i -> !(i[1].name =~ /${params.catac.mitoContig}/) }
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

    if (params.catac.barcodedTn5Config) {
        // Use the config file to reorganize the fastq-TI pairs into samples
        ch_ci_parse = VALIDATE_TI_CONFIG.out.config
            .splitCsv(header:true)
            .map { row -> tuple(row.fastq + '-' + row.sequence, row.sample) }

        ch_bead_fragments_parse = ch_ci_parse
            .combine(ASSEMBLE_FRAGMENTS.out.assemble_fragments, by:0)
            .map { i -> i.drop(1) }
            .groupTuple()
            .map { i -> [i[0], i[2]] }

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
        ch_bead_fragments_parse = ASSEMBLE_FRAGMENTS.out.assemble_fragments
            .map { i -> ['SuperloadedSample', i[2]] }
            .groupTuple()

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

    // Collect alignment stats
    SUMMARIZE_ALIGNMENTS(
        COMPILE_ALIGNMENTS.out.marked_bam,
        ch_images_pulled
    )

    // Compute TSS matrix
    COMPUTE_TSS_MATRIX(
        COMPILE_ALIGNMENTS.out.bam
            .combine(ch_tss),
        ch_images_pulled
    )

    // Call peaks of fragment coverage
    CALL_PEAKS(
        COMPILE_ALIGNMENTS.out.bam,
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
        COMPILE_QC_STATS.out.stats
            .join(COMPILE_ALIGNMENTS.out.bam, by:0)
            .combine(CLEAN_PEAKS.out.peaks, by:0),
        ch_images_pulled
    )

    // Calculate fraction of reads in called peaks
    FRACTION_OF_READS_IN_PEAKS(
        COMPILE_ALIGNMENTS.out.bam
            .join(CLEAN_PEAKS.out.peaks, by:0),
        ch_images_pulled
    )

    // Calculate fraction of reads in TSS
    FRACTION_OF_READS_IN_TSS(
        COMPILE_ALIGNMENTS.out.bam
            .combine(ch_tss),
        ch_images_pulled
    )

    // Calculate insert size metrics
    CALCULATE_INSERT_SIZE_METRICS(
        COMPILE_ALIGNMENTS.out.bam,
        ch_images_pulled
    )

    ch_report_crosstalk = Channel.empty()

    if (params.catac.mixed) {
        ch_species = ch_gtf
            .splitText()
            .filter { !(it =~ /#/) }
            .map { it.split('\\.')[0] } // groovylint-disable-line
            .unique()
            .collect()

        // Calculate mixed species stats
        SUMMARIZE_MIXED_SPECIES(
            COMPILE_QC_STATS.out.stats,
            ch_species,
            ch_images_pulled
        )

        ch_report_crosstalk = SUMMARIZE_MIXED_SPECIES.out.stats
    }

    // Calculate fragment saturation curve
    SEQUENCE_SATURATION(
        COMPILE_FRAGMENTS.out.fragments
            .map { i -> [i[0], i[1]] }
            .join(SUMMARIZE_ALIGNMENTS.out.stats
                .map { i -> [i[0], i[3]] }, by: 0),
        assay,
        ch_images_pulled
    )
    ch_seq_sat_out = SEQUENCE_SATURATION.out.results

    // Use ArchR
    ch_archr_output = Channel.empty()
    if (params.catac.archr) {
        ARCHR(
            COMPILE_ALIGNMENTS.out.bam.join(CLEAN_PEAKS.out.peaks, by:0)
                .combine(ch_gtf)
                .combine(ch_archr_ref)
                .combine(ch_blocklist),
            ch_images_pulled
        )
        ch_archr_output = ARCHR.out.clusterinfo
    }

    // Calculate TSS enrichment score for use in aggregate metrics
    TSS_ENRICHMENT(
        COMPUTE_TSS_MATRIX.out.tss_matrix,
        ch_images_pulled
    )

    ch_agg_metrics = Channel.empty()
            .mix(ch_fastqc_zip.map { i  -> i[1] },
                ch_fastqc_count.map { i -> i[1] },
                CUTADAPT_HEADCROP.out.readlength.map { i  -> i[1] },
                SUMMARIZE_ALIGNMENTS.out.stats.map { i  -> i[1, 2, 3] },
                COMPILE_QC_STATS.out.dedup.map { i  -> i[1] },
                DETERMINE_BARCODE_ALLOWLIST.out.allowlist.map { i  -> i[1] },
                DETERMINE_BARCODE_ALLOWLIST.out.quant.map { i  -> i[1] },
                DETERMINE_BARCODE_MERGES.out.metrics,
                TSS_ENRICHMENT.out.tss_metric,
                FRACTION_OF_READS_IN_PEAKS.out.frip.map { i  -> i[1] },
                FRACTION_OF_READS_IN_TSS.out.frit.map { i  -> i[1] },
                FINAL_QC_SE.out.qc_metrics,
                ch_report_crosstalk.map { i  -> i[1, 2] },
                CLEAN_PEAKS.out.summary.map { i  -> i[1] },
                FINAL_BAM_MERGE.out.count.map { i  -> i[1] },
                ch_debarcoder_count.map { i  -> i[1] },
                BWA_ALIGNMENT.out.count.map { i  -> i[1] },
                MERGE_REANN_READ_COUNTS.out.count.map { i  -> i[1] },
                ch_fastq_ti_pairs_read_counts .ifEmpty { [] },
                BEAD_FILT_SUMMARY.out.summary)
            .collect()

    ch_agg_metrics_input = ch_agg_metrics

    // Combine useful metrics into a single file
    AGGREGATE_METRICS(
        ch_agg_metrics_input,
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

    // Generate the report
    RENDER_REPORT(
        Channel.empty()
            .mix(ch_input_params_yaml,
                ch_params_yaml,
                ch_command_txt,
                BUILD_REPORT_CONTENTS.out.pipeline_summary,
                ch_fastqc_zip.map { i  -> i[1] },
                ch_debarcoder_count.map { i  -> i[1] },
                SUMMARIZE_ALIGNMENTS.out.stats.map { i  -> i[1, 2, 3] },
                COMPILE_QC_STATS.out.stats.map { i  -> i[1] },
                FINAL_QC_SE.out.qc_stats.map { i  -> i[1] },
                ch_determine_barcode_merges_params.map { i  -> i[1] },
                ch_determine_barcode_merges_implicated_barcodes.map { i  -> i[1] },
                DETERMINE_BARCODE_ALLOWLIST.out.quant.map { i  -> i[1] },
                DETERMINE_BARCODE_ALLOWLIST.out.allowlist.map { i  -> i[1] },
                COMPILE_QC_STATS.out.cells.map { i  -> i[1] },
                TSS_ENRICHMENT.out.tss_enrichment.map { i  -> i[1] },
                CALCULATE_INSERT_SIZE_METRICS.out.metrics.map { i  -> i[1] },
                ch_seq_sat_out.map { i  -> i[1] },
                ch_report_crosstalk.map { i  -> i[1, 2] },
                ch_archr_output.collect(),
                AGGREGATE_METRICS.out.summary,
                ch_fastq_ti_pairs_read_counts,
                ch_validated_ti_config ,
                TI_WARNING_MESSAGES.out.messages ,
                ch_sample_map)
            .collect(),
        ch_images_pulled
    )

    emit:
    MAKE_COUNT_MATRIX.out.matrix // tuple: [ sampleId, column_names.txt, row_names.txt, read_counts.txt]
}
