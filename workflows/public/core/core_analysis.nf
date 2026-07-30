/*
Analyze single-cell data with core modules
*/

include { MERGE_LANES                         } from "${params.omnition.modulesDir}/public/atac/merge_lanes.nf"
include { VALIDATE_FASTQS                     } from "${params.omnition.modulesDir}/public/core/validate_fastqs.nf"
include { COMPILE_FASTQ_VALIDATIONS           } from "${params.omnition.modulesDir}/public/core/compile_fastq_validations.nf"
include { PUBLISH_PARAMETERS                  } from "${params.omnition.modulesDir}/public/core/publish_parameters.nf"
include { FASTQC                              } from "${params.omnition.modulesDir}/public/core/fastqc.nf"
include { CUTADAPT_TRIM                       } from "${params.omnition.modulesDir}/public/core/cutadapt_trim.nf"
include { TI_DEBARCODER_CONFIG                } from "${params.omnition.modulesDir}/public/core/ti_debarcoder_config.nf"
include { PREPARE_ANTIBODY_FILE               } from "${params.omnition.modulesDir}/public/cytometry/prepare_antibody_file.nf"
include { DEBARCODER                          } from "${params.omnition.modulesDir}/public/core/debarcoder.nf"
include { DEBARCODER as DEBARCODER_TI         } from "${params.omnition.modulesDir}/public/atac/debarcoder.nf"
include { DEBARCODER_CITESEQ                  } from "${params.omnition.modulesDir}/public/citeseq/debarcoder_citeseq.nf"
include { GENERATE_BEAD_LIST                  } from "${params.omnition.modulesDir}/public/core/generate_bead_list.nf"

workflow CORE_ANALYSIS {
    take:
    ch_fastq // tuple: [ sampleId, raw R1 FASTQ files, raw R2 FASTQ files ]
    val_assay // which additional assay is being ran (string)
    ch_params_file // path: parameters file
    ch_images_pulled // boolean: true if Singularity pull completed

    main:
    ch_merge_lanes_input = ch_fastq
    ch_merge_lanes_complete = Channel.empty()
    ch_debarcoded_fastq = Channel.empty()
    ch_debarcoded_adt_fastq = Channel.empty()
    ch_debarcoded_edges = Channel.empty()
    ch_debarcoded_count = Channel.empty()

    // Merge lanesplit fastqs
    MERGE_LANES(
        ch_merge_lanes_input,
        ch_images_pulled
    )
    ch_merge_lanes_complete = MERGE_LANES.out.complete

    VALIDATE_FASTQS(
        ch_merge_lanes_complete
            .map { out -> out.flatten() }
            .map { out -> tuple(out[0], [out[1]], [out[2]]) },
        ch_images_pulled
    )

    compile_fastq_validations = VALIDATE_FASTQS.out.error_files
            .flatMap { it.getClass() == ArrayList ? it[1] : it }
            .collect()

    // Sort error files into passed and failed and halt if any failed
    COMPILE_FASTQ_VALIDATIONS(
        compile_fastq_validations,
        ch_images_pulled
    )

    // Pass the empty files if validations passed
    ch_val_pass = COMPILE_FASTQ_VALIDATIONS.out.val_pass.flatten()
        .map { item -> tuple(item.name - "_errors.txt", item) }

    // Publish parameters
    PUBLISH_PARAMETERS(
        ch_params_file.ifEmpty { [] },
        ch_images_pulled
    )

    ch_complete_fastq = ch_merge_lanes_complete
    ch_final_val_pass = ch_val_pass

    FASTQC(
        ch_complete_fastq
            .join(ch_val_pass),
        ch_images_pulled
    )

    ch_debarcoder_input_fastq = ch_complete_fastq

    ch_cutadapt_metrics = Channel.empty()
    ch_cutadapt_report = Channel.empty()

    if (!params.catac && !params.atac) {
        CUTADAPT_TRIM(
            ch_complete_fastq
                .join(ch_val_pass),
            ch_images_pulled
        )

        ch_cutadapt_metrics = CUTADAPT_TRIM.out.readlength
            .mix(CUTADAPT_TRIM.out.count)
        ch_cutadapt_report = CUTADAPT_TRIM.out.log
        ch_debarcoder_input_fastq = CUTADAPT_TRIM.out.fastq
    }

    if (params.cytometry) {
        ch_antibody_file   = Channel.fromPath(params.cytometry.adtFile)
        // Clean the antibody tag file by converting to unix and removing special characters
        PREPARE_ANTIBODY_FILE(
            ch_antibody_file,
            ch_images_pulled
        )
        ch_prepared_antibody_file = PREPARE_ANTIBODY_FILE.out.cleaned
    } else {
        ch_prepared_antibody_file = Channel.empty() // No ADT config for non-cytometry
    }

    ti_csv_ch = Channel.empty()

    if (params.catac) {
        // Create CSV file of TIs from Nextflow parameter
        assay_ti = params.catac.ti
        ti_csv_ch = Channel
        .from("${ assay_ti }")
        .map { it.replaceAll('\\s', '') } // groovylint-disable-line
        .map { it.replaceAll(',', '\n') } // groovylint-disable-line
        .map { it.replaceAll(':', ',') } // groovylint-disable-line
        .map { it.replaceAll('sequence,', '') } // groovylint-disable-line
        .map { it.replaceAll('\\[', '') } // groovylint-disable-line
        .map { it.replaceAll('\\]', '') } // groovylint-disable-line
        .collectFile(name:'TI.csv', seed: 'name,sequence\n', cache: 'lenient')
        .first()

        // Create a TI config for DEBARCODER
        TI_DEBARCODER_CONFIG(
            ch_debarcoder_input_fastq
                .combine(ti_csv_ch)
                .map { out -> tuple(out[0], out[2]) }
                .join(ch_val_pass),
            ch_images_pulled
        )
        if (params.catac) {
            ch_debarcoder_config = TI_DEBARCODER_CONFIG.out.config
        } else {
            ch_debarcoder_config = TI_DEBARCODER_CONFIG.out.config_r1
        }

        // Parse and correct capture oligos
        DEBARCODER_TI(
            ch_debarcoder_input_fastq
                .join(ch_debarcoder_config),
            ch_images_pulled
        )

        ch_debarcoded_fastq = DEBARCODER_TI.out.fastq
        ch_debarcoded_count = DEBARCODER_TI.out.count
        ch_debarcoded_edges = Channel.empty() // No edges output for TI debarcoder
        ch_debarcoded_adt_fastq = Channel.empty() // No ADT fastq
    } else {
        if (params.rna && params.cytometry) {
            // Call CITE-seq version
            DEBARCODER_CITESEQ(
                ch_debarcoder_input_fastq
                    .combine(PREPARE_ANTIBODY_FILE.out.adt_debarcoder_ref),
                ch_images_pulled
            )
            ch_debarcoded_fastq = DEBARCODER_CITESEQ.out.fastq
            ch_debarcoded_adt_fastq = DEBARCODER_CITESEQ.out.adt_fastq
            ch_debarcoded_edges = DEBARCODER_CITESEQ.out.edges
            ch_debarcoded_count = DEBARCODER_CITESEQ.out.count
        } else {
            // Parse and correct capture oligos (RNA-only)
            DEBARCODER(
                ch_debarcoder_input_fastq,
                ch_images_pulled
            )
            ch_debarcoded_fastq = DEBARCODER.out.fastq
            ch_debarcoded_adt_fastq = Channel.empty() // No ADT fastq for RNA-only
            ch_debarcoded_edges = DEBARCODER.out.edges
            ch_debarcoded_count = DEBARCODER.out.count
        }
    }

    ch_bead_list = Channel.empty()
    if (params.rna || params.cytometry) {
        // Generate a list of beads to be merged
        GENERATE_BEAD_LIST(
            ch_debarcoded_fastq,
            ch_images_pulled
        )
        ch_bead_list = GENERATE_BEAD_LIST.out.beads
    }

    if (params.atac) {
        ch_sample_map = ch_debarcoded_fastq
        .map { item  -> "${item[0]},${item[0]}" }
        .collectFile(name:'sample_map.csv', seed: 'sample,fastq', newLine: true, // groovylint-disable-line
            sort: true, cache: 'lenient', storeDir: "${params.atac.reportsDir}")
        .first()
    }

    emit:
    input_params_yaml       = PUBLISH_PARAMETERS.out.input_params_yaml
    params_yaml             = PUBLISH_PARAMETERS.out.params_yaml
    command_txt             = PUBLISH_PARAMETERS.out.command_txt
    fastqc_zip              = FASTQC.out.zip
    fastqc_count            = FASTQC.out.counts
    fastqc_sequence_traces  = FASTQC.out.sequence_traces
    fastqc_quality_scores   = FASTQC.out.quality_scores
    val_pass                = ch_val_pass
    complete_fastq          = ch_merge_lanes_complete
    cutadapt_metrics        = ch_cutadapt_metrics
    cutadapt_report         = ch_cutadapt_report
    prepared_antibody_file  = ch_prepared_antibody_file // path: antibodies.csv or empty
    debarcoder_fastq        = ch_debarcoded_fastq
    debarcoder_edges        = ch_debarcoded_edges
    debarcoder_count        = ch_debarcoded_count
    debarcoder_adt_fastq    = ch_debarcoded_adt_fastq
    ti_csv                  = ti_csv_ch
    bead_list               = ch_bead_list
}
