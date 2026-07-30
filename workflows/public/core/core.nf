/*
Validate parameters and execute core workflows
*/

include { samplesheetToList } from 'plugin/nf-schema'

// Import workflows
include { PULL_SINGULARITY_IMAGES } from "${params.omnition.workflowsDir}/public/core/pull_singularity_images.nf"
include { GET_FASTQ_FILES         } from "${params.omnition.workflowsDir}/public/core/get_fastq_files.nf"
include { CORE_ANALYSIS           } from "${params.omnition.workflowsDir}/public/core/core_analysis.nf"

workflow CORE {
    take:
    all_modules
    val_assay
    fastqLaneRegEx
    messages

    main:
    // Set global assay params
    params.core.assay                     = "Core"

    // Validate analysis workflow
    GET_FASTQ_FILES(
        params.core
    )

    fastqRegEx = PublicCore.fastqRegEx()

    params.core.sampleIds          = PublicCore.getSampleIds(params.core.assay, params.core.fastqFiles,
            fastqRegEx, log, messages)
    params.core.outputDir = PublicCore.validateOutput(params, params.core, log)

    params.core.resultsDir = "${params.core.outputDir}/Sample_Files/"
    params.core.reportsDir = "${params.core.outputDir}/report/"
    PublicCore.loadSampleLevelParams(params.core)

    if (workflow.containerEngine == 'singularity') {
        PULL_SINGULARITY_IMAGES(
            all_modules
        )
    }

    ch_images_pulled = workflow.containerEngine == 'singularity' ?
        PULL_SINGULARITY_IMAGES.out.images_pulled : Channel.of(true).first()

    ch_core_fastq      = Channel.empty()
    // Set workflow channels
    ch_core_fastq  = Channel.fromFilePairs("${params.core.input}/" + fastqRegEx, flat: true)
        .map { prefix, file1, file2 -> tuple(prefix - fastqLaneRegEx, file1, file2) }
        .groupTuple()
        .filter { !(it =~ /Undetermined/) }
        // groovylint-disable
        .ifEmpty {
            log.error("ERROR: [$params.core.assay] FASTQ read files are present in the input directory, "
                + "but the name format is incorrect. Check parameters and/or see " 
                + "documentation for file naming guidelines.")
            exit(1)
        }
        // groovylint-enable

    params.core.paramsFile = PublicCore.validateParamsFile(params)
    ch_core_params_file = params.core.paramsFile == null ?
        Channel.empty() : Channel.from(params.core.paramsFile)

    // Execute workflow
    CORE_ANALYSIS(
        ch_core_fastq,
        val_assay,
        ch_core_params_file,
        ch_images_pulled
    )

    emit:
    images_pulled           = ch_images_pulled
    input_params_yaml       = CORE_ANALYSIS.out.input_params_yaml
    params_yaml             = CORE_ANALYSIS.out.params_yaml
    command_txt             = CORE_ANALYSIS.out.command_txt
    fastqc_zip              = CORE_ANALYSIS.out.fastqc_zip
    fastqc_count            = CORE_ANALYSIS.out.fastqc_count
    fastqc_sequence_traces  = CORE_ANALYSIS.out.fastqc_sequence_traces
    fastqc_quality_scores   = CORE_ANALYSIS.out.fastqc_quality_scores
    val_pass                = CORE_ANALYSIS.out.val_pass
    complete_fastq          = CORE_ANALYSIS.out.complete_fastq
    cutadapt_metrics        = CORE_ANALYSIS.out.cutadapt_metrics
    prepared_antibody_file  = CORE_ANALYSIS.out.prepared_antibody_file
    debarcoder_fastq        = CORE_ANALYSIS.out.debarcoder_fastq
    debarcoder_edges        = CORE_ANALYSIS.out.debarcoder_edges
    debarcoder_count        = CORE_ANALYSIS.out.debarcoder_count
    debarcoder_adt_fastq    = CORE_ANALYSIS.out.debarcoder_adt_fastq
    bead_list               = CORE_ANALYSIS.out.bead_list
    ti_csv                  = CORE_ANALYSIS.out.ti_csv
}
