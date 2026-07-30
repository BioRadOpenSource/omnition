/*
Validate parameters and execute atac normal workflows
*/

include { samplesheetToList } from 'plugin/nf-schema'

// Import workflows
include { GET_BLOCKLIST_FILES     } from "${params.omnition.workflowsDir}/public/core/get_blocklist_files.nf"
include { ATAC_NORMAL_ANALYSIS    } from "${params.omnition.workflowsDir}/public/atac_normal/atac_normal_analysis.nf"
include { ATAC_NORMAL_REFERENCE   } from "${params.omnition.workflowsDir}/public/atac_normal/atac_normal_reference.nf"

workflow ATAC_NORMAL {
    take:
    ch_images_pulled           // boolean: images pulled indicator
    ch_core_input_params_yaml  // path: input params yaml file
    ch_core_params_yaml        // path: params yaml file
    ch_core_command_txt        // path: command txt file
    ch_core_fastqc_zip        // tuple: [ sampleId, fastqc zip files ]
    ch_core_fastqc_count      // tuple: [ sampleId, fastqc per base sequence content ]
    ch_core_fastqc_sequence_traces // tuple: [ sampleId, fastqc sequence traces ]
    ch_core_fastqc_quality_scores  // tuple: [ sampleId, fastqc quality scores ]
    ch_core_complete_fastq          // tuple: [ sampleId, complete R1 FASTQ files, complete R2 FASTQ files ]
    ch_core_debarcoded_fastq         // tuple: [ sampleId, debarcoded R1 FASTQ files, debarcoded R2 FASTQ files ]
    ch_core_debarcoder_count        // tuple: [ sampleId, debarcoder counts ]
    messages

    main:
    // Initialize mapped parameters if they were not set
    params.atac.reference   = params.atac.reference ?: [:]

    // Set global assay params
    params.atac.assay                     = 'ATAC'

    params.atac.sampleSheetMetaMap = [:]
    // Check if a sample sheet was provided
    if (params.atac.sampleSheet != null) {
        params.atac.sampleSheetMetaMap = samplesheetToList("${params.atac.sampleSheet}", "assets/public/sample_sheet_schema_input_atac.json")
    }

    GET_BLOCKLIST_FILES(
        params.atac
    )

    validator = new PublicValidate(workflow, params, params.atac, params.core, log, messages, "atac")
    validator.run()

    // Validate reference workflow
    if (params.atac.workflow in [ 'reference', 'full' ]) {
        log.info("INFO: [$params.atac.assay] Executing reference workflow.")
        messages.add("INFO: [$params.atac.assay] Executing reference workflow.")
    }

    // Validate analysis workflow
    if (params.atac.workflow in [ 'analysis', 'full' ]) {
        log.info("INFO: [$params.atac.assay] Executing analysis workflow.")
        messages.add("INFO: [$params.atac.assay] Executing analysis workflow.")
        validator.runAnalysisValidation()
    }

    // Create a file in the results direcotry of all messages passed during validation
    ch_atac_messages  = Channel.fromList(messages).collectFile(
        name: "${params.atac.reportsDir}/messages.txt", newLine: true, sort: 'index')

    // Generate references
    if (params.atac.workflow in [ 'reference', 'full' ]) {
        // Set workflow channels
        ch_atac_reference_fasta = Channel.fromPath(params.atac.reference.fasta.split(/[\s,]+/) as List) // groovylint-disable-line
        ch_atac_reference_gtf   = Channel.fromPath(params.atac.reference.gtf.split(/[\s,]+/) as List) // groovylint-disable-line

        // Execute reference workflow
        ATAC_NORMAL_REFERENCE(
            ch_atac_reference_fasta,
            ch_atac_reference_gtf,
            ch_images_pulled
        )
    }

    params.atac.prefix = params.prefix != null ? params.prefix.concat('-') : ''

    // Analyze data
    if (params.atac.workflow in [ 'analysis', 'full' ]) {
        ch_atac_reference_fasta = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.fasta  :
            Channel.fromPath(params.atac.reference.fasta.split(/[\s,]+/) as List) // groovylint-disable-line
        ch_atac_reference_gtf   = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.gtf    :
            Channel.fromPath(params.atac.reference.gtf.split(/[\s,]+/) as List) // groovylint-disable-line
        ch_atac_reference_blocklist = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.blocklist  :
            Channel.fromPath(params.atac.reference.blocklist.split(/[\s,]+/) as List) // groovylint-disable-line
        ch_atac_reference_index     = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.bwaIndex  :
            Channel.fromPath("${params.atac.reference.directory}/bwa-index").first()
        ch_atac_reference_sizes     = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.genomeSize :
            Channel.fromPath("${params.atac.reference.directory}/genome.sizes").first()
        /* groovylint-disable */
        ch_atac_archr_ref           = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.archrRef   :
            Channel.fromPath("${params.atac.reference.directory}/archr/BSgenome.ref.na.1.0_1.0.tar.gz")
        /* groovylint-enable */
        ch_atac_tss                 = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.tss        :
            Channel.fromPath("${params.atac.reference.directory}/TSS.${params.atac.tssWindowSize}.bed")

        // Catch missing reference files
        Channel.empty()
            .mix(ch_atac_reference_fasta,
                ch_atac_reference_gtf,
                ch_atac_reference_index,
                ch_atac_reference_sizes,
                ch_atac_archr_ref,
                ch_atac_reference_blocklist,
                ch_atac_tss)
            .flatten()
            .map { item -> [ item, item.exists() ] }
            .filter { !it[1] }
            .collect { it[0] }
            .subscribe {
                log.error('ERROR: [ATAC] Missing ATAC reference file(s). Check parameters and/or run ' + // groovylint-disable-line
                    'reference workflow.\n  ' + it.join('\n  ')) // groovylint-disable-line
                exit(1)
            }

        // Create channel for sample sheet
        ch_atac_sample_sheet_map = params.atac.sampleSheetMetaMap == [:] ?
            Channel.empty() : Channel.fromList(params.atac.sampleSheetMetaMap)

        // Execute workflow
        ATAC_NORMAL_ANALYSIS(
            ch_core_input_params_yaml,
            ch_core_params_yaml,
            ch_core_command_txt,
            ch_core_fastqc_zip,
            ch_core_fastqc_count,
            ch_core_fastqc_sequence_traces,
            ch_core_fastqc_quality_scores,
            ch_core_complete_fastq,
            ch_core_debarcoded_fastq,
            ch_core_debarcoder_count,
            ch_atac_reference_fasta,
            ch_atac_reference_gtf,
            ch_atac_reference_index,
            ch_atac_reference_sizes,
            ch_atac_archr_ref,
            ch_atac_reference_blocklist,
            ch_atac_tss,
            ch_atac_messages,
            params.atac.assay,
            ch_atac_sample_sheet_map,
            ch_images_pulled
        )
    }
}
