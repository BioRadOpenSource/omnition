/*
Validate parameters and execute atac combinatorial workflows
*/

include { samplesheetToList } from 'plugin/nf-schema'

// Import workflows
include { GET_BLOCKLIST_FILES          } from "${params.omnition.workflowsDir}/public/core/get_blocklist_files.nf"
include { ATAC_COMBINATORIAL_ANALYSIS  } from "${params.omnition.workflowsDir}/public/atac_combinatorial/atac_combinatorial_analysis.nf"
include { ATAC_COMBINATORIAL_REFERENCE } from "${params.omnition.workflowsDir}/public/atac_combinatorial/atac_combinatorial_reference.nf"

workflow ATAC_COMBINATORIAL {
    take:
    ch_images_pulled           // boolean: images pulled indicator
    ch_core_input_params_yaml  // path: input params yaml file
    ch_core_params_yaml        // path: params yaml file
    ch_core_command_txt        // path: command txt file
    ch_core_fastqc_zip        // tuple: [ sampleId, fastqc zip files ]
    ch_core_fastqc_count      // tuple: [ sampleId, fastqc per base sequence content ]
    ch_core_fastqc_sequence_traces // tuple: [ sampleId, fastqc sequence traces ]
    ch_core_fastqc_quality_scores  // tuple: [ sampleId, fastqc quality scores ]
    ch_core_debarcoded_fastq         // tuple: [ sampleId, debarcoded R1 FASTQ files, debarcoded R2 FASTQ files ]
    ch_core_debarcoder_count        // tuple: [ sampleId, debarcoder counts ]
    ch_core_ti_csv                // CSV config file of fastq and TIs
    messages

    main:
    // Initialize mapped parameters if they were not set
    params.catac.reference   = params.catac.reference ?: [:]

    // Set global assay params
    params.catac.assay                     = 'cATAC'

    params.catac.sampleSheetMetaMap = [:]
    // Check if a sample sheet was provided
    if (params.catac.sampleSheet != null) {
        params.catac.sampleSheetMetaMap = samplesheetToList("${params.catac.sampleSheet}", "assets/public/sample_sheet_schema_input_atac.json")
    }

    GET_BLOCKLIST_FILES(
        params.catac
    )

    validator = new PublicValidate(workflow, params, params.catac, params.core, log, messages, "catac")
    validator.run()

    // Validate reference workflow
    if (params.catac.workflow in [ 'reference', 'full' ]) {
        log.info("INFO: [$params.catac.assay] Executing reference workflow.")
        messages.add("INFO: [$params.catac.assay] Executing reference workflow.")
    }

    // Validate analysis workflow
    if (params.catac.workflow in [ 'analysis', 'full' ]) {
        log.info("INFO: [$params.catac.assay] Executing analysis workflow.")
        messages.add("INFO: [$params.catac.assay] Executing analysis workflow.")

        validator.runAnalysisValidation()
    }

    // Create a file in the results direcotry of all messages passed during validation
    ch_atac_messages  = Channel.fromList(messages).collectFile(
        name: "${params.catac.reportsDir}/messages.txt", newLine: true, sort: 'index')

    // Generate references
    if (params.catac.workflow in [ 'reference', 'full' ]) {
        // Set workflow channels
        ch_atac_reference_fasta = Channel.fromPath(params.catac.reference.fasta.split(/[\s,]+/) as List) // groovylint-disable-line
        ch_atac_reference_gtf   = Channel.fromPath(params.catac.reference.gtf.split(/[\s,]+/) as List) // groovylint-disable-line

        // Execute reference workflow
        ATAC_COMBINATORIAL_REFERENCE(
            ch_atac_reference_fasta,
            ch_atac_reference_gtf,
            ch_images_pulled
        )
    }

    params.catac.prefix             = params.prefix != null ? params.prefix.concat('-') : ''

    // Analyze data
    if (params.catac.workflow in [ 'analysis', 'full' ]) {
        // Set workflow channels
        ch_atac_TI_config = params.catac.barcodedTn5Config ? Channel.fromPath(params.catac.barcodedTn5Config) :
            Channel.empty()
        ch_atac_reference_fasta = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.fasta  :
            Channel.fromPath(params.catac.reference.fasta.split(/[\s,]+/) as List) // groovylint-disable-line
        ch_atac_reference_gtf   = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.gtf    :
            Channel.fromPath(params.catac.reference.gtf.split(/[\s,]+/) as List) // groovylint-disable-line
        ch_atac_reference_blocklist = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.blocklist :
            Channel.fromPath(params.catac.reference.blocklist.split(/[\s,]+/) as List) // groovylint-disable-line
        ch_atac_reference_index     = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.bwaIndex  :
            Channel.fromPath("${params.catac.reference.directory}/bwa-index").first()
        ch_atac_reference_sizes     = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.genomeSize :
            Channel.fromPath("${params.catac.reference.directory}/genome.sizes").first()
        /* groovylint-disable */
        ch_atac_archr_ref           = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.archrRef   :
            Channel.fromPath("${params.catac.reference.directory}/archr/BSgenome.ref.na.1.0_1.0.tar.gz")
        /* groovylint-enable */
        ch_atac_tss                 = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.tss        :
            Channel.fromPath("${params.catac.reference.directory}/TSS.${params.catac.tssWindowSize}.bed")

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

        // Create channel for sample sheet and TIs
        ch_catac_sample_sheet_map = params.catac.sampleSheetMetaMap == [:] ?
            Channel.empty() : Channel.fromList(params.catac.sampleSheetMetaMap)
        ch_catac_ti_map = params.catac.tiMetaMap == [:] ?
            Channel.empty() : Channel.fromList(params.catac.tiMetaMap)

        ATAC_COMBINATORIAL_ANALYSIS(
            ch_core_input_params_yaml,
            ch_core_params_yaml,
            ch_core_command_txt,
            ch_core_fastqc_zip,
            ch_core_fastqc_count,
            ch_core_fastqc_sequence_traces,
            ch_core_fastqc_quality_scores,
            ch_core_debarcoded_fastq,
            ch_core_debarcoder_count,
            ch_core_ti_csv,
            ch_atac_reference_fasta,
            ch_atac_reference_gtf,
            ch_atac_reference_index,
            ch_atac_reference_sizes,
            ch_atac_archr_ref,
            ch_atac_reference_blocklist,
            ch_atac_tss,
            ch_atac_TI_config,
            ch_atac_messages,
            params.catac.assay,
            ch_catac_sample_sheet_map,
            ch_catac_ti_map,
            ch_images_pulled
        )
    }
}
