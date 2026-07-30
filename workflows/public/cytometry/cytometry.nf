/*
Validate parameters and execute cytometry workflows
*/

include { samplesheetToList } from 'plugin/nf-schema'

// Import workflows
include { CYTOMETRY_REFERENCE     } from "${params.omnition.workflowsDir}/public/cytometry/cytometry_reference.nf"
include { CYTOMETRY_ANALYSIS      } from "${params.omnition.workflowsDir}/public/cytometry/cytometry_analysis.nf"

workflow CYTOMETRY {
    take:
    ch_images_pulled // boolean: images pulled indicator
    ch_core_input_params_yaml // path: input params yaml file
    ch_core_params_yaml // path: params yaml file
    ch_core_command_txt // path: command txt file
    ch_core_fastqc_zip // path: fastqc zip file
    ch_core_fastqc_count // path: fastqc count file
    ch_core_fastq // tuple: [ sampleId, raw R1 FASTQ files, raw R2 FASTQ files ]
    ch_core_prepared_antibody_file // path: formatted antibody file
    ch_core_debarcoded_fastq // tuple: [ sampleId, debarcoded R1 FASTQ files, debarcoded R2 FASTQ files ]
    ch_core_debarcoded_count // tuple: [ sampleId, debarcoded counts ]
    ch_do_merging_barcode_translate // tuple: [ sampleId, barcode translate ]
    messages

    main:
    // Initialize mapped parameters if they were not set
    params.cytometry.reference   = params.cytometry.reference ?: [:]

    // Set global assay params
    params.cytometry.assay = "Cytometry"

    params.cytometry.sampleSheetMetaMap = [:]
    // Check if a sample sheet was provided
    if (params.cytometry.sampleSheet != null) {
        params.cytometry.sampleSheetMetaMap = samplesheetToList("${params.cytometry.sampleSheet}", "assets/public/sample_sheet_schema_input_cytometry.json")
    }

    validator = new PublicValidate(workflow, params, params.cytometry, params.core, log, messages, "cytometry")
    validator.run()

    // Validate reference workflow
    if (params.cytometry.workflow in [ 'reference', 'full' ]) {
        log.info("INFO: [$params.cytometry.assay] Executing reference workflow.")
        messages.add("INFO: [$params.cytometry.assay] Executing reference workflow.")
        log.info("INFO: [$params.cytometry.assay] Executing analysis workflow.")
        messages.add("INFO: [$params.cytometry.assay] Executing analysis workflow.")
        validator.runAnalysisValidation()
    }

    // Create a file in the results direcotry of all messages passed during validation
    ch_cytometry_messages  = Channel.fromList(messages).collectFile(
        name: "${params.cytometry.reportsDir}/messages.txt", newLine: true, sort: 'index')

    // Generate references
    if (params.cytometry.workflow in [ 'reference', 'full' ]) {
        // Execute workflow
        CYTOMETRY_REFERENCE(
            ch_core_prepared_antibody_file,
            ch_images_pulled
        )
    }

    params.cytometry.prefix = params.prefix != null ? params.prefix.concat('-') : ''

    // Analyze data
    if (params.cytometry.workflow in [ 'analysis', 'full' ]) {
        // Set workflow channels
        ch_cytometry_reference_kite_t2g          = params.cytometry.workflow == 'full' ?
            CYTOMETRY_REFERENCE.out.reference_kite_t2g        :
            Channel.fromPath("${params.cytometry.reference.directory}/features_processed.t2g").first()
        ch_cytometry_reference_kallisto_index    = params.cytometry.workflow == 'full' ?
            CYTOMETRY_REFERENCE.out.reference_kallisto_index  :
            Channel.fromPath("${params.cytometry.reference.directory}/adt_index.idx").first()
        ch_cytometry_reference_seqsearch_index   = params.cytometry.workflow == 'full' ?
            CYTOMETRY_REFERENCE.out.reference_seqsearch_index :
            Channel.fromPath("${params.cytometry.reference.directory}/seqsearch_index.idx").first()

        // Catch missing reference files
        Channel.empty()
            .mix(ch_cytometry_reference_kite_t2g,
                ch_cytometry_reference_kallisto_index,
                ch_cytometry_reference_seqsearch_index)
            .flatten()
            .map { item -> [ item, item.exists() ] }
            .filter { !it[1] }
            .collect { it[0] }
            .subscribe {
                log.error("ERROR: [CYTOMETRY] Missing CYTOMETRY reference file(s). Check \
                parameters and/or run reference workflow.\n  " + it.join('\n  '))
                exit(1)
            }

        // Create channel for sample sheet
        ch_cytometry_sample_sheet_map = params.cytometry.sampleSheetMetaMap == [:] ?
            Channel.empty() : Channel.fromList(params.cytometry.sampleSheetMetaMap)

        // Execute workflow
        CYTOMETRY_ANALYSIS(
            ch_core_input_params_yaml,
            ch_core_params_yaml,
            ch_core_command_txt,
            ch_core_fastqc_zip,
            ch_core_fastqc_count,
            ch_core_fastq,
            ch_core_debarcoded_fastq,
            ch_core_debarcoded_count,
            ch_do_merging_barcode_translate,
            ch_core_prepared_antibody_file,
            ch_cytometry_reference_kite_t2g,
            ch_cytometry_reference_kallisto_index,
            ch_cytometry_reference_seqsearch_index,
            params.cytometry.assay,
            ch_cytometry_sample_sheet_map,
            ch_cytometry_messages,
            ch_images_pulled
        )
    }
}
