/*
Validate parameters and execute cytometry, and citeseq workflows
*/

include { samplesheetToList } from 'plugin/nf-schema'

// Import workflows
include { CYTOMETRY_REFERENCE     } from "${params.omnition.workflowsDir}/public/cytometry/cytometry_reference.nf"
include { CITESEQ_ANALYSIS        } from "${params.omnition.workflowsDir}/public/citeseq/citeseq_analysis.nf"

workflow CITESEQ {
    take:
    ch_images_pulled // boolean: images pulled indicator
    ch_core_prepared_antibody_file // path: formatted antibody file
    ch_core_debarcoded_fastq // tuple: [ sampleId, debarcoded R1 FASTQ files, raw R2 FASTQ files ]
    ch_core_params_yaml // path: params yaml file
    ch_rna_cell_calling_results // tuple: [ sampleId, csvs with cell calling results ]
    ch_rna_reference_allowlist // tuple: [ sampleId, allowlist file ]
    ch_rna_reference_translate // tuple: [ sampleId, translate file ]
    ch_rna_filtered_mtx // tuple: [ sampleId, mtx, barcodes, genes files ]
    ch_core_debarcoder_counts // tuple: [ sampleId, debarcoder counts csv ]
    ch_core_fastqc_zip // path: fastqc zip file
    ch_core_fastqc_count // path: fastqc count file
    ch_rna_batch_h5ad // path: rna batch h5ad file
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

    // Validate cytometry parameters
    cytometry_validator = new PublicValidate(workflow, params, params.cytometry, params.core, log, messages, "cytometry")
    cytometry_validator.run()

    if (params.cytometry.workflow in [ 'reference', 'full' ]) {
        log.info("INFO: [$params.cytometry.assay] Executing reference workflow.")
        messages.add("INFO: [$params.cytometry.assay] Executing reference workflow.")
    }

    // Generate cytometry references
    if (params.cytometry.workflow in [ 'reference', 'full' ]) {
        // Execute cytometry reference workflow
        CYTOMETRY_REFERENCE(
            ch_core_prepared_antibody_file,
            ch_images_pulled
        )
    }

    // Set cytometry, multiqc parameters
    params.cytometry.prefix = params.prefix != null ? params.prefix.concat('-') : ''

    // Analyze data - truncated cytometry
    if (params.cytometry.workflow in [ 'analysis', 'full' ]) {
        // Cytometry reference channels
        ch_cytometry_reference_kite_t2g                = params.cytometry.workflow == 'full' ?
            CYTOMETRY_REFERENCE.out.reference_kite_t2g        :
            Channel.fromPath("${params.cytometry.reference.directory}/features_processed.t2g").first()
        ch_cytometry_reference_kallisto_index          = params.cytometry.workflow == 'full' ?
            CYTOMETRY_REFERENCE.out.reference_kallisto_index  :
            Channel.fromPath("${params.cytometry.reference.directory}/adt_index.idx").first()
        ch_cytometry_reference_seqsearch_index         = params.cytometry.workflow == 'full' ?
            CYTOMETRY_REFERENCE.out.reference_seqsearch_index :
            Channel.fromPath("${params.cytometry.reference.directory}/seqsearch_index.idx").first()

        // Catch missing cytometry reference files
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
                parameters and/or run reference workflow.\n  " + it.join('\n  ')) // groovylint-disable-line
                exit(1)
            }

        // Create sample sheet channels
        ch_cytometry_sample_sheet_map = params.cytometry.sampleSheetMetaMap == [:] ?
            Channel.empty() : Channel.fromList(params.cytometry.sampleSheetMetaMap)

        log.info("INFO: [$params.cytometry.assay] Executing analysis workflow.")
        messages.add("INFO: [$params.cytometry.assay] Executing analysis workflow.")
        cytometry_validator.runAnalysisValidation()

        // Execute CITESEQ analysis workflow (modified cytometry that uses RNA outputs)
        CITESEQ_ANALYSIS(
            ch_core_debarcoded_fastq, // Pass RNA filtered fastq to avoid duplicating counts
            ch_core_params_yaml,
            ch_core_prepared_antibody_file,
            ch_cytometry_reference_kite_t2g,
            ch_cytometry_reference_kallisto_index,
            ch_cytometry_reference_seqsearch_index,
            ch_rna_cell_calling_results, // Pass RNA workflow outputs to use RNA cell calling results
            ch_rna_reference_allowlist, // Pass RNA workflow outputs to use RNA allowlist
            ch_rna_reference_translate, // Pass RNA workflow outputs to use RNA translation
            ch_rna_filtered_mtx,
            ch_core_debarcoder_counts,
            ch_core_fastqc_zip,
            ch_core_fastqc_count,
            ch_rna_batch_h5ad,
            params.cytometry.assay,
            ch_cytometry_sample_sheet_map,
            messages,
            ch_images_pulled
        )
    }
}
