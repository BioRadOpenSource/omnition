/*
Validate parameters and execute rna droplet workflows
*/

include { samplesheetToList } from 'plugin/nf-schema'

// Import workflows
include { RNA_ANALYSIS            } from "${params.omnition.workflowsDir}/public/rna/rna_analysis.nf"
include { RNA_REFERENCE           } from "${params.omnition.workflowsDir}/public/rna/rna_reference.nf"

workflow RNA {
    take:
    ch_images_pulled // boolean: images pulled indicator
    ch_core_val_pass // path: validation pass file
    ch_core_input_params_yaml // path: input params yaml file
    ch_core_params_yaml // path: params yaml file
    ch_core_command_txt // path: command txt file
    ch_core_fastqc_zip // path: fastqc zip file
    ch_core_fastqc_count // path: fastqc count file
    ch_core_fastqc_sequence_traces // path: fastqc sequence traces file
    ch_core_fastqc_quality_scores // path: fastqc quality scores file
    ch_core_fastq // tuple: [ sampleId, raw R1 FASTQ files, raw R2 FASTQ files ]
    ch_core_cutadapt_metrics // tuple: [ sampleId, cutadapt metrics ]
    ch_core_debarcoded_fastq // tuple: [ sampleId, debarcoded R1 FASTQ files, debarcoded R2 FASTQ files ]
    ch_core_debarcoded_count // tuple: [ sampleId, debarcoded counts ]
    ch_do_merging_edges_input // tuple: [ sampleId, merged edges ]
    ch_do_merging_corrected_edges // tuple: [ sampleId, corrected edges ]
    ch_do_merging_barcode_translate // tuple: [ sampleId, barcode translate ]
    messages

    main:
    // Initialize mapped parameters if they were not set
    params.rna.reference   = params.rna.reference ?: [:]

    // Set global assay params
    params.rna.assay                     = "RNA"

    params.rna.sampleSheetMetaMap = [:]
    // Check if a sample sheet was provided
    if (params.rna.sampleSheet != null) {
        params.rna.sampleSheetMetaMap = samplesheetToList("${params.rna.sampleSheet}", "assets/public/sample_sheet_schema_input_rna.json")
    }

    validator = new PublicValidate(workflow, params, params.rna, params.core, log, messages, "rna")
    validator.run()

    // Validate reference workflow
    if (params.rna.workflow in [ 'reference', 'full' ]) {
        log.info("INFO: [$params.rna.assay] Executing reference workflow.")
        messages.add("INFO: [$params.rna.assay] Executing reference workflow.")
    }

    // Validate analysis workflow
    if (params.rna.workflow in [ 'analysis', 'full' ]) {
        log.info("INFO: [$params.rna.assay] Executing analysis workflow.")
        messages.add("INFO: [$params.rna.assay] Executing analysis workflow.")

        validator.runAnalysisValidation()
    }

    // Create a file in the results direcotry of all messages passed during validation
    ch_rna_messages  = Channel.fromList(messages).collectFile(
        name: "${params.rna.reportsDir}/messages.txt", newLine: true, sort: 'index')

    // Generate references
    if (params.rna.workflow in [ 'reference', 'full' ]) {
        // Set workflow channels
        // if params.rna.reference.fasta can be a comma or space separated string, split on either
        ch_rna_reference_fasta   = Channel.fromList(params.rna.reference.fasta.split(/[\s,]+/) as List) // groovylint-disable-line
        ch_rna_reference_gtf     = Channel.fromList(params.rna.reference.gtf.split(/[\s,]+/) as List) // groovylint-disable-line

        // Execute workflow
        RNA_REFERENCE(
            ch_rna_reference_fasta,
            ch_rna_reference_gtf,
            ch_images_pulled
        )
    }

    params.rna.prefix             = params.prefix != null ? params.prefix.concat('-') : ''
    params.rna.multiqc            = [ report: params.rna.prefix ?  // groovylint-disable-line
            "${params.prefix}-omnition-qc_report.html" : 'omnition-qc_report.html',  // groovylint-disable-line
            config:'{top_modules: ["fastqc","fastq_screen","cutadapt","star","picard"], fastqscreen_simpleplot: true}' ]

    // Create empty channels to emit if not running analysis
    ch_cell_calling_results = Channel.empty()
    ch_allowlist           = Channel.empty()
    ch_translate           = Channel.empty()
    ch_batch_h5ad          = Channel.empty()

    // Analyze data
    if (params.rna.workflow in [ 'analysis', 'full' ]) {
        // Set workflow channels
        ch_rna_reference_index         = params.rna.workflow == 'full' ?
            RNA_REFERENCE.out.reference_index         :
            Channel.fromPath("${params.rna.reference.directory}/star-index").first()
        ch_rna_reference_saf           = params.rna.workflow == 'full' ?
            RNA_REFERENCE.out.reference_saf           :
            Channel.fromPath("${params.rna.reference.directory}/annotation.saf").first()
        ch_rna_reference_symbols       = params.rna.workflow == 'full' ?
            RNA_REFERENCE.out.reference_symbols       :
            Channel.fromPath("${params.rna.reference.directory}/gene_symbols.txt").first()
        ch_rna_reference_refflat       = params.rna.workflow == 'full' ?
            RNA_REFERENCE.out.reference_refflat       :
            Channel.fromPath("${params.rna.reference.directory}/annotation.refflat").first()
        ch_rna_reference_interval_list = params.rna.workflow == 'full' ?
            RNA_REFERENCE.out.reference_interval_list :
            Channel.fromPath("${params.rna.reference.directory}/ribosomal.interval_list").first()

        if (params.rna.mixed) {
            ch_rna_mixed_symbols       = params.rna.workflow == 'full' ?
                RNA_REFERENCE.out.reference_mixed_symbols :
                Channel.fromPath("${params.rna.reference.directory}/*_gene_symbols.txt")
            ch_rna_reference_fasta = params.rna.workflow == 'full' ?
                RNA_REFERENCE.out.reference_fasta               :
                Channel.fromPath("${params.rna.reference.directory}/mixed.fa").first()
            ch_rna_reference_gtf   = params.rna.workflow == 'full' ?
                RNA_REFERENCE.out.reference_gtf                 :
                Channel.fromPath("${params.rna.reference.directory}/mixed.filtered.gtf").first()
        } else {
            ch_rna_mixed_symbols       = Channel.empty()
        }

        // Catch missing reference files
        Channel.empty()
            .mix(ch_rna_reference_index,
                ch_rna_reference_saf,
                ch_rna_reference_symbols,
                ch_rna_reference_refflat,
                ch_rna_reference_interval_list,
                ch_rna_mixed_symbols)
            .flatten()
            .map { item -> [ item, item.exists() ] }
            .filter { !it[1] }
            .collect { it[0] }
            .subscribe {
                log.error("ERROR: [RNA] Missing RNA reference file(s). Check \
                parameters and/or run reference workflow.\n  " + it.join('\n  '))
                exit(1)
            }

        // Create channel for sample sheet
        ch_rna_sample_sheet_map = params.rna.sampleSheetMetaMap == [:] ?
            Channel.empty() : Channel.fromList(params.rna.sampleSheetMetaMap)

        // Execute workflow
        RNA_ANALYSIS(
            ch_core_val_pass,
            ch_core_input_params_yaml,
            ch_core_params_yaml,
            ch_core_command_txt,
            ch_core_fastqc_zip,
            ch_core_fastqc_count,
            ch_core_fastqc_sequence_traces,
            ch_core_fastqc_quality_scores,
            ch_core_fastq,
            ch_core_cutadapt_metrics,
            ch_core_debarcoded_fastq,
            ch_core_debarcoded_count,
            ch_do_merging_edges_input,
            ch_do_merging_corrected_edges,
            ch_do_merging_barcode_translate,
            ch_rna_reference_index,
            ch_rna_reference_saf,
            ch_rna_reference_symbols,
            ch_rna_reference_refflat,
            ch_rna_reference_interval_list,
            ch_rna_mixed_symbols,
            params.rna.assay,
            ch_rna_sample_sheet_map,
            ch_rna_messages,
            ch_images_pulled
        )
        ch_cell_calling_results = RNA_ANALYSIS.out.cell_calling_results
        ch_allowlist = RNA_ANALYSIS.out.allowlist
        ch_translate = RNA_ANALYSIS.out.translate
        ch_batch_h5ad = RNA_ANALYSIS.out.batch_h5ad
        ch_filtered_mtx = RNA_ANALYSIS.out.filtered_matrix
    }

    emit:
    cell_calling_results = ch_cell_calling_results
    allowlist            = ch_allowlist
    translate            = ch_translate
    batch_h5ad           = ch_batch_h5ad
    filtered_mtx         = ch_filtered_mtx
}
