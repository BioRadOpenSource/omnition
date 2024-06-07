/*
Validate parameters and execute rna droplet workflows
*/

// Import workflows
include { PULL_IMAGES } from '../core/singularity.nf'
include { RNA_ANALYSIS } from './rna_analysis.nf'
include { RNA_REFERENCE } from './rna_reference.nf'
include { GET_FASTQ_FILES } from '../core/fastq_files_array.nf'

workflow RNA {
    take:
    messages

    main:
    // Initialize mapped parameters if they were not set
    params.rna.reference   = params.rna.reference ?: [:]

    // Set global assay params
    params.rna.assay                     = "3' RNA Droplet"

    validator = new Validate(workflow, params, params.rna, log, messages)
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
        GET_FASTQ_FILES(
            params.rna
        )
        validator.runAnalysisValidation()
    }
    if (workflow.containerEngine == 'singularity') {
        workflow_modules = Channel.fromPath("${baseDir}/modules/rna/*.nf")
        core_modules = Channel.fromPath("${baseDir}/modules/core/*.nf")

        PULL_IMAGES(
            workflow_modules.mix(core_modules)
        )
    }

    ch_images_pulled = workflow.containerEngine == 'singularity' ?
        PULL_IMAGES.out.images_pulled : Channel.of(true).first()

    // Create a file in the results direcotry of all messages passed during validation
    ch_rna_messages  = Channel.fromList(messages).collectFile(
        name: "${params.rna.reportsDir}/messages.txt", newLine: true, sort: 'index')

    // Generate references
    if (params.rna.workflow in [ 'reference', 'full' ]) {
        // Set workflow channels
        ch_rna_reference_fasta   = Channel.fromPath(params.rna.reference.fasta)
        ch_rna_reference_gtf     = Channel.fromPath(params.rna.reference.gtf)

        // Execute workflow
        RNA_REFERENCE(
            ch_rna_reference_fasta,
            ch_rna_reference_gtf,
            ch_images_pulled
        )
    }

    // Analyze data
    if (params.rna.workflow in [ 'analysis', 'full' ]) {
        // Set workflow channels
        ch_rna_fastq  = Channel.fromFilePairs("${params.rna.input}/" + Core.fastqRegEx(), flat: true)
            .map { prefix, file1, file2 -> tuple(prefix - ~/(_L[0-9][0-9][0-9])$/, file1, file2) }
            .groupTuple()
            .filter { !(it =~ /Undetermined/) }
            .ifEmpty {
                log.error("ERROR: [$params.rna.assay] FASTQ read files are present in the input directory, " \
                    + "but the name format is incorrect. Check parameters and/or see " \
                    + "documentation for file naming guidelines.")
                exit(1)
            }

        ch_rna_params_file                     = params.rna.paramsFile == null ?
            Channel.empty()                                   :
            Channel.from(params.rna.paramsFile)
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
                log.error("ERROR: [3' RNA] Missing 3' RNA Droplet reference file(s). Check \
                parameters and/or run reference workflow.\n  " + it.join('\n  '))
                exit(1)
            }

        // Execute workflow
        RNA_ANALYSIS(
            ch_rna_fastq,
            ch_rna_params_file,
            ch_rna_reference_index,
            ch_rna_reference_saf,
            ch_rna_reference_symbols,
            ch_rna_reference_refflat,
            ch_rna_reference_interval_list,
            ch_rna_mixed_symbols,
            params.rna.assay,
            ch_rna_messages,
            ch_images_pulled
        )
    }
}
