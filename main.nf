#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

log.info(PublicCore.bioradLogo(params.monochrome_logs))
messages = []
fastqLaneRegEx = ~/(_L[0-9][0-9][0-9])$/

// Check if help flag was given
if (params.help) {
    message = PublicCore.helpMessage(params.monochrome_logs)
    log.info(message)
    System.exit(0)
}

validateParameters()
log.info(paramsSummaryLog(workflow))

// Import workflows
include { RNA                } from "${params.omnition.workflowsDir}/public/rna/rna.nf"                               addParams(options: params.rna)
include { ATAC_NORMAL        } from "${params.omnition.workflowsDir}/public/atac_normal/atac_normal.nf"               addParams(options: params.atac)
include { ATAC_COMBINATORIAL } from "${params.omnition.workflowsDir}/public/atac_combinatorial/atac_combinatorial.nf" addParams(options: params.catac)
include { CORE               } from "${params.omnition.workflowsDir}/public/core/core.nf"                             addParams(options: params.core)
include { DO_MERGING         } from "${params.omnition.workflowsDir}/public/do_merging/do_merging.nf"                 addParams(options: params.do_merging)
include { CYTOMETRY          } from "${params.omnition.workflowsDir}/public/cytometry/cytometry.nf"                   addParams(options: params.cytometry)
include { CITESEQ            } from "${params.omnition.workflowsDir}/public/citeseq/citeseq.nf"                       addParams(options: params.cytometry)

workflow {
    core_modules = Channel.fromPath("${projectDir}/modules/public/core/*.nf")
    rna_modules = Channel.fromPath("${projectDir}/modules/public/rna/*.nf")
    atac_modules = Channel.fromPath("${projectDir}/modules/public/atac/*.nf")
    cytometry_modules = Channel.fromPath("${projectDir}/modules/public/cytometry/*.nf")
    citeseq_modules = Channel.fromPath("${projectDir}/modules/public/citeseq/*.nf")

    // Analyze 3' RNA Data
    if ((params.rna && !params.cytometry) || (params.cytometry && !params.rna)) {
        // Load presets for do_merging
        PublicCore.loadDOMergingPresetParams(params)

        if (params.rna) {
            // Load RNA presets
            PublicCore.loadAssayPresetParams(params, params.rna, params.preset.rna, "rna")
            all_modules = core_modules.mix(rna_modules)

            CORE(
                all_modules,
                "RNA",
                fastqLaneRegEx,
                messages
            )

            DO_MERGING(
                CORE.out.images_pulled,
                CORE.out.val_pass,
                CORE.out.debarcoder_fastq,
                CORE.out.debarcoder_edges,
                CORE.out.bead_list,
                messages
            )

            RNA(
                CORE.out.images_pulled,
                CORE.out.val_pass,
                CORE.out.input_params_yaml,
                CORE.out.params_yaml,
                CORE.out.command_txt,
                CORE.out.fastqc_zip,
                CORE.out.fastqc_count,
                CORE.out.fastqc_sequence_traces,
                CORE.out.fastqc_quality_scores,
                CORE.out.complete_fastq,
                CORE.out.cutadapt_metrics,
                CORE.out.debarcoder_fastq,
                CORE.out.debarcoder_count,
                DO_MERGING.out.edges_input,
                DO_MERGING.out.corrected_edges,
                DO_MERGING.out.barcode_translate,
                messages
            )
        }
        if (params.cytometry) {
            // Load Cytometry presets
            PublicCore.loadAssayPresetParams(params, params.cytometry, params.preset.cytometry, "cytometry")
            all_modules = core_modules
                .mix(cytometry_modules)
                .mix(rna_modules)

            CORE(
                all_modules,
                "Cytometry",
                fastqLaneRegEx,
                messages
            )

            DO_MERGING(
                CORE.out.images_pulled,
                CORE.out.val_pass,
                CORE.out.debarcoder_fastq,
                CORE.out.debarcoder_edges,
                CORE.out.bead_list,
                messages
            )

            // Analyze Cytometry data
            CYTOMETRY(
                CORE.out.images_pulled,
                CORE.out.input_params_yaml,
                CORE.out.params_yaml,
                CORE.out.command_txt,
                CORE.out.fastqc_zip,
                CORE.out.fastqc_count,
                CORE.out.complete_fastq,
                CORE.out.prepared_antibody_file,
                CORE.out.debarcoder_fastq,
                CORE.out.debarcoder_count,
                DO_MERGING.out.barcode_translate,
                messages
            )
        }
    }

    // Analyze CITE-seq Data (combined RNA + Cytometry)
    if (params.rna && params.cytometry) {
        // Load presets for do_merging
        PublicCore.loadDOMergingPresetParams(params)

        // Load RNA and cytometry presets
        PublicCore.loadAssayPresetParams(params, params.rna, params.preset.rna, "rna")
        PublicCore.loadAssayPresetParams(params, params.cytometry, params.preset.cytometry, "cytometry")
        all_modules = core_modules
            .mix(rna_modules)
            .mix(cytometry_modules)
            .mix(citeseq_modules)

        CORE(
            all_modules,
            "RNA", // groovylint-disable-line
            fastqLaneRegEx,
            messages
        )

        DO_MERGING(
            CORE.out.images_pulled,
            CORE.out.val_pass,
            CORE.out.debarcoder_fastq,
            CORE.out.debarcoder_edges,
            CORE.out.bead_list,
            messages
        )

        RNA(
            CORE.out.images_pulled,
            CORE.out.val_pass,
            CORE.out.input_params_yaml,
            CORE.out.params_yaml,
            CORE.out.command_txt,
            CORE.out.fastqc_zip,
            CORE.out.fastqc_count,
            CORE.out.fastqc_sequence_traces,
            CORE.out.fastqc_quality_scores,
            CORE.out.complete_fastq,
            CORE.out.cutadapt_metrics,
            CORE.out.debarcoder_fastq,
            CORE.out.debarcoder_count,
            DO_MERGING.out.edges_input,
            DO_MERGING.out.corrected_edges,
            DO_MERGING.out.barcode_translate,
            messages
        )

        CITESEQ(
            CORE.out.images_pulled,
            CORE.out.prepared_antibody_file,
            CORE.out.debarcoder_adt_fastq, // Pass RNA filtered fastq to avoid duplicating counts
            CORE.out.params_yaml,
            RNA.out.cell_calling_results, // Pass RNA workflow outputs to use RNA cell calling results
            RNA.out.allowlist, // Pass RNA workflow outputs to use RNA allowlist
            RNA.out.translate, // Pass RNA workflow outputs to use RNA translation
            RNA.out.filtered_mtx,
            CORE.out.debarcoder_count,
            CORE.out.fastqc_zip,
            CORE.out.fastqc_count,
            RNA.out.batch_h5ad,
            messages
        )
    }

    // Analyze ATAC data
    if (params.atac) {
        // Load ATAC presets
        PublicCore.loadAssayPresetParams(params, params.atac, params.preset.atac, "atac")
        all_modules = core_modules.mix(atac_modules)

        CORE(
            all_modules,
            "ATAC",
            fastqLaneRegEx,
            messages
        )

        ATAC_NORMAL(
            CORE.out.images_pulled,
            CORE.out.input_params_yaml,
            CORE.out.params_yaml,
            CORE.out.command_txt,
            CORE.out.fastqc_zip,
            CORE.out.fastqc_count,
            CORE.out.fastqc_sequence_traces,
            CORE.out.fastqc_quality_scores,
            CORE.out.complete_fastq,
            CORE.out.debarcoder_fastq,
            CORE.out.debarcoder_count,
            messages
        )
    }

    // Analyze cATAC data
    if (params.catac) {
        // Load cATAC presets
        PublicCore.loadAssayPresetParams(params, params.catac, params.preset.atac, "catac")

        // Prepare the TI meta map
        params.catac.tiMetaMap = [:]
        tis = [:]
        // Check if a TI file was provided
        if (params.catac.tiSheet != null) {
            params.catac.tiMetaMap = samplesheetToList("${params.catac.tiSheet}", "assets/public/sample_sheet_schema_input_ti.json")

            // Remove the empty TIs from the tiMetaMap
            for (int i = 0; i < params.catac.tiMetaMap.size(); i++) {
                currentMap = params.catac.tiMetaMap[i][0]
                currentMap = currentMap.findAll { key, value -> !(value instanceof List && value.isEmpty()) } // groovylint-disable-line
                params.catac.tiMetaMap[i][0] = currentMap
                // add the TIs to the tis map
                tis_to_add = currentMap.keySet() - ['sampleId']
                for (ti in tis_to_add) {
                    // grab the TI and the TI.sequence from the current map
                    tis[ti] = tis[ti] ?: [:]
                    tis[ti].sequence = currentMap[ti].sequence
                }
            }
            // Set the params.ti to tis
            params.catac.ti = tis
        }
        all_modules = core_modules.mix(atac_modules)

        CORE(
            all_modules,
            "cATAC",
            fastqLaneRegEx,
            messages
        )

        ATAC_COMBINATORIAL(
            CORE.out.images_pulled,
            CORE.out.input_params_yaml,
            CORE.out.params_yaml,
            CORE.out.command_txt,
            CORE.out.fastqc_zip,
            CORE.out.fastqc_count,
            CORE.out.fastqc_sequence_traces,
            CORE.out.fastqc_quality_scores,
            CORE.out.debarcoder_fastq,
            CORE.out.debarcoder_count,
            CORE.out.ti_csv,
            messages
        )
    }
}
