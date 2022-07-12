/*
Generate output reports
*/

params.options = [:]

process GENERATE_REPORT {
    container "bioraddbg/omnition-report:${workflow.manifest.version}"
    if (workflow.profile == 'aws') {
        label 'xlarge'
    } else {
        label 'cpu_small'
        label 'memory_small'
    }

    publishDir "${params.reportsDir}", mode: 'copy', overwrite: true

    input:

    path metric_summary
    path pipeline_summary
    path fastqc
    path debarcoding
    path alignments
    path qc_stats
    path determine_barcode_merges_params
    path determine_barcode_merges_implicated_barcodes
    path deconv_barcode_quant
    path deconv_barcode_allowlist
    path cell_data
    path tss_metrics
    path insert_size
    path sequence_saturation
    path crosstalk
    path bead_filt_index_summary
    path bead_filt_sample_summary
    path cluster_info
    path metric_summary
    path ti_metrics
    path ti_map
    path messages
    path sample_map
    val images_pulled

    output:
    path('atac_*.html')

    script:
    """
    mkdir tmp
    mkdir -p results/sampleId/plots
    echo ${params} | sed 's/,/\\n/g' | sed 's/:/,/g' | tr -d '[' | tr -d ']' | tail -n +3 > ./tmp/params.csv

    # due to errors with file permissions in singularity copy out report generation dir
    cp -r /opt/biorad/app/* ./

    # The report generator looks for validated.csv to check for TIs, we dont generate that file for superloaded data
    # So I am making that file here when I see that it is superloaded (presence of fastqTIreadcounts_superloaded.csv)
    if [ -f "fastqTIreadcounts_superloaded.csv" ]; then
        awk 'BEGIN{FS=OFS=","} {print \$1,\$2,\$3,\$4}' fastqTIreadcounts_superloaded.csv >> \
        SuperloadedSample.validated.csv
    fi

    python processing/process_data.py -i . -o app/assets/data/
    npm run build -- --output-directory=. --filename=atac_\$(date "+%y%m%d-%H%M").html
    """
}
