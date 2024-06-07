/*
Compile qc stats for TI samples
*/

params.options = [:]

process COMPILE_QC_STATS {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/alignments",
        pattern: '*{.dedup_stats.txt,.barcodeTranslate.tsv}', mode: 'copy',
        overwrite: true
    label 'cpu_large'
    label 'memory_xxsmall'

    input:
    tuple val(sampleId), path(stats), path(dict), path(dedup), path(cells)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}.QCstats.csv"), emit:stats
    tuple val(sampleId), path("${sampleId}.barcodeTranslate.tsv"), emit:translate
    tuple val(sampleId), path("${sampleId}.dedup_stats.txt"), emit:dedup
    tuple val(sampleId), path("${sampleId}.cell_data.csv"), emit:cells

    script:
    """
    # Take the header of the first file and add all files
    head -1 ${stats[0]} > ${sampleId}.QCstats.csv
    tail -n +2 -q ${stats} >> ${sampleId}.QCstats.csv

    # Take the header of the first file and add all files
    # Sort by third column, descending
    head -1 ${cells[0]} > ${sampleId}.cell_data.csv
    tail -n +2 -q ${cells} | sort -t, -nk3 -r -T '.' >> ${sampleId}.cell_data.csv

    cat ${dict} > ${sampleId}.barcodeTranslate.tsv

    atacCompileDedupStats.py -s ${sampleId} -i ./
    """
}
