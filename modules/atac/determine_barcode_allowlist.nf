/*
Bead filtration with knee-calling with call_nums_cells.R
*/

params.options = [:]

process DETERMINE_BARCODE_ALLOWLIST {
    tag "${sampleId}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_small'

    input:
    tuple val(sampleId), path(count)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}.barcodeQuantSimple.csv"), emit: quant
    tuple val(sampleId), path("${sampleId}_barcode_allowlist.csv"), emit: allowlist
    tuple val(sampleId), path("${sampleId}.deconvolutionParams.orig.csv"), emit: params

    script:
    if (params.options.barcode.force.get(sampleId)) {
        """
        # Combine counts
        > ${sampleId}_tempCounts.tsv
        for FILE in ${count}; do
            if [[ "\$FILE" != *"${params.options.mitoContig}"* ]];then
                cat \$FILE >> ${sampleId}_tempCounts.tsv
            fi
        done
        awk '{ s[\$1] += \$2} END { for (i in s) { print i","s[i]} }' \
            ${sampleId}_tempCounts.tsv > ${sampleId}.barcodeQuantSimple.csv

        atacCallNumCells.R \
            --count_matrix ${sampleId}.barcodeQuantSimple.csv \
            --sample_name ${sampleId} \
            --override ${params.options.barcode.force.get( sampleId )} \
            --results_folder ./

        coreSubsampleByDistance.py -i ${sampleId}_allalgos_loglog.csv -n 1000 -x 2 -y 3
        coreSubsampleByDistance.py -i ${sampleId}_allalgos_cumfrac.csv -n 1000

        # Identify knee cutoffs
        CELLCALLER=\$(< ${sampleId}_barcode_allowlist.csv wc -l)
        MINFRAG=\$(awk NR=="\$CELLCALLER + 1" ${sampleId}_allalgos_loglog.csv |  cut -d, -f2)

        > ${sampleId}.deconvolutionParams.orig.csv
        echo "bead_threshold_nosafety,\$MINFRAG" >> ${sampleId}.deconvolutionParams.orig.csv
        echo "bead_threshold,\$MINFRAG" >> ${sampleId}.deconvolutionParams.orig.csv
        echo "above_knee,\$CELLCALLER" >> ${sampleId}.deconvolutionParams.orig.csv
        """
    } else {
        """
        # Combine counts
        > ${sampleId}_tempCounts.tsv
        for FILE in ${count}; do
            if [[ "\$FILE" != *"${params.options.mitoContig}"* ]];then
                cat \$FILE >> ${sampleId}_tempCounts.tsv
            fi
        done
        awk '{ s[\$1] += \$2} END { for (i in s) { print i","s[i]} }' \
            ${sampleId}_tempCounts.tsv > ${sampleId}.barcodeQuantSimple.csv

        atacCallNumCells.R \
            --count_matrix ${sampleId}.barcodeQuantSimple.csv \
            --sample_name ${sampleId} \
            --results_folder ./

        coreSubsampleByDistance.py -i ${sampleId}_allalgos_loglog.csv -n 1000 -x 2 -y 3
        coreSubsampleByDistance.py -i ${sampleId}_allalgos_cumfrac.csv -n 1000

        # Identify knee cutoffs
        CELLCALLER=\$(< ${sampleId}_barcode_allowlist.csv wc -l)
        MINFRAG=\$(awk NR=="\$CELLCALLER + 1" ${sampleId}_allalgos_loglog.csv |  cut -d, -f2)

        > ${sampleId}.deconvolutionParams.orig.csv
        echo "bead_threshold_nosafety,\$MINFRAG" >> ${sampleId}.deconvolutionParams.orig.csv
        echo "bead_threshold,\$MINFRAG" >> ${sampleId}.deconvolutionParams.orig.csv
        echo "above_knee,\$CELLCALLER" >> ${sampleId}.deconvolutionParams.orig.csv
        """
    }
}
