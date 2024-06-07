/*
Collecting metadata on edges
*/

params.options = [:]

process EDGE_METADATA {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xlarge'

    input:
    tuple val(sampleId), path(raw_edges), path(corrected_edges), path(allowlist)
    val images_pulled

    output:
    tuple val(sampleId), path("*.bead_merging_metadata.csv"), emit: metadata

    script:
    """
        # Create an above knee only edge file
        rnaFilterBeads.py -e ${corrected_edges} \
            -b ${allowlist} \
            -s ${sampleId}_corrected_filtered

        # Sort the beads lexicographically
        # If its a self-dimer, sort the UMIs lexicographically
        # count the different raw edge types in a pipe
        awk -F'\t' 'BEGIN {OFS = FS} {
            if (\$1 == \$3 && \$2 >= \$4) print \$1, \$2, \$3, \$4;
            else if (\$1 == \$3 && \$2 < \$4) print \$3, \$4, \$1, \$2;
            else if (\$1>\$3) print \$1, \$2, \$3, \$4;
            else print \$3, \$4, \$1, \$2}' ${raw_edges} | tee >(wc -l > TOTAL_EDGES) |
            awk '{!seen[\$0]++}END{for (i in seen) print seen[i] "\t" i}' |
            tee >(wc -l > R_UNIQUE_EDGES) | tee >(grep -c N > R_AMBIGUOUS_EDGES) |
            grep -v N | cut -f2,4 -d\$'\t' | awk '{!seen[\$0]++}END{for (i in seen) print seen[i] "\t" i}' |
            wc -l > R_FINAL_EDGES

        # bash math to count the number of raw redundant edges
        R_REDUNDANT_EDGES=\$(( \$(cat R_UNIQUE_EDGES) - \$(cat R_AMBIGUOUS_EDGES) - \$(cat R_FINAL_EDGES)))

        # count the different corrected edge types in a pipe
        cat ${corrected_edges} |
            awk '{!seen[\$0]++}END{for (i in seen) print seen[i] "\t" i}' |
            tee >(wc -l > C_UNIQUE_EDGES) | tee >(grep -c N > C_AMBIGUOUS_EDGES) |
            grep -v N | cut -f2,4 -d\$'\t' | awk '{!seen[\$0]++}END{for (i in seen) print seen[i] "\t" i}' |
            wc -l > C_FINAL_EDGES

        # bash math to count the number of corrected redundant edges
        C_REDUNDANT_EDGES=\$(( \$(cat C_UNIQUE_EDGES) - \$(cat C_AMBIGUOUS_EDGES) - \$(cat C_FINAL_EDGES)))

        # count the different above knee edge types in a pipe
        cat ${sampleId}_corrected_filtered_edges.tsv | tee >(wc -l > ABOVE_KNEE_EDGES) |
            awk '{!seen[\$0]++}END{for (i in seen) print seen[i] "\t" i}' |
            tee >(wc -l > A_UNIQUE_EDGES) | tee >(grep -c N > A_AMBIGUOUS_EDGES) |
            grep -v N | cut -f2,4 -d\$'\t' | awk '{!seen[\$0]++}END{for (i in seen) print seen[i] "\t" i}' |
            wc -l > A_FINAL_EDGES

        # bash math to count the number of above knee redundant edges
        A_REDUNDANT_EDGES=\$(( \$(cat A_UNIQUE_EDGES) - \$(cat A_AMBIGUOUS_EDGES) - \$(cat A_FINAL_EDGES)))

        # write edge counts to metadata csv
        printf "total_edges,above_knee_edges," \
            > ${sampleId}.bead_merging_metadata.csv
        printf "raw_unique_edges,raw_ambiguous_umis,raw_redundant_edges," \
            >> ${sampleId}.bead_merging_metadata.csv
        printf "corrected_unique_edges,corrected_ambiguous_umis,corrected_redundant_edges," \
            >> ${sampleId}.bead_merging_metadata.csv
        printf "above_knee_unique_edges,above_knee_ambiguous_umis,above_knee_redundant_edges\n" \
            >> ${sampleId}.bead_merging_metadata.csv
        printf "\$(cat TOTAL_EDGES),\$(cat ABOVE_KNEE_EDGES)," \
            >> ${sampleId}.bead_merging_metadata.csv
        printf "\$(cat R_UNIQUE_EDGES),\$(cat R_AMBIGUOUS_EDGES),\$R_REDUNDANT_EDGES," \
            >> ${sampleId}.bead_merging_metadata.csv
        printf "\$(cat C_UNIQUE_EDGES),\$(cat C_AMBIGUOUS_EDGES),\$C_REDUNDANT_EDGES," \
            >> ${sampleId}.bead_merging_metadata.csv
        printf "\$(cat A_UNIQUE_EDGES),\$(cat A_AMBIGUOUS_EDGES),\$A_REDUNDANT_EDGES\n" \
            >> ${sampleId}.bead_merging_metadata.csv
    """
}
