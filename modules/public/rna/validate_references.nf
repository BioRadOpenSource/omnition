/*
Checking reference formats
*/

process VALIDATE_REFERENCES {
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_small'

    input:
    path fasta
    path gtf
    val images_pulled

    output:
    path fasta, emit: fasta
    path gtf, emit: gtf

    script:
    """
    # Extracting FASTA contig names
    awk '/^>/ { print substr(\$1, 2, length(\$1)) }' ${fasta} | sort -T '.' > fasta_contigs.txt

    # Extracting GTF contig names
    awk '!/^#/ { print \$1 }' ${gtf} | sort -u -T '.' > gtf_contigs.txt

    # Verifying that all of the contigs in gtf_contigs.txt are present in fasta_contigs.txt
    common_contigs=\$(comm -13 <(sort -T '.' fasta_contigs.txt) <(sort -T '.' gtf_contigs.txt) | wc -l)
    if [[ \$common_contigs -ne 0 ]]; then
        echo "[ERROR] GTF file contains contigs not present in the reference FASTA file."
        exit 1
    fi
    """

    stub:
    """
    # Create expected output files
    outputs=()

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("fasta" "workflow.manifest.version" "gtf")

    parameters=("${fasta}" "${workflow.manifest.version}" "${gtf}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "VALIDATE_REFERENCES" "\$param_name" "\$parameter"
    done
    """
}
