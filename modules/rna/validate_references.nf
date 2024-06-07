/*
Checking reference formats
*/

params.options = [:]

process VALIDATE_REFERENCES {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
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
}
