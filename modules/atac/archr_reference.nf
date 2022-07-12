/*
Generate ArchR references
*/

params.options = [:]

process ARCHR_REFERENCE {
    container "bioraddbg/omnition-archr:${workflow.manifest.version}"
    publishDir "${params.options.reference.directory}/archr", mode: 'copy'
    if (workflow.profile == 'aws') {
        label 'xlarge'
    } else {
        label 'cpu_small'
        label 'memory_large'
    }

    input:
    path fasta
    path gtf
    val images_pulled

    output:
    path 'BSgenome.ref.na.1.0_1.0.tar.gz', emit: pkg

    script:
    """
    mkdir refchrs

    export TMPDIR=\$(pwd)/tmp
    mkdir -p \$TMPDIR

    fapath=\$(realpath ${fasta})
    (cd refchrs; faidx -x \$fapath)

    atacArchrReference.R -g ${gtf}
    """
}
