/*
Run ArchR pipeline
*/

process ARCHR {
    tag "${sampleId}"
    container "${params.container_image.archr_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}", pattern:'ArchR/ArrowFiles/*.{arrow}',
        mode: 'copy', overwrite: true
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}", pattern:'ArchR/*.gz',
        mode: 'copy', overwrite: true
    label 'cpu_small'
    label 'memory_xlarge'

    input:
    tuple val(sampleId), path(bam), path(bai), path(clean_peaks), path(gtf), path(archr_ref), path(blocklist)
    val images_pulled

    output:
    tuple path('ArchR/*.rds'), path('ArchR/ArrowFiles/*.arrow'), emit: archrproj
    path 'ArchR/*_umap_and_clusterID.rda', emit: clusterinfo
    tuple path('ArchR/Plots/*-peak-heatmap.pdf'), path('ArchR/*-peak-heatmap.tsv'), optional: true, emit: heatmap
    tuple path("ArchR/${sampleId}.mtx.gz"), path("ArchR/${sampleId}.rownames.tsv.gz"),
        path("ArchR/${sampleId}.colnames.tsv.gz"), optional: true, emit: matrix

    script:
    """
    export TMPDIR=\$(pwd)/tmp
    mkdir -p \$TMPDIR

    export R_LIBS_USER='r_libs_user/'
    mkdir \$R_LIBS_USER

    if [[ "$workflow.profile" =~ "demo" ]]
    then
        DEMO=True
    else
        DEMO=False
    fi

    R CMD INSTALL BSgenome.ref.na.1.0_1.0.tar.gz

    publicAtacArchr.R -c ${task.cpus} -g ${gtf} -b ${blocklist} -p ${clean_peaks} -d \$DEMO

    pigz -p ${task.cpus} ArchR/${sampleId}.rownames.tsv
    pigz -p ${task.cpus} ArchR/${sampleId}.colnames.tsv
    pigz -p ${task.cpus} ArchR/${sampleId}.mtx
    """
}
