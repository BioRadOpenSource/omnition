/*
Reannotate bam files based on merging scheme
*/

params.options = [:]

process REANNOTATE_BAM {
    tag "${sampleId}, ${chr}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_xxsmall'

    input:
    tuple val(sampleId), val(chr), path(bam), path(index), path(bedpe), path(dict)
    val images_pulled

    output:
    tuple val(sampleId), val(chr), path("${sampleId}.*.drop.bam"),
        path("${sampleId}.*.drop.bam.bai"), emit: reannotate_bam
    tuple val(sampleId), path('*reannotate_bam_read_counts.csv'), emit: count

    script:
    """
    BEAD_TAG=XB
    DROP_TAG=XC
    atacPaintDropBarcode.py \
        --input ${bam} \
        --output ${sampleId}.${chr}.drop.bam \
        --bead-barcode \$BEAD_TAG \
        --drop-barcode \$DROP_TAG \
        --dict-file ${dict} \
        --hq-frags ${bedpe}
    samtools index -@ ${task.cpus} ${sampleId}.${chr}.drop.bam

    # Count input Proper Pairs
    input_proper_pairs=\$(samtools view -c -f 3 -F 4 -F 256 ${bam[0]})
    in_count=\$((\${input_proper_pairs} / 2))

    # Count output proper pairs
    output_proper_pairs=\$(samtools view -c -f 3 -F 4 -F 256 ${sampleId}.${chr}.drop.bam)
    out_count=\$((\${output_proper_pairs} / 2))

    # Write number of output reads to file
    READCOUNTFILE="${sampleId}_${chr}_reannotate_bam_read_counts.csv"
    printf "sample,process,metric,count\n" > \$READCOUNTFILE
    printf "${sampleId},reannotate,input,\${in_count}\n" >> \$READCOUNTFILE
    printf "${sampleId},reannotate,output,\${out_count}\n" >> \$READCOUNTFILE
    """
}
