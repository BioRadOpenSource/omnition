/*
Generate output reports
*/

params.options = [:]

process RENDER_REPORT {
    container "bioraddbg/omnition-report:${workflow.manifest.version}"
    publishDir "${params.options.reportsDir}", mode: 'copy', overwrite: true
    label 'cpu_small'
    label 'memory_medium'

    input:
    path input
    val images_pulled

    output:
    path("${params.options.prefix}atac_*.html")

    script:
    params.external = null
    """
    mkdir tmp
    mkdir -p results/sampleId/plots
    echo ${params} | sed 's/,/\\n/g' | sed 's/:/,/g' | tr -d '[' | tr -d ']' | tail -n +3 > ./tmp/params.csv

    # If using Singularity, copy report files to the working directory
    if [[ "$workflow.profile" =~ "standard" ]]
    then
        cp -r /opt/biorad/app/* ./
        APP_PREFIX="."
    else
        APP_PREFIX="/opt/biorad/app"
    fi

    # The report generator looks for validated.csv to check for TIs, we dont generate that file for superloaded data
    # So I am making that file here when I see that it is superloaded (presence of fastqTIreadcounts_superloaded.csv)
    if [ -f "fastqTIreadcounts_superloaded.csv" ]; then
        awk 'BEGIN{FS=OFS=","} {print \$1,\$2,\$3,\$4}' fastqTIreadcounts_superloaded.csv >> \
        SuperloadedSample.validated.csv
    fi

    # Creating the file that contains the information necessary to display release version
    printf %s "${workflow.manifest.version}" >> version.txt

    # Perform preprocessing of data
    python \${APP_PREFIX}/processing/process_data.py -i . -o \${APP_PREFIX}/app/assets/data/ -t atac

    # If using Docker or Tower, cd to the app directory
    if [[ $workflow.profile =~ (docker|tower|awsbatch) ]]
    then
        BASE_DIR=\$PWD
        cd \${APP_PREFIX}
    fi

    npm run build:atac -- --output-directory=. \
    --filename="${params.options.prefix}atac_\$(date "+%y%m%d-%H%M").html"

    # If using Docker or Tower, move the report to the working directory
    if [[ $workflow.profile =~ (docker|tower|awsbatch) ]]
    then
        mv \${APP_PREFIX}/*.html \${BASE_DIR}
    fi
    """
}
