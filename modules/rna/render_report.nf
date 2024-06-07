/*
Generate customer report
*/

params.options = [:]

process RENDER_REPORT {
    container "bioraddbg/omnition-report:${workflow.manifest.version}"
    publishDir "${params.options.reportsDir}", mode: 'copy', overwrite: true
    label 'cpu_small'
    label 'memory_medium'

    input:
    path samples
    path settings
    path input
    val images_pulled

    output:
    path("${params.options.prefix}rna_*.html")

    script:
    """
    mkdir tmp
    mkdir -p results/sampleId/plots
    echo "${params}" | sed 's/,/\\n/g' | sed 's/:/,/g' | tr -d '[' | tr -d ']' | tail -n +3 > ./tmp/params.csv

    cat messages.txt *_messages.txt |
    grep -v '^\$' > "messages.temp" && mv "messages.temp" "messages.txt"

    # If using Singularity, copy report files to the working directory
    if [[ "${workflow.profile}" =~ "standard" ]]
    then
        cp -r /opt/biorad/app/* ./
        APP_PREFIX="."
    else
        APP_PREFIX="/opt/biorad/app"
    fi

    # Creating the file that contains the information necessary to display release version
    printf %s "${workflow.manifest.version}" >> version.txt

    # Perform preprocessing of data
    python \${APP_PREFIX}/processing/process_data.py -i . -o \${APP_PREFIX}/app/assets/data/ -t rna

    # If using Docker or Tower, cd to the app directory
    if [[ ${workflow.profile} =~ (docker|tower|awsbatch) ]]
    then
        BASE_DIR=\$PWD
        cd \${APP_PREFIX}
    fi

    # Build the report
    npm run build:rna -- --output-directory=. --filename=${params.options.prefix}rna_\$(date "+%y%m%d-%H%M").html

    # If using Docker or Tower, move the report to the working directory
    if [[ ${workflow.profile} =~ (docker|tower|awsbatch) ]]
    then
        mv \${APP_PREFIX}/*.html \${BASE_DIR}
    fi
    """
}
