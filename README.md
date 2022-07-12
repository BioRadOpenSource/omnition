# Bio-Rad Laboratories, Inc. Omnition Single-Cell Analysis Software

This analysis pipeline is designed to analyze single-cell ATAC-Seq data and combinatorial ATAC-Seq data. The pipeline is built using a combination of [Nextflow](https://www.nextflow.io/docs/latest/index.html), [Conda](https://docs.conda.io/en/latest/index.html) environments, and [Docker](https://docs.docker.com/) containers. To accommodate user environments, the pipeline containers are also compatible with [Singularity](https://sylabs.io/guides/3.6/user-guide/index.html).

## Dependencies

* This software has been tested on the following Linux operating systems: 64-bit CentOS 7 and 8, and Ubuntu 18.04.6, 20.04 LTS, 21.04, and 21.10.
    * The software may run on additional Linux distrubitions and/or versions beyond these if they are able to run the dependencies and versions listed below, but these are not officially supported.
* Internet connection
* [Nextflow](https://github.com/nextflow-io/nextflow#quick-start) (>=21.04.0)
* [Docker](https://docs.docker.com/engine/install/) (>=20.10.7) or [Singularity](https://sylabs.io/guides/3.6/user-guide/quick_start.html#quick-installation-steps) (>=3.6.4)
> **NOTE:** If using Docker, your `USER` must be added to the [`docker` root user group](https://docs.docker.com/engine/install/linux-postinstall/) before executing the pipeline. On shared systems, such as HPC clusters, this may not be possible due to security risks and the pipeline should be executed using the Singularity profile (default) instead.

## Quick Start

### Configure Credentials

> **NOTE:** Configuring credentials is only necessary while working with the private repository.

1. Create a GitHub account and make sure you have access to this repository.

1. Create a GitHub [personal access token](https://docs.github.com/en/github/authenticating-to-github/keeping-your-account-and-data-secure/creating-a-personal-access-token). This will provide command-line access to any private GitHub repositories that you have access to. Only check the `repo` box while selecting the scope, this will also check a handful of boxes below it. Create the token and copy the token value.
    > **NOTE:** Make sure to copy the token value as you will not be able to view it again after leaving the page. If you didn't copy the value, you will have to delete the token and make a new one.

1. Open your command-line terminal.

1. Connect to the computer where you'll be using the pipeline.

1. Add your GitHub username and token as temporary environment variables. Replace the text after `=` with the appropriate values.

    ```
    GITHUB_USER=<GitHub username>
    GITHUB_TOKEN=<GitHub personal access token>
    ```

1. Create a Nextflow directory. Nextflow will look here for your credentials and it will also serve as the default storage location for Nextflow pipelines.

    ```
    mkdir -p ~/.nextflow
    ```

1. Create a source code management (SCM) configuration file. This will allow Nextflow to access pipelines in your private repositories without requiring you to log in each time. 

    ```
    cat << EOF > ~/.nextflow/scm
    providers {
    
        github {
            user = "$GITHUB_USER"
            password = "$GITHUB_TOKEN"
        }
    
    }
    EOF
    ```

### Configure Parameters

The pipeline requires a YAML- or JSON-formatted file listing assay-specific parameters when not running the test data. A detailed description of the file contents, including the structure, are provided below. Assay parameters must be nested under the appropriate assay name (e.g. `atac`) using indentation and only one value may be listed per line. File paths must be relative to the directory where the analysis is being executed. If specifying multiple files, etc., you must create a new nested level of parameters. See below for an example of how to nest parameters as part of mixed species analysis.

#### Parameters

> **NOTE:** When configuring bead and cell calling settings, priority will be as follows:
>   1) Sample-specific settings within an assay-specific section
>   2) Assay-specific settings
>   3) Default settings

- **outputDir (default: `"./results"`):** Location to deposit final results and reports.
- **workflow:** Type of workflow to be executed. Must be "reference", "analysis", or "full".
    - **"reference":** Only perform the reference-generation portion of the given assay
        - Required parameters:
            - reference
                - directory
                - fasta
                - gtf
        - Optional parameters:
            - blocklist
            - mixed
            - barcodedTn5
            - tssWindowSize
    - **"analysis":** Only perform the analysis portion of the given assay
        - Required parameters:
            - input
            - reference
            - blocklist
        - Optional parameters:
            - mixed
            - barcodedTn5
            - ti
            - barcodedTn5Config
            - i7asti
            - tiread
            - tssWindowSize
            - mergeMethod
            - qualityThreshold
            - barcode
                - force
            - trim
            - sortSize
            - rounding
            - maxInsertSize
            - tierroroverride
            - overrides
        - **NOTE**: This workflow assumes that references have already been generated and will error if they are not found.
    - **"full":** Perform both reference-generation and analysis portions of the given assay
        - Required parameters:
            - reference
                - directory
                - fasta
                - gtf
            - input
        - Optional parameters:
            - blocklist
            - mixed
            - barcodedTn5
            - ti
            - barcodedTn5Config
            - i7asti
            - tiread
            - tssWindowSize
            - mergeMethod
            - qualityThreshold
            - barcode
                - force
            - trim
            - sortSize
            - maxInsertSize
            - tierroroverride
            - overrides


- **reference:** 
    - **directory:** Directory path to write newly-generated reference files to and/or containing pre-generated references (must be quoted)
    - **fasta:** File path(s) to reference FASTA file(s) (may be gzipped or decompressed; must be quoted)
    - **gtf:** File path(s) to reference GTF file(s) (may be gzipped or decompressed; must be quoted)
    - **blocklist:** File path(s) to reference blocklist BED file(s). ( the parameter must be quoted; see below for details)
- **input:** Directory path containing raw FASTQ files (may be gzipped or decompressed; must be quoted)
- **mixed:** Boolean value indicating if the mixed species workflow should be used (must be true or false)
    - **false:** Only allows a single FASTA and GTF file to be provided for reference workflow.
    - **true:** Requires two FASTA and GTF files to be provided for reference workflow. 
- **barcodedTn5:** Boolean value indicating if the samples have tagmentation indexes (must be true or false)
    - **false:** Samples are assumed to not have tagmentation indexes.
    - **true:** All samples are assumed to have tagmentation indexes and read 2 is parsed for the index.  Index reports are generated.
- **ti:** A list of user defined TI sequences.  Presets are present in conf/atac_preset.config.  See atac combinatorial config file below for example of formatting.
- **barcodedTn5Config:** File path to a CSV file that specifies the TIs and FASTQs assigned to each sample. If barcodedTn5 is set to true, and no config file is provided, the default behavior is to merge all TIs and FASTQs as a single sample. Only include TIs present in a given experiment as unused TIs that are included in the config can lead to errors. Must be formatted with these three columns:
    - **Sample** The name of the sample.
    - **Fastq** The name given to the fastq file pair.
    - **TI** The name of the TI used (see conf/atac_preset.config for preset list).
- **i7asti:** Boolean value indicating whether the I7 index should be treated as a tagmentation index. Requires **barcodedTn5** to be `true`.
- **tiread:** The read that contains the TI sequence, assuming the TI is not in the i7 sequence. Requires **barcodedTn5** to be `true` and **i7asti** to be `false`.(DEFAULT: "r1"; also accepts "r2")
- **tssWindowSize:** The full window size in bases around TSS for TSS enrichment score calculations (DEFAULT: 4000; must be an even integer)
- **mergeMethod:** The read(s) in a pair to use for determining the transposase insertion site used in bead merging. The default is "both", which uses the transposition site from both ends of the fragment. (DEFAULT: "both"; also accepts "r1" or "r2")
- **qualityThreshold:** Minimum mapq for a read to be kept.
- **barcode:** Assay-specific settings related to knee calling configuration.
    - **force:** The number of barcodes to return. This will bypass all knee calling algorithms and return the top `n` number of barcodes based on unique reads (must be an integer).
- **trim:** Number of bases to trim from the 5' end of R2 reads (DEFAULT: 0; must be an integer).
- **sortSize:** Set the sort collection size ratio for Mark Duplicates.  This parameter is to be used when the module runs out of memory due to high duplications.  Lower the number to make memory footprint size more manageable (DEFAULT: 0.25, must be float).  Recommended trial value is 0.01 in case of memory errors.
- **rounding:** Rounds insert sites to the nearest 10, 100, or 1000 bases before performing deconvolution. (DEFAULT: 0; no rounding is performed)
- **maxInsertSize:** The largest insert that the pipeline will recognize (DEFAULT:2000; must be at least 100)
- **tierroroverride:** If set to true, the pipeline will ignore errors in the TI config file.  Only fastq-TI combinations specified in the config file will be used.  All others will be ignored. (DEFAULT: false)
- **overrides: (non-combinatorial)** Sample-specific settings for overriding barcode calling and deconvolution configuration.
    - **<sample_id>:** The `sample_id` of the sample to be overriden.
        - **barcode:** See above for description.
            - **force**
        - **trim:** See above for description.
        - **mergeMethod:** See above for description.
- **overrides: (combinatorial)** FASTQ+TI-specific settings for overriding barcode calling and deconvoulution configuration. Overrides can be set on the FASTQ level or on the FASTQ+TI level. Overrides set on the FASTQ level will be applied to all FASTQ+TI combinations coming from that FASTQ file.
    - **<sample_id>:** The `FASTQ name` of the fastq to be overriden.
        - **barcode:** See above for description.
            - **force**
        - **trim:** See above for description.
        - **mergeMethod:** See above for description. Example: If the FASTQ files were named SampleA_S1_R1_001.fastq.gz and SampleA_S1_R1_001.fastq.gz, the FASTQ name would be SampleA_S1.
    - **<sample_id>:** The `FASTQ name` of the fastq to be overriden.
        - **<TI_name>:** The `TI name` of the TI to be overriden. Example: TI12
            - **barcode:** See above for description.
                - **force**
            - **trim:** See above for description.
            - **mergeMethod:** See above for description.

#### Examples

parameters_normal.yaml (atac normal)
```
atac:
  workflow: "full"
  input: "test/data/atac/normal/"
  reference: 
    directory: "test/references/atac/"
    fasta:
        species1: "test/references/atac/final_hg38.fa.gz"
        species2: "test/references/atac/final_mm10.fa.gz"
    gtf:
        species1: "test/references/atac/final_hg38.gtf.gz"
        species2: "test/references/atac/final_mm10.gtf.gz"
  mixed: true 
```

parameters_optional.yaml (atac combinatorial)
```
atac:
  workflow: "full"
  input: "test/data/atac/combinatorial/"
  reference: 
    directory:"test/references/atac/"
    fasta:
        species1: "test/references/atac/final_hg38.fa.gz"
        species2: "test/references/atac/final_mm10.fa.gz"
    gtf:
        species1: "test/references/atac/final_hg38.gtf.gz"
        species2: "test/references/atac/final_mm10.gtf.gz"
  mixed: true
  barcodedTn5: true
  ti:
    mycustomti1: "AGGTTATA"
    mycustomti2: "CTTGGTAT"
  barcodedTn5Config: "test/config/atac/atac_mixed_TIs.csv"
  tssWindowSize: 2000
  barcode:
    force: 1000
  trim: 5
  overrides:
    Cycler_S3_newindex:
      mergeMethod: r2
      mycustomti1:
        mergeMethod: r1
```

atac_mixed_TIs.csv
```
Sample,Fastq,TI
SampleA,Cycler_S3_newindex,mycustomti1
SampleB,Cycler_S3_newindex,mycustomt12
```

### Custom TI definition for combinatorial ATAC

To use custom TIs, add the TI sequences to your yaml config file (see atac combinatorial yaml for example).  Custom TIs must be of equal length and must be located in the 5' end of R2, immediately after the bead barcode in R1, or in the i7 sequencing barcode.  For TIs located in the 5' end of R2, if a constant sequence is located 3' of the TI sequence, please use the trim parameter to remove this sequence.  This will prevent the pipeline from attempting to align this sequence to the reference genome.  Set the trim parameter to the length of your constant sequence, not including the length of the TIs.  

### References

The pipline is only compatible with Ensembl formatted references.  It is not compatible with NCBI references at the present time.

### Execution

> **NOTE:** For pre-release versions or until the Docker containers are made public, you MUST login to Docker via `docker login` using credentials that have access to the bioraddbg Docker Hub organization before executing the pipeline and regardless of whether you are using Docker or Singularity. This only needs to be done once per `USER` per machine.

The pipeline may be executed using either Singularity (default) or Docker. Upon execution, Nextflow will download the pipeline from GitHub (the code may be found in `~/.nextflow/assets/BioRadOpenSource/omnition-test/`) and will download the Singularity or Docker containers as needed.

Please note the use of `-` and `--` in the execution commands. Arguments with a single `-` in front are Nextflow arguments and those with `--` are user-defined parameters. You may also use the Nextflow `-r` flag to specify a git tag (i.e. release), branch, or hash to execute the pipeline at that point in the git history.

> **NOTE:** When launching a nextflow run there will be a randomly generated [adjective_name] from the lists in the following code (https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/util/NameGenerator.groovy) that appears in the terminal. This is an inherent feature of nextflow not something that anyone at Bio-Rad added.

After completion, reports can be found in the `report/` subdirectory and intermediate files can be found in the `Sample_Files/` subdirectory. Additionally, a Nextflow cache directory (`.nextflow/`) and working directory (`work/`) will be created in the directory where Nextflow was executed from. This allows interrupted/failed analyses to resume from their point of failure. If you would like to continue a run, add `-resume` to the execution commands below. You can delete the cache and working directories to save space after completion though this will inhibit the `-resume` feature.


#### Singularity

    nextflow run BioRadOpenSource/omnition-test -params-file <path to parameters.yaml> --outputDir <output directory path>

#### Docker

    nextflow run BioRadOpenSource/omnition-test -params-file <path to parameters.yaml> --outputDir <output directory path> -profile docker
    
#### Test Data 

There are two different test profiles available to run with the pipeline. The `demo_atac` profile runs with the settings specified in the `parameters_normal.yaml` example. Similarly, the `demo_catac` profile runs with the settings specified in the `parameters_optional.yaml` example. The parameters used in these different testing profiles can be found in `conf/test.config`.

Below are examples of how to run the different testing profiles, while executing the pipeline with Docker or Singularity (specified as `standard` under the -profile flag. The test is specified following the `profile` flag.

Docker:

    nextflow run BioRadOpenSource/omnition-test --outputDir <output directory path> -profile demo_atac,docker
    nextflow run BioRadOpenSource/omnition-test --outputDir <output directory path> -profile demo_catac,docker

Singularity:

    nextflow run BioRadOpenSource/omnition-test --outputDir <output directory path> -profile demo_atac,standard
    nextflow run BioRadOpenSource/omnition-test --outputDir <output directory path> -profile demo_catac,standard

## Appendix

### ATAC BAM File Tags

This pipeline places a number of tags on the bam files produced in the workflow. Below is a table defining them.
| Tag | Type | First Annotated | Description |
|-----|------|-----------------|-------------|
| XB  | Z | `process tag_barcodes` | The cell barcode as determined by DEAD. Barcode and UMI are temporarily stored in the read name after `process dead` through `process bwa_alignment` and ultimately annotated during `process tag_barcodes`.
| PG  | Z | `process mark_duplicates` | Indicates MarkDuplicates was run. Mark duplicates updates the flag to indicate if a read is a duplicate.
| DB  | Z | `process decovolute (soon to be reannotate_bam)` | Deconvolute adds the cell name for merged barcodes as well as barcodes that do not need to be merged
| RG  | Z | `process tag_cells_above_knee` | If the cell name (see above) is above the knee cutoff, that tag is copied to the RG tag.
| AS | i | `process bwa_alignment` | Assigned by BWA and conforms to [SAM specifications](https://samtools.github.io/hts-specs/SAMtags.pdf): Alignment score from the aligner. Lay definition: See previous sentence.
| nM | i | `process bwa_alignment` | Assigned by BWA. The edit distance in the alignment.
| MC | Z | `process bwa_alignment` | Compressed representation of the alignment in the [CIGAR format].(https://www.drive5.com/usearch/manual/cigar.html)
| MD | Z | `process bwa_alignment` | Mismatching positions and bases.
| XS | i | `process bwa_alignment` | Alignment scores for suboptimal alignments.
| XA | Z | `process bwa_alignment` | Suboptimal alignment hits. Format: (chr,pos,CIGAR,NM;).

### Blocklist

The ATAC analysis workflow will ignore alignments to blocklisted regions within the genome during bead multiplet deconvolution. Bio-Rad does not provide blocklists. File name must be in the format SPECIESNAME.blocklist.bed where SPECIESNAME is the same species name used for the fasta and gtf files. Blocklists must be formatted as a three column BED file with feature names matching those in the reference genome FASTA and reference GTF:
```bash
chr10      0       45700
chr10      38481300        38596500
chr10      38782600        38967900
chr10      39901300        41712900
chr10      41838900        42107300
chr10      42279400        42322500
chr10      126946300       126953400
chr10      133625800       133797400
```
A blocklist is not required for the pipeline to execute. If one is not provided, an empty one is generated at runtime.

When doing a mixed species run, up to two blocklists may be supplied. 

Publicly available blocklists are available [here](https://github.com/Boyle-Lab/Blacklist). The feature names need to be updated to match those of the reference being used.
