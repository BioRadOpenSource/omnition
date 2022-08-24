# Bio-Rad Laboratories, Inc. Omnition Single-Cell Analysis Software

This analysis pipeline is designed to analyze single-cell ATAC-Seq data and combinatorial single-cell ATAC-Seq data. The pipeline is built using a combination of [Nextflow](https://www.nextflow.io/docs/latest/index.html), [Conda](https://docs.conda.io/en/latest/index.html) environments, and [Docker](https://docs.docker.com/) containers. To accommodate user environments, the pipeline containers are also compatible with [Singularity](https://sylabs.io/guides/3.6/user-guide/index.html).

## Dependencies

* This software has been tested on the following Linux operating systems: 64-bit CentOS 7 and 8, and Ubuntu 18.04.6, 20.04 LTS, 21.04, and 21.10.
    * The software may run on additional Linux distrubitions and/or versions beyond these if they are able to run the dependencies and versions listed below, but these are not officially supported.
* Internet connection
* [Nextflow](https://github.com/nextflow-io/nextflow#quick-start) (>=21.04.0)
* [Docker](https://docs.docker.com/engine/install/) (>=20.10.7) or [Singularity](https://sylabs.io/guides/3.6/user-guide/quick_start.html#quick-installation-steps) (>=3.6.4)

> **NOTE:** If using Docker, your `USER` must be added to the [`docker` root user group](https://docs.docker.com/engine/install/linux-postinstall/) before executing the pipeline. On shared systems, such as HPC clusters, this may not be possible due to security risks and the pipeline should be executed using the Singularity profile (default) instead.

> **NOTE:** The user must verify with their system administrator that Docker or Singularity is available before using Omnition.

## Quick Start

### Configure Parameters

The pipeline requires a YAML- or JSON-formatted file listing assay-specific parameters when not running the test data. A detailed description of the file contents, including the structure, are provided below. Assay parameters must be nested under the appropriate assay name (`atac`) using indentation and only one value may be listed per line. File paths must be relative to the directory where the analysis is being executed. If specifying multiple files, you must create a new nested level of parameters. See below for an example of how to nest parameters as part of a mixed species analysis.

#### Parameters

> **NOTE:** When configuring bead and cell calling settings, priority will be as follows:
>   1) Sample-specific settings within an assay-specific section
>   2) Assay-specific settings
>   3) Global settings
>   4) Default settings

#### Global

- **outputDir (default: `"./results"`):** Location to deposit final results and reports.

#### ATAC

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
            - mitoContig
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
            - blocklist
        - Optional parameters:
            - mitoContig
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
    - **directory:** Directory path to write newly-generated reference files to and/or containing pre-generated references. 
        - Required.
        - Must be quoted.
    - **fasta:** File path(s) to reference FASTA file(s). 
        - Required. 
        - Optionally gzip compressed; must be quoted.
    - **gtf:** File path(s) to reference GTF file(s). 
        - Required. 
        - Optionally gzip compressed; must be quoted.
    - **blocklist:** File path(s) to reference blocklist BED file(s). 
        - Optional. Must be quoted.
- **input:** Directory path containing raw FASTQ files. 
    - Required. 
    - Must be gzip compressed; must be quoted.
- **mitoContig:** The name of the mitochondrial contig in the reference. 
    - Optional.
    - DEFAULT: "MT".
- **mixed:** Boolean value indicating if the mixed species workflow should be used. 
    - Optional.
    - If `false`: Only allows a single FASTA and GTF file to be provided for reference workflow.
    - If `true:` Requires two FASTA and GTF files to be provided for reference workflow.
    - DEFAULT: `false`.
- **barcodedTn5:** Boolean value indicating if the samples have tagmentation indexes. 
    - Optional.
    - If `false`: Samples are assumed to not have tagmentation indexes.
    - If `true`: All samples are assumed to have tagmentation indexes and read 2 is parsed for the index.  Index reports are generated.
    - DEFAULT: `false`.
- **ti:** A list of user defined TI sequences. 
    - Optional. 
    - See ATAC Combinatorial Config file below for example of formatting.
    - DEFAULT: Presets are present in `conf/atac_preset.config`.
- **barcodedTn5Config:**  File path to a CSV file that specifies the TIs and FASTQs assigned to each sample.
    - Optional.
    - DEFAULT: None.
    - If barcodedTn5 is set to `true`, and no config file is provided, the default behavior is to merge all TIs and FASTQs as a single sample. Only include TIs present in a given experiment as unused TIs that are included in the config can lead to errors. Must be formatted with these three columns:
    - **Sample** The name of the sample.
    - **Fastq** The name given to the fastq file pair.
    - **TI**: The name of the TI used (see `conf/atac_preset.config` for list of preset TIs).
- **i7asti:** Boolean value indicating whether the I7 index should be treated as a tagmentation index. 
    - Optional. 
    - Requires **barcodedTn5** to be `true`.
    - DEFAULT: `false`.
- **tiread:** The read that contains the TI sequence, assuming the TI is not in the i7 sequence. 
    - Optional. 
    - Requires **barcodedTn5** to be `true` and **i7asti** to be `false`.
    - One of "r1" or "r2".
    - DEFAULT: "r1".
- **tssWindowSize:** The full window size in bases around TSS for TSS enrichment score calculations. 
    - Optional.
    - Must be an even integer greater than 0. 
    - DEFAULT: 4000.
- **mergeMethod:** The read(s) in a pair to use for determining the transposase insertion site used in bead merging. 
    - Optional.
    - One of "r1", "r2", or "both".
    - DEFAULT: "both".
- **qualityThreshold:** Minimum MAPQ for a read to be included in the analysis. 
    - Optional.
    - DEFAULT: 30.
- **barcode:** Settings related to knee calling configuration. 
    - Optional.
    - **force:** The number of barcodes to return. This will bypass all knee calling algorithms and return the top `n` number of barcodes based on unique reads.
        - An integer greater than zero.
        - DEFAULT: None. Omnition will determine the number of cells in each sample.
- **trim:** Number of bases to trim from the 5' end of R2 read. 
    - Optional.
    - An integer greater than or equal to 0.
    - DEFAULT: 0.
- **sortSize:** Set the sort collection size ratio for MarkDuplicates. 
    - Optional.
    - A floating point number greater than 0.  
    - This parameter is to be used when the module runs out of memory due to very high duplication rate.  
    - Lower the number to make memory footprint size more manageable. Recommended trial value is 0.01 in case of memory errors.
    - DEFAULT: 0.25.  
- **rounding:** Rounds insert sites to the nearest 10, 100, or 1000 bases before performing bead merging. 
    - Optional. 
    - One of 10, 100, or 1000.
    - DEFAULT: 0; no rounding is performed.
- **maxInsertSize:** The largest insert that Omnition will recognize.
    - Optional.
    - An integer greater than 100. 
    - DEFAULT: 2000.
- **tierroroverride:** Boolean value. 
    - If set to `true`, the pipeline will ignore errors in the TI configuration file. Only fastq-TI combinations specified in the config file will be used.  All others observed in the sequencing data will be ignored.
    - DEFAULT: `false`.
- **overrides: (non-combinatorial)** Sample-specific settings for overriding barcode calling and bead merging configuration.
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
        - **<TI_name>:** The `TI name` of the TI to be overriden. Example: TI12.
            - **barcode:** See above for description.
                - **force**
            - **trim:** See above for description.
            - **mergeMethod:** See above for description.

#### Examples

`parameters_atac.yaml` (ATAC-Seq)
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

`parameters_catac.yaml` (Combinatorial ATAC-Seq Superloading)
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
    ti1: "AAAGAA"
    ti2: "TTTGGG"
  barcodedTn5Config: "test/config/atac/atac_mixed_TIs.csv"
  tiread: "r1"
```
> NOTE: Removing the `barcodedTn5Config` parameter from above will change the configuration to Combinatorial ATAC-Seq Superloading.

atac_mixed_TIs.csv
```
Sample,Fastq,TI
SampleA,Cycler_S3_newindex,mycustomti1
SampleB,Cycler_S3_newindex,mycustomt12
```

### Custom Transposition Index (TI) definition for combinatorial ATAC-Seq

To use custom TIs, add the TI sequences to your yaml config file. See `parameters_catac.yaml`, above, for an example.  

Custom TIs must be of equal length and must be located in one of the following locations:
- The 5' end of R2 (see parameter `tiread: r2`)
- Immediately after the bead barcode in R1 (set parameter `tiread: r1`)
- In the i7 sequencing barcode (see parameter `i7asti: true`).  

For TIs located in the 5' end of R2, use the `trim` parameter to remove all bases upstream of the TI sequence. 

### References

Omnition is compatible with [Human/GRCh38](https://uswest.ensembl.org/Homo_sapiens/Info/Index) and [Mouse/GRCm39](http://uswest.ensembl.org/Mus_musculus/Info/Index) references from ENSEMBL.  It is not compatible with other species or references produced by NCBI.

### Execution

The pipeline may be executed using either Singularity (default) or Docker. Upon execution, Nextflow will download the pipeline from GitHub. The code may be found in `~/.nextflow/assets/BioRadOpenSource/omnition/` after download. Omnition will download the Singularity images or Docker containers as needed.

Please note the use of single `-` and double `--` in the execution commands. Arguments with a single `-` in front are Nextflow arguments and those with double  `--` are user-defined parameters.

> **NOTE:** Nextflow gives each pipeline run a randomly generated name formatted as `[adjective_name]`. Names are taken from those listed in the following code (https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/util/NameGenerator.groovy). These run names are generated by Nextflow and not endorsed by Bio-Rad.

After completion, reports can be found in the `report/` subdirectory and intermediate files can be found in the `Sample_Files/` subdirectory. Additionally, a Nextflow cache directory (`.nextflow/`) and working directory (`work/`) will be created in the directory where Nextflow was executed. 

#### Nextflow Functionality

##### Resuming an Analysis
Interrupted or failed analyses can be restarted from their stopping point by using the `-resume` flag in the pipeline execution commands. Deleting the cache (`.nextflow/`) or working directory (`work/`) will clear the cache and disable the `-resume` functionality.

#### Diagnostics
Pipeline diagnostic reports generated by Nextflow will be produced by default by Omnition. For each analysis, a `pipeline_info` directory will be generated in the `--outputDir` and populated with a Nextflow execution report, trace report, and timeline report. For more information on these reports, see the [Nextflow documentation](https://www.nextflow.io/docs/latest/tracing.html).

#### Singularity

    nextflow run BioRadOpenSource/omnition -params-file <path to parameters.yaml> --outputDir <output directory path>

#### Docker

    nextflow run BioRadOpenSource/omnition -params-file <path to parameters.yaml> --outputDir <output directory path> -profile docker
    
#### Test Data 

There are two test profiles available to run with the pipeline. The profile `demo_atac` runs with the settings specified in `parameters_atac.yaml`, above. The profile `demo_catac` runs with the settings specified in the `parameters_catac.yaml` example, above. The parameters used in these different testing profiles can be found in `conf/test.config`.

Below are examples of how to run the different testing profiles while executing the pipeline with Docker (`-profile docker`) or Singularity (`-profile standard`). The test type is then specified separated by a comma in the `-profile` flag.

Docker:

    nextflow run BioRadOpenSource/omnition --outputDir <output directory path> -profile docker,demo_atac
    nextflow run BioRadOpenSource/omnition --outputDir <output directory path> -profile docker,demo_catac

Singularity:

    nextflow run BioRadOpenSource/omnition --outputDir <output directory path> -profile standard,demo_atac
    nextflow run BioRadOpenSource/omnition --outputDir <output directory path> -profile standard,demo_catac

## Appendix

### Resource Allocation
Omnition is composed of a series of “processes” that connect to form a pipeline. Each process is allocated a portion of the host system’s CPU and RAM when Nextflow’s scheduler launches a processing task. Nextflow’s scheduler launches processes in parallel when the system’s capacity permits it. 

Bio-Rad has provided default resource allocations for the processes in Omnition to optimize performance on a system with 16 CPU and 64 GB RAM. If a process fails due to insufficient resources, it will be retried a maximum of two times with each retry requesting double the resources of the previous one. Should the process fail after a second retry, no further retries will be carried out. Each process is assigned two labels that dictate the CPU and memory allocation to that process. The default labels are defined in the tables below. 

#### CPU Allocation Labels
| Label | CPUs |
|-------|------|
| `cpu_xsmall` | 1 |
| `cpu_small` | 2 |
| `cpu_medium` | 4 |
| `cpu_large` | 8 |
| `cpu_xlarge` | 16 |

#### Memory Allocation Labels
| Label | RAM (GB) |
|-------|----------|
| `memory_xxsmall` | 0.9375 |
| `memory_xsmall` | 3.75 |
| `memory_small` | 7.5 |
| `memory_medium` | 15 |
| `memory_large` | 30 |
| `memory_xlarge` | 60 |

#### Default Process Labels
| Process | CPUs | RAM |
|---------|------|-----|
| AGGREGATE_METRICS | `cpu_small` | `memory_xsmall` |
| ANNOTATE_FRAGMENTS| `cpu_medium` | `memory_medium` |
| ARCHR | `cpu_small` | `memory_xlarge` |
| ARCHR_REFERENCE | `cpu_small` | `memory_large` |
| ASSEMBLE_BASIC_QC | `cpu_medium` | `memory_xxsmall` |
| ASSEMBLE_FRAGMENTS | `cpu_xsmall` | `memory_xsmall` |
| BEAD_FILT_SUMMARY | `cpu_small` | `memory_xxsmall` |
| BUILD_REPORT_CONTENTS | `cpu_medium` | `memory_xxsmall` |
| BWA_ALIGNMENT | `cpu_xlarge` | `memory_xlarge` |
| BWA_INDEX | `cpu_small` | `memory_large` |
| CALCULATE_BEADS_PER_DROP | `cpu_xsmall` | `memory_xxsmall` |
| CALCULATE_INSERT_SIZE_METRICS | `cpu_small` | `memory_xsmall` |
| CALL_PEAKS | `cpu_xsmall` | `memory_medium` |
| CHECK_MITO_CONTIG | `cpu_xsmall` | `memory_xsmall` |
| CHECK_TI_COUNTS | `cpu_xsmall` | `memory_xxsmall` |
| CHECK_TI_COUNTS_SUPERLOADED | `cpu_xsmall` | `memory_xxsmall` |
| CLEAN_PEAKS | `cpu_xsmall` | `memory_xsmall` |
| COMBINE_BLOCKLISTS | `cpu_xsmall` | `memory_xxsmall` |
| COMBINE_REFERENCES | `cpu_small` | `memory_xxsmall` |
| COMPILE_ALIGNMENTS | `cpu_large` | `memory_xsmall` |
| COMPILE_FRAGMENTS | `cpu_large` | `memory_medium` |
| COMPILE_QC_STATS | `cpu_large` | `memory_xxsmall` |
| COMPUTE_DECONVOLUTION_STAT_CHR | `cpu_medium` | `memory_medium` |
| COMPUTE_TSS_MATRIX | `cpu_small` | `memory_xxsmall` |
| CUTADAPT_HEADCROP | `cpu_medium` | `memory_medium` |
| DEAD | `cpu_small` | `memory_xxsmall` |
| DETERMINE_BARCODE_ALLOWLIST | `cpu_medium` | `memory_small` |
| DETERMINE_BARCODE_MERGES | `cpu_medium` | `memory_small` |
| FASTQC | `cpu_small` | `memory_small` |
| FILTER_BLOCKLISTS | `cpu_xsmall` | `memory_xxsmall` |
| FILTER_REFERENCES | `cpu_xsmall` | `memory_xxsmall` |
| FINAL_BAM_MERGE | `cpu_xlarge` | `memory_xxsmall` |
| FINAL_FRAG_MERGE | `cpu_medium` | `memory_xsmall` |
| FINAL_QC_SE | `cpu_medium` | `memory_xlarge` |
| FORMAT_BLOCKLIST | `cpu_xsmall` | `memory_xxsmall` |
| FRACTION_OF_READS_IN_PEAKS | `cpu_xsmall` | `memory_xsmall` |
| FRACTION_OF_READS_IN_TSS | `cpu_xsmall` | `memory_xsmall` |
| GENERATE_EMPTY_BLOCKLIST | `cpu_xsmall` | `memory_xxsmall` |
| GENERATE_GENOME_SIZES | `cpu_medium` | `memory_xxsmall` |
| GENERATE_REPORT | `cpu_small` | `memory_small` |
| GENERATE_TSS_WINDOWS | `cpu_medium` | `memory_xxsmall` |
| GUNZIP | `cpu_xsmall` | `memory_xxsmall` |
| MAKE_COUNT_MATRIX | `cpu_xsmall` | `memory_xlarge` |
| MARK_DUPLICATES | `cpu_medium` | `memory_xlarge` |
| MERGE_LANES | `cpu_xsmall` | `memory_xxsmall` |
| MERGE_REANN_READ_COUNTS | `cpu_small` | `memory_xxsmall` |
| PUBLISH_PARAMETERS | `cpu_medium` | `memory_xxsmall` |
| REANNOTATE_BAM | `cpu_medium` | `memory_xxsmall` |
| REANNOTATE_FRAGMENTS | `cpu_medium` | `memory_medium` |
| SEQUENCE_SATURATION | `cpu_medium` | `memory_medium` |
| SPLIT_BAM | `cpu_xlarge` | `memory_xsmall` |
| SPLIT_FASTQ | `cpu_xsmall` | `memory_xxsmall` |
| SUMMARIZE_ALIGNMENTS | `cpu_medium` | `memory_xsmall` |
| SUMMARIZE_MIXED_SPECIES | `cpu_xsmall` | `memory_xxsmall` |
| TAG_BARCODES | `cpu_xlarge` | `memory_xlarge` |
| TI_DEAD_CONFIG | `cpu_xsmall` | `memory_xxsmall` |
| TI_ERROR_CHECK | `cpu_xsmall` | `memory_xxsmall` |
| TI_WARNING_MESSAGES | `cpu_xsmall` | `memory_xsmall` |
| TSS_ENRICHMENT | `cpu_xsmall` | `memory_xxsmall` |
| VALIDATE_TI_CONFIG | `cpu_xsmall` | `memory_xxsmall` |

#### Adjusting Resource Allocation
To change the resource allocation for a given resource tier (e.g., increase memory for the memory_xlarge label), create a configuration file like the one below and save it as `resources.config`:

```
process{ 
    withLabel: memory_xlarge {  
        memory = 128.GB  
        time = 24.h  
    }
} 
```

This file should then be saved and passed to Omnition when running the software, for example: 

`nextflow run BioRadOpenSource/omnition -params-file analysis.yaml -profile standard –c resources.config` 

The resource allocations in `resources.config` will then override the default settings. 

### ATAC BAM File Tags

This pipeline places a number of tags on the bam files produced in the workflow. Below is a table defining them.
| Tag | Type | First Annotated | Description |
|-----|------|-----------------|-------------|
| AS | i | `process BWA_ALIGNMENT` | Assigned by BWA and conforms to [SAM specifications](https://samtools.github.io/hts-specs/SAMtags.pdf): Alignment score from the aligner. Lay definition: See previous sentence.
| NM | i | `process BWA_ALIGNMENT` | Assigned by BWA. The edit distance in the alignment.
| MC | Z | `process BWA_ALIGNMENT` | Compressed representation of the alignment's mate in the [CIGAR format].(https://www.drive5.com/usearch/manual/cigar.html)
| MD | Z | `process BWA_ALIGNMENT` | Mismatching positions and bases.
| XS | i | `process BWA_ALIGNMENT` | Alignment scores for suboptimal alignments.
| XA | Z | `process BWA_ALIGNMENT` | Suboptimal alignment hits. Format: (chr,pos,CIGAR,NM).
| XB  | Z | `process BWA_ALIGNMENT` | The bead barcode sequence for a read.
| PG  | Z | `process MARK_DUPLICATES` | Indicates MarkDuplicates was executed.
| DI | i | `process MARK_DUPLICATES` | The duplicate set index as assigned by MarkDuplicates.
| DB  | Z | `process REANNOTATE_BAM` | The cell barcode.

### Blocklist

> NOTE: Bio-Rad has adopted the term "blocklist” as a replacement for “blacklist.”

The ATAC analysis workflow will ignore alignments to blocklisted regions within the genome during bead merging. Bio-Rad does not provide blocklists. File name must be in the format SPECIESNAME.blocklist.bed where  SPECIESNAME is the same species name used for the fasta and gtf files. Blocklists must be formatted as a three column BED file with feature names matching those in the reference genome FASTA and reference GTF:
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

### Benchmarks

The dataset used for quantification was a mixed species run with 1,262,814,294 total paired reads which were downsampled to 4 different levels (100%, 75%, 50%, 25%, 10%) when run through the analysis workflow. Without downsampling the pipeline called 32,812 total cells.

A system with 16 CPUs (Intel Xeon Cascade Lake) and 128 GB RAM was (AWS `r5d.4xlarge`) used for the runs used to generate these results.

Reference workflow
| name | realtime (run time) | peak_rss (peak RAM usage) |
|------|---------------------|---------------------------|
| ARCHR_REFERENCE | 9m 2s | 11.3 GB |
| BWA_INDEX | 2h 25s | 8.2 GB |
| COMBINE_REFERENCES | 24s | 41.5 MB |
| FILTER_BLOCKLISTS | 18.5s | 66.7 MB |
| FILTER_REFERENCES | 26.3s | 63.8 MB |
| FORMAT_BLOCKLIST | 4s | 57.3 MB |
| GENERATE_GENOME_SIZES | 23.1s | 31.2 MB |
| GENERATE_TSS_WINDOWS | 9.8s | 137.6 MB |
| GUNZIP_FASTA | 27.2s | 58.4 MB |
| GUNZIP_GTF | 9.3s | 59.1 MB |

Analysis workflow

| % of sample reads used | total paired reads |
|------------------------|--------------------|
| 100 | 1,262,814,294 |
| 75 | 947,110,720 |
| 50 | 631,407,146 |
| 25 | 315,703,572 |
| 10 | 126,281,428 |

| name | realtime (run time) | peak_rss (peak RAM usage) | % of initial sample reads used |
|------|---------------------|---------------------------|--------------------------------|
| AGGREGATE_METRICS | 3m 12s | 767.2 MB | 100 |
| AGGREGATE_METRICS | 2m 47s | 759.1 MB | 75 |
| AGGREGATE_METRICS | 2m 6s | 706.3 MB | 50 |
| AGGREGATE_METRICS | 1m 39s | 643.3 MB | 25 |
| AGGREGATE_METRICS | 2m 29s | 932.1 MB | 10 |
| ANNOTATE_FRAGMENTS | 9s | 957.4 MB | 100 |
| ANNOTATE_FRAGMENTS | 9s | 996.6 MB | 75 |
| ANNOTATE_FRAGMENTS | 9s | 982.8 MB | 50 |
| ANNOTATE_FRAGMENTS | 9s | 650.1 MB | 25 |
| ANNOTATE_FRAGMENTS | 7s | 353.3 MB | 10 |
| ARCHR | 1h 51m 47s | 68.4 GB | 100 |
| ARCHR | 1h 50m 11s | 61.2 GB | 75 |
| ARCHR | 1h 42m 24s | 76.3 GB | 50 |
| ARCHR | 1h 31m 43s | 70.1 GB | 25 |
| ARCHR | 1h 25m 14s | 51.1 GB | 10 |
| ASSEMBLE_BASIC_QC | 4s | 352 MB | 100 |
| ASSEMBLE_BASIC_QC | 4s | 359.7 MB | 75 |
| ASSEMBLE_BASIC_QC | 4s | 351.9 MB | 50 |
| ASSEMBLE_BASIC_QC | 4s | 360.4 MB | 25 |
| ASSEMBLE_BASIC_QC | 4s | 361.9 MB | 10 |
| ASSEMBLE_FRAGMENTS | 9s | 998 MB | 100 |
| ASSEMBLE_FRAGMENTS | 9s | 989.6 MB | 75 |
| ASSEMBLE_FRAGMENTS | 9s | 99.1 MB | 50 |
| ASSEMBLE_FRAGMENTS | 9s | 99.8 MB | 25 |
| ASSEMBLE_FRAGMENTS | 8s | 98.9 MB | 10 |
| BEAD_FILT_SUMMARY | 11.4s | 14.8 MB | 100 |
| BEAD_FILT_SUMMARY | 9.1s | 14.6 MB | 75 |
| BEAD_FILT_SUMMARY | 6.3s | 14.5 MB | 50 |
| BEAD_FILT_SUMMARY | 3.2s | 14.6 MB | 25 |
| BEAD_FILT_SUMMARY | 2.4s | 14.7 MB | 10 |
| BUILD_REPORT_CONTENTS | 2.5s | 122.4 MB | 100 |
| BUILD_REPORT_CONTENTS | 2.6s | 121.9 MB | 75 |
| BUILD_REPORT_CONTENTS | 2.5s | 123.2 MB | 50 |
| BUILD_REPORT_CONTENTS | 2.5s | 122.4 MB | 25 |
| BUILD_REPORT_CONTENTS | 2.6s | 121.8 MB | 10 |
| BWA_ALIGNMENT | 38m 22s | 29.6 GB | 100 |
| BWA_ALIGNMENT | 43m 50s | 29.6 GB | 75 |
| BWA_ALIGNMENT | 9m 54s | 29.2 GB | 50 |
| BWA_ALIGNMENT | 9m 44s | 23 GB | 25 |
| BWA_ALIGNMENT | 33m 19s | 16.5 GB | 10 |
| CALCULATE_BEADS_PER_DROP | 2s | 37.3 MB | 100 |
| CALCULATE_BEADS_PER_DROP | 2s | 36.3 MB | 75 |
| CALCULATE_BEADS_PER_DROP | 2s | 36.4 MB | 50 |
| CALCULATE_BEADS_PER_DROP | 2s | 35.5 MB | 25 |
| CALCULATE_BEADS_PER_DROP | 3s | 74.8 MB | 10 |
| CALCULATE_INSERT_SIZE_METRICS | 25m 58s | 1.1 GB | 100 |
| CALCULATE_INSERT_SIZE_METRICS | 20m 8s | 1.1 GB | 75 |
| CALCULATE_INSERT_SIZE_METRICS | 13m 54s | 972.2 MB | 50 |
| CALCULATE_INSERT_SIZE_METRICS | 8m 29s | 667.4 MB | 25 |
| CALCULATE_INSERT_SIZE_METRICS | 4m 5s | 1.4 GB | 10 |
| CALL_PEAKS | 44m 6s | 3.3 GB | 100 |
| CALL_PEAKS | 37m 25s | 2.7 GB | 75 |
| CALL_PEAKS | 29m 26s | 2.2 GB | 50 |
| CALL_PEAKS | 17m 51s | 1.4 GB | 25 |
| CALL_PEAKS | 9m 9s | 736.3 MB | 10 |
| CHECK_TI_COUNTS_SUPERLOADED | 41m 20s | 14.8 MB | 100 |
| CHECK_TI_COUNTS_SUPERLOADED | 32m 12s | 28 MB | 75 |
| CHECK_TI_COUNTS_SUPERLOADED | 21m 27s | 33.4 MB | 50 |
| CHECK_TI_COUNTS_SUPERLOADED | 10m 35s | 32.6 MB | 25 |
| CHECK_TI_COUNTS_SUPERLOADED | 1m 57s | 14.9 MB | 10 |
| CLEAN_PEAKS | 13.1s | 372.7 MB | 100 |
| CLEAN_PEAKS | 12.4s | 415.2 MB | 75 |
| CLEAN_PEAKS | 12s | 409.5 MB | 50 |
| CLEAN_PEAKS | 10.8s | 372.9 MB | 25 |
| CLEAN_PEAKS | 9.5s | 355.7 MB | 10 |
| COMPILE_ALIGNMENTS | 2h 11m 50s | 403.7 MB | 100 |
| COMPILE_ALIGNMENTS | 1h 56m 13s | 435.4 MB | 75 |
| COMPILE_ALIGNMENTS | 1h 26m 7s | 419.7 MB | 50 |
| COMPILE_ALIGNMENTS | 50m 8s | 385.8 MB | 25 |
| COMPILE_ALIGNMENTS | 18m 34s | 630.5 MB | 10 |
| COMPILE_QC_STATS | 1.1s | 13.1 MB | 100 |
| COMPILE_QC_STATS | 1.1s | 11.8 MB | 75 |
| COMPILE_QC_STATS | 1s | 12 MB | 50 |
| COMPILE_QC_STATS | 1.4s | 12.4 MB | 25 |
| COMPILE_QC_STATS | 958ms | 12 MB | 10 |
| COMPUTE_DECONVOLUTION_STAT_CHR | 9s | 996.6 MB | 100 |
| COMPUTE_DECONVOLUTION_STAT_CHR | 9s | 99.5 MB | 75 |
| COMPUTE_DECONVOLUTION_STAT_CHR | 9s | 99.9 MB | 50 |
| COMPUTE_DECONVOLUTION_STAT_CHR | 9s | 99.6 MB | 25 |
| COMPUTE_DECONVOLUTION_STAT_CHR | 9.9s | 99.9 MB | 10 |
| COMPUTE_TSS_MATRIX | 1h 56s | 1.7 GB | 100 |
| COMPUTE_TSS_MATRIX | 51m 54s | 1.5 GB | 75 |
| COMPUTE_TSS_MATRIX | 44m 37s | 1.2 GB | 50 |
| COMPUTE_TSS_MATRIX | 34m 42s | 802.5 MB | 25 |
| COMPUTE_TSS_MATRIX | 26m 26s | 653.2 MB | 10 |
| COUNT_BEADS_PER_PARTITION | 3.8s | 98.5 MB | 100 |
| COUNT_BEADS_PER_PARTITION | 7.6s | 98.3 MB | 75 |
| COUNT_BEADS_PER_PARTITION | 6s | 98 MB | 50 |
| COUNT_BEADS_PER_PARTITION | 7.1s | 99.7 MB | 25 |
| COUNT_BEADS_PER_PARTITION | 5.6s | 98.6 MB | 10 |
| CUTADAPT_HEADCROP | 8m 31s | 4.3 GB | 100 |
| CUTADAPT_HEADCROP | 7m 49s | 4.1 GB | 75 |
| CUTADAPT_HEADCROP | 5m 10s | 4 GB | 50 |
| CUTADAPT_HEADCROP | 2m 5s | 3.3 GB | 25 |
| CUTADAPT_HEADCROP | 36.7s | 991.6 MB | 10 |
| DEAD | 3h 31m 53s | 135.7 MB | 100 |
| DEAD | 2h 35m 53s | 135.6 MB | 75 |
| DEAD | 1h 43m 36s | 135.7 MB | 50 |
| DEAD | 51m 51s | 135.6 MB | 25 |
| DEAD | 12m 42s | 135.7 MB | 10 |
| DETERMINE_BARCODE_ALLOWLIST | 3m 5s | 1.2 GB | 100 |
| DETERMINE_BARCODE_ALLOWLIST | 2m 13s | 1.2 GB | 75 |
| DETERMINE_BARCODE_ALLOWLIST | 60s | 1.1 GB | 50 |
| DETERMINE_BARCODE_ALLOWLIST | 59s | 1.1 GB | 25 |
| DETERMINE_BARCODE_ALLOWLIST | 41.1s | 996.7 MB | 10 |
| DETERMINE_BARCODE_MERGES | 34.3s | 658 MB | 100 |
| DETERMINE_BARCODE_MERGES | 9.9s | 570.5 MB | 75 |
| DETERMINE_BARCODE_MERGES | 9.7s | 369.7 MB | 50 |
| DETERMINE_BARCODE_MERGES | 9.1s | 198.4 MB | 25 |
| DETERMINE_BARCODE_MERGES | 5s | 163 MB | 10 |
| FASTQC | 2h 19m 22s | 596.2 MB | 100 |
| FASTQC | 1h 44m 49s | 609.9 MB | 75 |
| FASTQC | 1h 10m 15s | 602.5 MB | 50 |
| FASTQC | 36m 25s | 577.1 MB | 25 |
| FASTQC | 8m 50s | 568.3 MB | 10 |
| FINAL_BAM_MERGE | 57.1s | 44.2 MB | 100 |
| FINAL_BAM_MERGE | 59.1s | 47.8 MB | 75 |
| FINAL_BAM_MERGE | 58s | 46.8 MB | 50 |
| FINAL_BAM_MERGE | 58.9s | 45.1 MB | 25 |
| FINAL_BAM_MERGE | 9.9s | 44.1 MB | 10 |
| FINAL_FRAG_MERGE | 59s | 43.3 MB | 100 |
| FINAL_FRAG_MERGE | 58.7s | 45.1 MB | 75 |
| FINAL_FRAG_MERGE | 49.9s | 8.7 MB | 50 |
| FINAL_FRAG_MERGE | 24.8s | 43 MB | 25 |
| FINAL_FRAG_MERGE | 8.1s | 55.7 MB | 10 |
| FINAL_QC_SE | 59.5s | 4.7 GB | 100 |
| FINAL_QC_SE | 48.9s | 3 GB | 75 |
| FINAL_QC_SE | 39.6s | 2.6 GB | 50 |
| FINAL_QC_SE | 36.2s | 988.2 MB | 25 |
| FINAL_QC_SE | 9.8s | 762.6 MB | 10 |
| FRACTION_OF_READS_IN_PEAKS | 2h 2m 32s | 522.3 MB | 100 |
| FRACTION_OF_READS_IN_PEAKS | 1h 39m 30s | 553 MB | 75 |
| FRACTION_OF_READS_IN_PEAKS | 1h 15m 56s | 528.1 MB | 50 |
| FRACTION_OF_READS_IN_PEAKS | 45m 3s | 481.4 MB | 25 |
| FRACTION_OF_READS_IN_PEAKS | 20m 45s | 732.7 MB | 10 |
| FRACTION_OF_READS_IN_TSS | 2h 20m 34s | 512.1 MB | 100 |
| FRACTION_OF_READS_IN_TSS | 1h 54m 39s | 548.5 MB | 75 |
| FRACTION_OF_READS_IN_TSS | 1h 26m 45s | 532.8 MB | 50 |
| FRACTION_OF_READS_IN_TSS | 52m 11s | 500.9 MB | 25 |
| FRACTION_OF_READS_IN_TSS | 24m 9s | 775.7 MB | 10 |
| GENERATE_REPORTS | 11m 55s | 8.1 GB | 100 |
| GENERATE_REPORTS | 9m 32s | 7.1 GB | 75 |
| GENERATE_REPORTS | 6m 22s | 5.9 GB | 50 |
| GENERATE_REPORTS | 2m 48s | 1.6 GB | 25 |
| GENERATE_REPORTS | 2m 11s | 4.2 GB | 10 |
| MAKE_COUNT_MATRIX | 14m 41s | 20 GB | 100 |
| MAKE_COUNT_MATRIX | 12m 56s | 18.6 GB | 75 |
| MAKE_COUNT_MATRIX | 9m 37s | 16.4 GB | 50 |
| MAKE_COUNT_MATRIX | 6m 20s | 12.9 GB | 25 |
| MAKE_COUNT_MATRIX | 3m 31s | 13.5 GB | 10 |
| MARK_DUPLICATES | 28m 31s | 48.8 GB | 100 |
| MARK_DUPLICATES | 9m 1s | 42.2 GB | 75 |
| MARK_DUPLICATES | 9m 49s | 37.6 GB | 50 |
| MARK_DUPLICATES | 7m 4s | 35.8 GB | 25 |
| MARK_DUPLICATES | 1m 9s | 31.7 GB | 10 |
| MERGE_LANES | 7m 30s | 65.3 MB | 100 |
| MERGE_LANES | 8m 22s | 53.9 MB | 75 |
| MERGE_LANES | 2m 28s | 40.1 MB | 50 |
| MERGE_LANES | 2m 56s | 41.5 MB | 25 |
| MERGE_LANES | 3.7s | 55.6 MB | 10 |
| MERGE_REANN_READ_COUNTS | 3.9s | 126.8 MB | 100 |
| MERGE_REANN_READ_COUNTS | 5.9s | 127 MB | 75 |
| MERGE_REANN_READ_COUNTS | 5.5s | 126.9 MB | 50 |
| MERGE_REANN_READ_COUNTS | 6.3s | 126.8 MB | 25 |
| MERGE_REANN_READ_COUNTS | 5.4s | 126.9 MB | 10 |
| REANNOTATE_BAM | 9s | 98 MB | 100 |
| REANNOTATE_BAM | 9s | 97.8 MB | 75 |
| REANNOTATE_BAM | 9s | 97.2 MB | 50 |
| REANNOTATE_BAM | 9s | 90.3 MB | 25 |
| REANNOTATE_BAM | 9s | 6.1 MB | 10 |
| REANNOTATE_FRAGMENTS | 9s | 998.2 MB | 100 |
| REANNOTATE_FRAGMENTS | 9s | 985.5 MB | 75 |
| REANNOTATE_FRAGMENTS | 9s | 895.8 MB | 50 |
| REANNOTATE_FRAGMENTS | 9.5s | 583.7 MB | 25 |
| REANNOTATE_FRAGMENTS | 5s | 340.7 MB | 10 |
| SEQUENCE_SATURATION | 2m 5s | 7.5 GB | 100 |
| SEQUENCE_SATURATION | 1m 52s | 6.7 GB | 75 |
| SEQUENCE_SATURATION | 1m 32s | 5.6 GB | 50 |
| SEQUENCE_SATURATION | 1m 4s | 4 GB | 25 |
| SEQUENCE_SATURATION | 37.6s | 2.2 GB | 10 |
| SPLIT_BAM | 3m 51s | 943.4 MB | 100 |
| SPLIT_BAM | 3m 34s | 944 MB | 75 |
| SPLIT_BAM | 41.5s | 933.8 MB | 50 |
| SPLIT_BAM | 59.9s | 924.6 MB | 25 |
| SPLIT_BAM | 9.7s | 798.1 MB | 10 |
| SPLIT_FASTQ | 2h 9m 36s | 67.1 MB | 100 |
| SPLIT_FASTQ | 1h 42m 17s | 64.3 MB | 75 |
| SPLIT_FASTQ | 59m 57s | 65.7 MB | 50 |
| SPLIT_FASTQ | 28m 28s | 65.6 MB | 25 |
| SPLIT_FASTQ | 6m 12s | 63.1 MB | 10 |
| SUMMARIZE_ALIGNMENTS | 2h 54m 18s | 1.3 GB | 100 |
| SUMMARIZE_ALIGNMENTS | 2h 28m | 1.2 GB | 75 |
| SUMMARIZE_ALIGNMENTS | 1h 43m 13s | 940.6 MB | 50 |
| SUMMARIZE_ALIGNMENTS | 59m 45s | 1.4 GB | 25 |
| SUMMARIZE_ALIGNMENTS | 32m 56s | 1.4 GB | 10 |
| SUMMARIZE_MIXED_SPECIES | 3.1s | 96.3 MB | 100 |
| SUMMARIZE_MIXED_SPECIES | 2.9s | 117.3 MB | 75 |
| SUMMARIZE_MIXED_SPECIES | 3.1s | 105.7 MB | 50 |
| SUMMARIZE_MIXED_SPECIES | 3.1s | 100.1 MB | 25 |
| SUMMARIZE_MIXED_SPECIES | 3.8s | 121.9 MB | 10 |
| TAG_BARCODES | 5m 19s | 14.2 GB | 100 |
| TAG_BARCODES | 7m 48s | 14.1 GB | 75 |
| TAG_BARCODES | 4m 34s | 9.7 GB | 50 |
| TAG_BARCODES | 59.7s | 9.5 GB | 25 |
| TAG_BARCODES | 44.1s | 2.1 GB | 10 |
| TI_DEAD_CONFIG | 3.4s | 12.6 MB | 100 |
| TI_DEAD_CONFIG | 6.7s | 59.8 MB | 75 |
| TI_DEAD_CONFIG | 4.9s | 47.9 MB | 50 |
| TI_DEAD_CONFIG | 5s | 54.1 MB | 25 |
| TI_DEAD_CONFIG | 6.8s | 55 MB | 10 |
| TSS_ENRICHMENT | 4m 23s | 16.4 MB | 100 |
| TSS_ENRICHMENT | 4m 21s | 16.4 MB | 75 |
| TSS_ENRICHMENT | 4m 19s | 16.5 MB | 50 |
| TSS_ENRICHMENT | 4m 17s | 16.3 MB | 25 |
| TSS_ENRICHMENT | 4m 15s | 16.3 MB | 10 |
