# Bio-Rad Laboratories, Inc. Omnition Single-Cell Analysis Software



## **Introduction to Omnition**
Omnition is a pipeline designed to process data generated with Bio-Rad's ddSEQ™ Single-Cell 3’ RNA-Seq Kit or ddSEQ™ SureCell ATAC-Seq Library Prep Kit and the dsciATAC protocol. It performs debarcoding, alignment, bead merging, cell calling, and feature counting. Output is an HTML report and files for downstream biological analysis. 

This page is intended to be a quick start. For the complete Omnition user guides, please see: 

[Omnition Single-cell 3’ RNA-Seq](https://www.bio-rad.com/webroot/web/html/general/omnition/Omnition_HTML5_v1.1.0_RNA/Content/Intro/1.1_Online_Introduction_SW.htm) 

[Omnition Single-cell ATAC-Seq](https://www.bio-rad.com/webroot/web/html/general/omnition/Omnition_HTML5_v1.1.0_ATAC/Content/Intro/1.1_Online_Introduction_SW.htm)
<br><hr>
<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#introduction-to-omnition">Introduction to Omnition</a></li>
    <li><a href="#installation">Installation</a>
      <ul>
        <li><a href="#system-requirements">System requirements</a></li>
        <li><a href="#hardware-requirements">Hardware requirements</a></li>
        <li><a href="#software-requirements">Software requirements</a>
        <li><a href="#download-omnition">Download Omnition</a></li>
        <li><a href="#verify-installation">Verify installation</a></li>
      </ul>
    </li>
    <li><a href="#getting-started">Getting started</a>
      <ul>
        <li><a href="#references">References</a>
          <ul>
            <li><a href="#human-genome">Human genome</a></li>
            <li><a href="#mouse-genome">Mouse genome</a></li>
          </ul>
        </li>
        <li><a href="#running-omnition">Running Omnition</a>
          <ul>
            <li><a href="#yaml-configuration">YAML configuration</a></li>
            <li><a href="#run-atac-pipeline">Run ATAC pipeline</a></li>
            <li><a href="#run-rna-pipeline">Run RNA pipeline</a></li>
          </ul>
        </li>
        <li><a href="#outputs">Outputs</a>
      </ul>
    </li>
  </ol>
</details>

# **Installation**
Omnition Analysis Software utilizes the Nextflow framework to connect individual processes, and runs the
processes in virtual environments called containers using the Docker or Singularity container programs. Once Omnition is installed on a system meeting minimum requirements, Nextflow integration with GitHub and Docker will prepare workflows without any additional configuration.

## **System requirements**
* Omnition is designed to run on a local Linux server, high performance computing (HPC) cluster, or cloud virtual machine, and has been tested on the 64-bit CentOS 7 and 8, Amazon Linux 2, and Ubuntu 18.04.6, 20.04 LTS, 21.04, and 21.10 Linux operating systems.
  > **NOTE**: Although they might be functional, Bio-Rad does not support additional Linux variants or other versions of the specified operating systems.

* Internet connection
  > **NOTE**: For installation without direct internet access, please see user manual.
## **Hardware requirements**

| Requirement      | ATAC seq analysis | Combinatorial ATAC seq analysis | RNA seq analysis | Recommended for >=12x samples |
|:-----------------|:-----------------:|:-------------------------------:|:----------------:|:-----------------------------:|
| CPU              | 16                | 16                              | 16               | 64                            |
| RAM              | 64 GB             | 128 GB                          | 64 GB            | 512 GB                        |
| IOPS*            | 3,000             | 3,000                           | 3,000            | 16,000                        |
| I/O throughput*  | 125 mbps          | 125 mbps                        | 125 mbps         | 1,000 mbps                    |
| EBS volume type* | gp3               | gp3                             | gp3              | gp3                           |

*AWS specific cloud computing specifications 

## **Software requirements**

* [Nextflow](https://github.com/nextflow-io/nextflow#quick-start) (v22.04.0 to v23.10.1) or [Nextflow with Conda](https://anaconda.org/bioconda/nextflow)
  > **NOTE:** Specify the Nextflow version in the Conda installation by using the following command: `conda install –c bioconda nextflow=<version>`

* Only one of the following container programs are needed: Either [Docker](https://docs.docker.com/get-docker/) (>=20.10.7) or [Singularity](https://sylabs.io/guides/3.6/user-guide/quick_start.html#quick-installation-steps) (>=3.6.4)
  > **NOTE:** If using Docker, your `USER` must be added to the [docker root user group](https://docs.docker.com/engine/install/linux-postinstall/) before executing the pipeline. On shared systems, such as HPC clusters, this may not be possible due to security risks and the pipeline should be executed using the Singularity profile (default) instead. The user must verify with their system administrator that Docker or Singularity is available before using Omnition.

## **Download Omnition**
To download and run Omnition, use nextflow to retrieve the latest Omnition version from GitHub.
```
nextflow pull BioRadOpenSource/omnition
```

## **Verify installation**
Omnition includes small demonstration datasets to verify that the environment has been properly built and all software dependencies are in place. To verify the success of the installation for each analysis type, run the Nextflow command for the container system that is installed on your computer (Singularity or Docker).

To verify each analysis workflow is installed correctly, run the following commands with either Docker or Singularity.
> **NOTE:** Working files (and output files, unless otherwise specified) will be generated in the same directory the pipeline was run from on the command line.

```
mkdir /home/ubuntu/demo_data
cd /home/ubuntu/demo_data
```
> **NOTE:** Specified file paths are for example purposes

**Docker:**

    # Verify single and mixed species 3’ RNA workflows
    nextflow run BioRadOpenSource/omnition -profile demo_rna_single,docker --rna.output /home/ubuntu/demo_data/rna
    nextflow run BioRadOpenSource/omnition -profile demo_rna_mixed_options,docker --rna.output /home/ubuntu/demo_data/rna_mixed

    # Verify ATAC-seq and combinatorial ATAC-seq workflows
    nextflow run BioRadOpenSource/omnition -profile demo_atac,docker --atac.output /home/ubuntu/demo_data/atac
    nextflow run BioRadOpenSource/omnition -profile demo_catac,docker --catac.output /home/ubuntu/demo_data/catac

**Singularity:**

    # Verify single and mixed species 3’ RNA workflows
    nextflow run BioRadOpenSource/omnition -profile demo_rna_single,standard --rna.output /home/ubuntu/demo_data/rna
    nextflow run BioRadOpenSource/omnition -profile demo_rna_mixed_options,standard --rna.output /home/ubuntu/demo_data/rna_mixed
    
    # Verify ATAC-seq and combinatorial ATAC-seq workflows
    nextflow run BioRadOpenSource/omnition -profile demo_atac,standard --atac.output /home/ubuntu/demo_data/atac
    nextflow run BioRadOpenSource/omnition -profile demo_catac,standard --catac.output /home/ubuntu/demo_data/catac

> **NOTE:** Please note the use of `-` and `--` in the execution commands. Arguments with a single `-` in front are Nextflow arguments and those with `--` are user-defined parameters. You may also use the Nextflow -r flag to specify a git tag (i.e. release), branch, or hash to execute the pipeline at that point in the git history.

> **NOTE:** When launching a nextflow run there will be a randomly generated [adjective_names] that appear in the terminal. These names are from a [list prepared by Nextflow](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/util/NameGenerator.groovy). This is an inherent feature of nextflow, not something added by Bio-Rad.

When the run finishes, you should see a message that it has completed without any failed tasks.


# **Getting started**

Before running a workflow, Omnition requires the following files to run:
- Genome FASTA and GTF files from ENSEMBL.
- Genome reference sequences must be formatted as FASTA files.
- Annotations must be formatted as GTF files.
- Sequence names in the FASTA and GTF files must match.
- Directory with input FASTQs

> **NOTE:** Input reference files can be compressed as gzip (.gz) files

## **References**
Omnition is compatible with references from ENSEMBL. The ENSEMBL [Human/GRCh38](https://ensembl.org/Homo_sapiens/Info/Index) and [Mouse/GRCm39](http://ensembl.org/Mus_musculus/Info/Index) references are supported by Omnition. 
Other species from ENSEMBL are not supported and references produced by other sources (NCBI, GENCODE, etc.) are not compatible.

Human and mouse references can be obtained with the following commands.

### **Human genome**
```
mkdir ~/references/human
cd ~/references/human
```
> **NOTE:** Specified file paths are for example purposes

```
curl -o Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

curl -o Homo_sapiens.GRCh38.106.gtf.gz \
http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz

```
### **Mouse Genome**
```
mkdir ~/references/mouse
cd ~/references/mouse
```
```
curl -o Mus_musculus.GRCm39.dna.primary_assembly.fa.gz \
http://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

curl -o Mus_musculus.GRCm39.106.gtf.gz \
http://ftp.ensembl.org/pub/release-106/gtf/mus_musculus/Mus_musculus.GRCm39.106.gtf.gz

```

## **Running Omnition**

### **YAML configuration**

The pipeline requires a YAML-formatted file with input/output paths and assay-specific parameters when not running the test data. A detailed description of the file contents, available parameters, and how to format them can be found in the user manual. Example YAML configs can be found under the [example-yamls](./example-yamls) folder of this repo.

### **Run RNA pipeline**
With Singularity:
```
nextflow run BioRadOpenSource/omnition -params-file <path to rna_example.yaml>
```
With Docker:
```
nextflow run BioRadOpenSource/omnition -params-file <path to rna_example.yaml> -profile docker
```

### **Run ATAC pipeline**

With Singularity:
```
nextflow run BioRadOpenSource/omnition -params-file <path to atac_example.yaml>
```
With Docker:
```
nextflow run BioRadOpenSource/omnition -params-file <path to atac_example.yaml> -profile docker
```

## **Outputs**
After completion, reports can be found in the `results/report/` subdirectory and intermediate files can be found in the `results/Sample_Files/` subdirectory. Additionally, a Nextflow cache directory (`.nextflow/`) and working directory (`work/`) will be created in the directory where Nextflow was executed from. This allows interrupted/failed analyses to resume from their point of failure. If you would like to continue a run, add `-resume` to the execution commands below. You can delete the cache and working directories to save space after completion though this will inhibit the `-resume` feature.



<br><hr>
[Return to top](#introduction-to-omnition)

