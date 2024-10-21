# hg38_based_wes_pipeline
Somatic variant calling pipeline based on GATK best practice guide for hg38# WES Workflow Automation Script

This repository contains a shell script that automates a whole-exome sequencing (WES) workflow. The script processes FASTQ files, performs alignment, recalibration, variant calling, and filtering, and outputs the final set of high-quality variants. The workflow is designed to minimize manual input and can be easily configured using a text-based input file.

## Features
- **Automated processing** of raw sequencing data (FASTQ) to filtered variants (VCF).
- Supports alignment, duplicate marking, base quality recalibration, and variant calling.
- Uses widely adopted bioinformatics tools such as **GATK**, **SAMTOOLS**, **PICARD**, and **BOWTIE2**.
- Simple configuration via an `input_paths.txt` file for easy reuse and modification.
- Outputs are generated in a structured directory based on the sample name for better organization.

## Prerequisites
Before using this script, ensure you have the following software installed:
- **Conda** with the `Mamba` package manager.
- Bioinformatics tools required for the workflow, such as GATK, SAMTOOLS, PICARD, and BOWTIE2, installed in a Conda environment.

### Required Tools:
- GATK
- SAMTOOLS
- PICARD
- BOWTIE2
- BCFTOOLS

## Setup

### Step 1: Install the Mamba Package Manager
The first step is to ensure `mamba` is installed in your base conda environment to manage packages efficiently.

```bash
conda install -n base --override-channels -c conda-forge mamba 'python_abi=*=*cp*'
```

### Step 2: Create and Set Up the Environment
Create a Conda environment named `wes_env` and install all the necessary tools.

```bash
mamba create -n wes_env gatk4 bwa samtools picard bowtie2
```

Activate the environment:
```bash
conda activate wes_env
```

### Step 3: Clone the Repository
Clone this GitHub repository to your local machine.

```bash
git clone https://github.com/yourusername/wes_workflow.git
cd wes_workflow
```

### Step 4: Prepare the Input Configuration File
The workflow requires an input configuration file (`input_paths.txt`) where you specify the paths for input data, reference genome, known sites, and output directory.

Here’s an example of what the `input_paths.txt` file should look like:

```txt
SAMPLE_NAME=HD701_s1_S5
INPUT_DIR=/media/afroz/NEW_VOLUME1/GATK_practice/runs/HD701_s1/hg38
OUTPUT_DIR=/media/afroz/NEW_VOLUME1/GATK_practice/outputs
REF_GENOME_DIR=/media/afroz/NEW_VOLUME1/GATK_practice/ref/genome_38
INDEX_PREFIX=/media/afroz/NEW_VOLUME1/GATK_practice/ref/genome/hg38_index
KNOWN_SITES_DIR=/media/afroz/NEW_VOLUME1/GATK_practice/ref/known_sites_38
```

You will need to modify this file according to your specific input files, reference genome, and known sites directories.

### Step 5: Run the Script
To run the workflow, first ensure the `wes_env` environment is activated. The script will prompt you to activate it as well.

```bash
./wes_workflow.sh
```

The script will:
1. Read the input paths from `input_paths.txt`.
2. Process the FASTQ files into BAM, perform alignment, base quality score recalibration (BQSR), and call variants.
3. Output the results (e.g., recalibrated BAM, VCF files) in a directory structure organized by the sample name.

### Output
- The output files will be generated in the directory specified by `OUTPUT_DIR` in the configuration file.
- For each sample, a folder will be created with the sample name, containing:
  - BAM files
  - Recalibration reports
  - VCF files (SNPs, Indels, merged variants)
  - Final filtered variants in VCF format

## Input Configuration File
The `input_paths.txt` configuration file is crucial for running the script smoothly. Below is an explanation of each field:

| Field | Description |
|-------|-------------|
| `SAMPLE_NAME` | Name of the sample (used to locate the FASTQ files and name the output files) |
| `INPUT_DIR` | Directory where the input FASTQ files are stored |
| `OUTPUT_DIR` | Parent directory where the output files should be stored |
| `REF_GENOME_DIR` | Directory containing the reference genome (FASTA file) |
| `INDEX_PREFIX` | Prefix of the Bowtie2 index files for the reference genome |
| `KNOWN_SITES_DIR` | Directory containing the known sites files for BQSR (such as hapmap, 1000G, Mills, etc.) |

## Workflow Details
The script performs the following steps:

1. **Convert FASTQ to BAM** using `gatk FastqToSam`.
2. **Mark Illumina Adaptors** in the unaligned BAM file using `gatk MarkIlluminaAdapters`.
3. **Convert BAM back to FASTQ** using `gatk SamToFastq` for alignment.
4. **Align reads** using `bowtie2` and sort them using `samtools`.
5. **Mark duplicates** in the aligned BAM using `picard MarkDuplicates`.
6. **Add read group information** using `picard AddOrReplaceReadGroups`.
7. **Prepare reference dictionary, fasta index, and BAM index** using `picard` and `samtools`.
8. **Perform base quality score recalibration** using `gatk BaseRecalibrator` and `gatk ApplyBQSR`.
9. **Call variants** using `gatk Mutect2`.
10. **Filter SNPs and Indels** using `gatk VariantFiltration`.
11. **Merge and extract variants** using `gatk MergeVcfs` and `bcftools`.

## Example Usage
Once the environment is set up and the `input_paths.txt` is properly configured, simply run the script:
```bash
./wes_workflow.sh
```

This will automate the entire WES analysis workflow for the sample specified in `input_paths.txt`.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributions
Contributions, issues, and feature requests are welcome! Feel free to fork this repository and submit pull requests.

## Contact
For any queries, reach out to [Afroz Shaik](afrozshaik2157@gmail.com).



