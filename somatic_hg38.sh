#!/bin/bash

# Step 1: Prompt user to activate "wes_env" environment
read -p "Please activate the 'wes_env' environment before running the script. Press Enter to continue..."

# Step 2: Read input from a file
CONFIG_FILE="input_paths.txt"

# Check if the config file exists
if [ ! -f "$CONFIG_FILE" ]; then
  echo "Configuration file '$CONFIG_FILE' not found!"
  exit 1
fi

# Load the configuration file
source $CONFIG_FILE

# Define paths and file names based on inputs from the configuration file
FASTQ_R1="${INPUT_DIR}/${SAMPLE_NAME}_L001_R1_001.fastq.gz"
FASTQ_R2="${INPUT_DIR}/${SAMPLE_NAME}_L001_R2_001.fastq.gz"
REF_GENOME="${REF_GENOME_DIR}/Homo_sapiens_assembly38.fasta"
OUTPUT_SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE_NAME}"

# Ensure output directory exists
mkdir -p $OUTPUT_SAMPLE_DIR

# Define other file paths based on sample name
BAM_OUTPUT="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}.bam"
ADAPTOR_BAM="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}_illumina_marked.bam"
FASTQ_OUTPUT="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}.fastq.gz"
ALIGNED_BAM="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}_sorted_aligned_reads.bam"
DUPLICATE_BAM="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}_flagged_reads.bam"
FINAL_BAM="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}_with_RG.bam"
RECAL_REPORT="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}_BaseRecalReport.grp"
BQSR_BAM="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}_BQSR.bam"
VCF_OUTPUT="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}.vcf.gz"
LEFT_ALIGN_VCF="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}_LeftAlign.vcf"
SNP_VCF="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}_SNP.vcf"
INDEL_VCF="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}_INDEL.vcf"
SELECTED_SNP_VCF="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}_Select_SNP.vcf"
SELECTED_INDEL_VCF="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}_Select_INDEL.vcf"
MERGED_VCF="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}_Merged.vcf"
FINAL_VCF="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}_All_Pass.vcf"

# Step 3: Convert FASTQ to BAM
echo "Converting FASTQ to BAM..."
gatk FastqToSam \
  FASTQ=$FASTQ_R1 \
  FASTQ2=$FASTQ_R2 \
  OUTPUT=$BAM_OUTPUT \
  READ_GROUP_NAME="${SAMPLE_NAME}_reads" \
  SAMPLE_NAME=$SAMPLE_NAME \
  LIBRARY_NAME=HG38_based \
  PLATFORM_UNIT=NextSeq \
  PLATFORM=illumina \
  SEQUENCING_CENTER=4GEN.bio

# Step 4: Mark Illumina Adaptors
echo "Marking Illumina adaptors..."
gatk MarkIlluminaAdapters \
  I=$BAM_OUTPUT \
  O=$ADAPTOR_BAM \
  M="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}_illumina_metrices.txt"

# Step 5: Convert BAM to FASTQ
echo "Converting BAM to FASTQ..."
gatk SamToFastq \
  I=$ADAPTOR_BAM \
  FASTQ=$FASTQ_OUTPUT \
  CLIPPING_ATTRIBUTE=XT \
  CLIPPING_ACTION=2 \
  INTERLEAVE=true

# Step 6: Alignment using Bowtie2
echo "Performing alignment..."
bowtie2 -x $INDEX_PREFIX -U $FASTQ_OUTPUT | samtools sort -@ 8 -o $ALIGNED_BAM

# Step 7: Mark Duplicates
echo "Marking duplicates..."
picard MarkDuplicates \
  INPUT=$ALIGNED_BAM \
  OUTPUT=$DUPLICATE_BAM \
  METRICS_FILE="${OUTPUT_SAMPLE_DIR}/${SAMPLE_NAME}_flagged_reads_metrics.txt"

# Step 8: Add or Replace Read Groups
echo "Adding header information..."
picard AddOrReplaceReadGroups \
  I=$DUPLICATE_BAM \
  O=$FINAL_BAM \
  RGID=1 \
  RGLB=library \
  RGPL=illumina \
  RGPU=machine \
  RGSM=$SAMPLE_NAME

# Step 9: Prepare reference dictionary and indexes
echo "Preparing reference dictionary and indexes..."
picard CreateSequenceDictionary \
  R=$REF_GENOME \
  O="${REF_GENOME_DIR}/Homo_sapiens_assembly38.dict"

samtools faidx $REF_GENOME
samtools index $FINAL_BAM

# Step 10: Base Quality Score Recalibration (BQSR)
echo "Performing BQSR..."
gatk BaseRecalibrator \
  -R $REF_GENOME \
  -I $FINAL_BAM \
  --known-sites ${KNOWN_SITES_DIR}/hapmap_3.3.hg38.vcf.gz \
  --known-sites ${KNOWN_SITES_DIR}/1000G_omni2.5.hg38.vcf.gz \
  --known-sites ${KNOWN_SITES_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz \
  --known-sites ${KNOWN_SITES_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -O $RECAL_REPORT

gatk ApplyBQSR \
  -R $REF_GENOME \
  -I $FINAL_BAM \
  --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
  -O $BQSR_BAM \
  --bqsr $RECAL_REPORT

# Step 11: Variant Calling using Mutect2
echo "Performing variant calling..."
gatk Mutect2 \
  -R $REF_GENOME \
  -I $BQSR_BAM \
  -O $VCF_OUTPUT

# Step 12: Left Align Variants
echo "Left aligning variants..."
gatk LeftAlignAndTrimVariants \
  -R $REF_GENOME \
  -V $VCF_OUTPUT \
  -O $LEFT_ALIGN_VCF

# Step 13: Separate SNPs and INDELs
echo "Separating SNPs..."
gatk SelectVariants \
  -R $REF_GENOME \
  -V $LEFT_ALIGN_VCF \
  --select-type-to-include SNP \
  -O $SNP_VCF

echo "Separating INDELs..."
gatk SelectVariants \
  -R $REF_GENOME \
  -V $LEFT_ALIGN_VCF \
  --select-type-to-include INDEL \
  -O $INDEL_VCF

# Step 14: Filtering SNPs and INDELs
echo "Filtering INDELs..."
gatk VariantFiltration \
  -R $REF_GENOME \
  -V $INDEL_VCF \
  --filter-name "my_filter_INDEL" \
  --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || QUAL < 30" \
  -O $SELECTED_INDEL_VCF

echo "Filtering SNPs..."
gatk VariantFiltration \
  -R $REF_GENOME \
  -V $SNP_VCF \
  --filter-name "my_filter_SNP" \
  --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || QUAL < 30" \
  -O $SELECTED_SNP_VCF

# Step 15: Merging SNPs and INDELs
echo "Merging SNP and INDEL files..."
gatk MergeVcfs \
  -I $SELECTED_SNP_VCF \
  -I $SELECTED_INDEL_VCF \
  -O $MERGED_VCF

# Step 16: Extracting Passed Variants
echo "Extracting variants that passed filters..."
bcftools view -i 'FILTER="PASS"' $MERGED_VCF > $FINAL_VCF

echo "Workflow completed successfully!"

