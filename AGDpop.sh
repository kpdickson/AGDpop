
#!/bin/bash

set -euo pipefail

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Code available at: https://github.com/kpdickson/AGDpop
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ==============================================================================

# AGD Pipeline: FastQ --> Filtered VCFs
# Description: Maps paired-end reads to reference, generates BAMs, performs
#              joint variant calling followed by filtering.

# ==============================================================================

# --- 1. SETUP PATHS & PARAMETERS ---
REF="./AGDref.fasta"
FASTQ_DIR="./paired_fastq"
BAM_DIR="./processed_bams"
ID_DIR="./ref" 
OUT_DIR="./final_VCFs"
THREADS=64

# --- VCF FILTER PARAMETERS ---
QUAL=30
MISS=0.3
MIN_DP=3

mkdir -p $BAM_DIR $OUT_DIR

# --- 2. INDEX REFERENCE ---
if [ ! -f "${REF}.bwt" ]; then
    echo "Reference index not found. Indexing now..."
    bwa index $REF
    samtools faidx $REF
fi

# --- 3. MAPPING ---
echo "Starting Mapping..."
for R1 in ${FASTQ_DIR}/*_R1*.paired.fq; do
    R2="${R1/_R1/_R2}"
    SAMPLE=$(basename "$R1" | sed 's/_R1.*//')
    
    if [ ! -f "$R2" ]; then
        echo "Warning: R2 for $SAMPLE not found ($R2). Skipping."
        continue
    fi

    # Skip mapping if BAM already exists
    if [ -f "${BAM_DIR}/${SAMPLE}_sorted.bam" ]; then
        echo "BAM for $SAMPLE already exists, skipping."
        continue
    fi

    echo "Processing: $SAMPLE"
    bwa mem -t $THREADS -M -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" $REF "$R1" "$R2" | \
    samtools view -h -q 30 -Sb - | \
    samtools sort -@ $THREADS -o ${BAM_DIR}/${SAMPLE}_sorted.bam -
    
    samtools index ${BAM_DIR}/${SAMPLE}_sorted.bam
done

# --- 4. JOINT VARIANT CALLING ---
echo -e "\nGenerating BAM list..."
ls ${BAM_DIR}/*_sorted.bam > ${OUT_DIR}/bam.list
echo "Running joint calling on $(wc -l < ${OUT_DIR}/bam.list) viable samples..."

# Loop through 4 split categories
for TYPE in para_nuc perk_nuc para_mito perk_mito; do
    echo -e "\n--- Processing: $TYPE ---"
    
    ID_FILE="${ID_DIR}/${TYPE}.ids"
    if [ ! -f "$ID_FILE" ]; then 
        echo "Error: Region file $ID_FILE not found! Skipping."
        continue 
    fi

    RAW_VCF="${OUT_DIR}/cohort_${TYPE}_raw.vcf.gz"
    FINAL_VCF="${OUT_DIR}/final_${TYPE}_DP${MIN_DP}_miss${MISS}.vcf.gz"

    # Joint Calling
    bcftools mpileup -Ou -f $REF -R $ID_FILE -a FORMAT/AD,FORMAT/DP -b ${OUT_DIR}/bam.list | \
    bcftools call -f GQ,GP -mv -Oz -o $RAW_VCF

    # Filtering
    bcftools filter -S . -e "FMT/DP < ${MIN_DP}" $RAW_VCF -Ou | \
    bcftools filter -i "QUAL >= ${QUAL} && F_MISSING < ${MISS}" -Oz -o $FINAL_VCF
    
    bcftools index $FINAL_VCF
    
    echo "Results for $TYPE:"
    echo "  Final Filtered SNPs: $(bcftools view -H $FINAL_VCF | wc -l)"
done

echo -e "\nPipeline Complete."




