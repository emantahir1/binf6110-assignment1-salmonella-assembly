#!/bin/bash
# Salmonella enterica Genome Assembly and Variant Analysis
# BINF6110 Assignment 1 Part 2
# Author: Eman Tahir
# Date: February 12, 2026

set -e  # Exit on error

echo "======================================"
echo "Salmonella enterica Genome Analysis"
echo "======================================"

# Project directory
PROJECT_DIR=~/binf6110_assign1
cd $PROJECT_DIR

# ============================================
# STEP 1: DATA ACQUISITION
# ============================================
echo "Step 1: Data already downloaded"
echo "  - Raw reads: salmonella_reads.fastq (196,031 reads)"
echo "  - Reference: salmonella_ref.fasta (NC_003197.2 + NC_003277.2)"

# ============================================
# STEP 2: GENOME ASSEMBLY
# ============================================
echo "Step 2: De novo assembly with Flye"
flye --nano-hq salmonella_reads.fastq \
     --out-dir flye_output \
     --threads 8

echo "Assembly complete!"
cat flye_output/assembly_info.txt

# ============================================
# STEP 3: REFERENCE ALIGNMENT
# ============================================
echo "Step 3: Aligning reads to reference genome"

# Align with minimap2
minimap2 -ax lr:hq -t 8 \
  salmonella_ref.fasta \
  salmonella_reads.fastq > reads_aligned.sam

# Convert to sorted BAM
samtools view -bS reads_aligned.sam | samtools sort -o reads_sorted.bam
samtools index reads_sorted.bam

# Index reference
samtools faidx salmonella_ref.fasta

# Get coverage statistics
samtools coverage reads_sorted.bam > detailed_coverage.txt

echo "Alignment complete!"

# Clean up large SAM file
rm reads_aligned.sam

# ============================================
# STEP 4: VARIANT CALLING
# ============================================
echo "Step 4: Variant calling with Bcftools"

bcftools mpileup -Ou -f salmonella_ref.fasta reads_sorted.bam | \
  bcftools call -mv -Ob -o reads_variants.bcf

bcftools view reads_variants.bcf > reads_variants.vcf

# Get variant statistics
echo "Total variants:"
bcftools view -H reads_variants.vcf | wc -l

echo "SNPs:"
bcftools view -H reads_variants.vcf | \
  awk '{if(length($4)==1 && length($5)==1) print}' | wc -l

echo "Insertions:"
bcftools view -H reads_variants.vcf | \
  awk '{if(length($4)<length($5)) print}' | wc -l

echo "Deletions:"
bcftools view -H reads_variants.vcf | \
  awk '{if(length($4)>length($5)) print}' | wc -l

# Extract variant data for plotting
bcftools view -H reads_variants.vcf | \
  awk '{print $1, $2, length($4), length($5), $6}' > variant_data.txt

# Get genome depth for visualization
samtools depth reads_sorted.bam > genome_depth.txt

echo "Variant calling complete!"

# ============================================
# STEP 5: VISUALIZATION
# ============================================
echo "Step 5: Creating visualizations with R"

# Run R scripts to generate figures
Rscript figure1_assembly.R
Rscript figure2_alignment.R
Rscript figure3_variants.R

echo "Figures generated!"

# ============================================
# ANALYSIS COMPLETE
# ============================================
echo "======================================"
echo "Analysis Complete!"
echo "======================================"
echo ""
echo "Results:"
echo "  - Assembly: flye_output/assembly.fasta (3 contigs, 5.10 Mb)"
echo "  - Variants: reads_variants.vcf (11,465 variants)"
echo "  - Figures: figure1_assembly_final.png"
echo "            figure2_alignment_final.png"
echo "            figure3_variants_final.png"
echo ""
echo "For IGV visualization:"
echo "  1. Load salmonella_ref.fasta"
echo "  2. Load reads_sorted.bam"
echo "  3. Load reads_variants.vcf"
