#!/bin/bash
# ==============================================================================
# Salmonella enterica Genome Assembly and Variant Analysis
# Complete Workflow - BINF6110 Assignment 1 Part 2
# ==============================================================================
# Author: Eman Tahir
# Date: February 12, 2026
# Description: Complete pipeline from environment setup to final visualization
# ==============================================================================

set -e  # Exit on any error

echo "======================================"
echo "SALMONELLA GENOME ANALYSIS PIPELINE"
echo "======================================"

# ==============================================================================
# PREREQUISITES
# ==============================================================================
# Required software (install via conda):
# - flye
# - minimap2
# - samtools
# - bcftools
# - seqtk
# - R with tidyverse, ggplot2, patchwork

# ==============================================================================
# STEP 0: ENVIRONMENT SETUP
# ==============================================================================
echo "Step 0: Setting up conda environment"

# Create conda environment (if not already created)
# conda create -n binf6110_env python=3.12
# conda activate binf6110_env

# Install bioinformatics tools
# conda install -c bioconda flye minimap2 samtools bcftools seqtk -y
# conda install -c conda-forge r-base r-essentials -y

# In R, install packages:
# install.packages("tidyverse", repos="http://cran.rstudio.com/")
# install.packages("ggplot2", repos="http://cran.rstudio.com/")
# install.packages("patchwork", repos="http://cran.rstudio.com/")

echo "Conda environment: binf6110_env"
echo "Ensure environment is activated: conda activate binf6110_env"

# ==============================================================================
# STEP 1: DATA ACQUISITION
# ==============================================================================
echo ""
echo "Step 1: Data Acquisition"
echo "========================"

# Create project directory
PROJECT_DIR=~/binf6110_assign1
mkdir -p $PROJECT_DIR
cd $PROJECT_DIR

# Download raw sequencing reads from NCBI SRA
echo "Downloading raw sequencing data from NCBI SRA..."
# wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR32410565/SRR32410565

# Convert SRA to FASTQ (if needed)
# fastq-dump --outdir . --gzip SRR32410565
# gunzip SRR32410565.fastq.gz
# mv SRR32410565.fastq salmonella_reads.fastq

# For this analysis, data files are already present:
echo "Using existing data files:"
echo "  - salmonella_reads.fastq (196,031 reads, 809 Mb)"

# Download reference genome
echo "Downloading reference genome..."
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz
# gunzip GCF_000006945.2_ASM694v2_genomic.fna.gz
# mv GCF_000006945.2_ASM694v2_genomic.fna salmonella_ref.fasta

echo "Reference: Salmonella enterica serovar Typhimurium LT2 (GCF_000006945.2)"
echo "  - Chromosome NC_003197.2 (4,857,450 bp)"
echo "  - Plasmid NC_003277.2 (93,933 bp)"

# ==============================================================================
# STEP 2: QUALITY CONTROL
# ==============================================================================
echo ""
echo "Step 2: Quality Control"
echo "======================="

# Get read statistics
echo "Calculating read statistics..."
seqkit stats salmonella_reads.fastq -a -T > read_statistics.txt

echo "Read Statistics:"
cat read_statistics.txt

# Extract read quality and length for visualization
seqkit fx2tab -q -l salmonella_reads.fastq > read_quality_length.txt

echo "Quality control complete!"

# ==============================================================================
# STEP 3: GENOME ASSEMBLY
# ==============================================================================
echo ""
echo "Step 3: De Novo Genome Assembly"
echo "================================"

echo "Running Flye assembler..."
flye --nano-hq salmonella_reads.fastq \
     --out-dir flye_output \
     --threads 8 \
     --genome-size 5m

echo "Assembly complete!"
echo ""
echo "Assembly Statistics:"
cat flye_output/assembly_info.txt

# ==============================================================================
# STEP 4: REFERENCE ALIGNMENT
# ==============================================================================
echo ""
echo "Step 4: Reference Genome Alignment"
echo "==================================="

# Index reference genome
echo "Indexing reference genome..."
samtools faidx salmonella_ref.fasta

# Align raw reads to reference with minimap2
echo "Aligning reads to reference genome with Minimap2..."
minimap2 -ax lr:hq -t 8 \
  salmonella_ref.fasta \
  salmonella_reads.fastq > reads_aligned.sam

echo "Alignment complete! Converting to BAM format..."

# Convert SAM to sorted BAM
samtools view -bS reads_aligned.sam | samtools sort -o reads_sorted.bam

# Index BAM file
samtools index reads_sorted.bam

# Calculate coverage statistics
samtools coverage reads_sorted.bam > detailed_coverage.txt

echo "Alignment Statistics:"
cat detailed_coverage.txt

# Clean up large SAM file to save space
rm reads_aligned.sam

# ==============================================================================
# STEP 5: VARIANT CALLING
# ==============================================================================
echo ""
echo "Step 5: Variant Calling"
echo "======================="

echo "Calling variants with Bcftools..."

# Call variants using bcftools mpileup and call
bcftools mpileup -Ou -f salmonella_ref.fasta reads_sorted.bam | \
  bcftools call -mv -Ob -o reads_variants.bcf

# Convert BCF to VCF for easier inspection
bcftools view reads_variants.bcf > reads_variants.vcf

echo "Variant calling complete!"
echo ""
echo "Variant Statistics:"

# Count total variants
TOTAL_VARIANTS=$(bcftools view -H reads_variants.vcf | wc -l)
echo "Total variants: $TOTAL_VARIANTS"

# Count SNPs
SNP_COUNT=$(bcftools view -H reads_variants.vcf | \
  awk '{if(length($4)==1 && length($5)==1) print}' | wc -l)
echo "SNPs: $SNP_COUNT"

# Count insertions
INS_COUNT=$(bcftools view -H reads_variants.vcf | \
  awk '{if(length($4)<length($5)) print}' | wc -l)
echo "Insertions: $INS_COUNT"

# Count deletions
DEL_COUNT=$(bcftools view -H reads_variants.vcf | \
  awk '{if(length($4)>length($5)) print}' | wc -l)
echo "Deletions: $DEL_COUNT"

# Extract variant data for visualization
bcftools view -H reads_variants.vcf | \
  awk '{print $1, $2, length($4), length($5), $6}' > variant_data.txt

# Extract genome-wide depth for coverage visualization
samtools depth reads_sorted.bam > genome_depth.txt

# ==============================================================================
# STEP 6: DATA VISUALIZATION
# ==============================================================================
echo ""
echo "Step 6: Generating Visualizations"
echo "=================================="

echo "Creating Figure 1: Assembly Quality Assessment..."
Rscript scripts/figure1_assembly.R

echo "Creating Figure 2: Reference Alignment Quality..."
Rscript scripts/figure2_alignment.R

echo "Creating Figure 3: Variant Landscape Analysis..."
Rscript scripts/figure3_variants.R

echo "All figures generated successfully!"

# ==============================================================================
# STEP 7: IGV VISUALIZATION (MANUAL STEP)
# ==============================================================================
echo ""
echo "Step 7: IGV Visualization (Manual)"
echo "=================================="
echo "To visualize alignments and variants in IGV:"
echo "1. Download and install IGV: https://software.broadinstitute.org/software/igv/"
echo "2. Load genome: File → Load Genome from File → salmonella_ref.fasta"
echo "3. Load alignment: File → Load from File → reads_sorted.bam"
echo "4. Load variants: File → Load from File → reads_variants.vcf"
echo "5. Navigate to regions of interest and take screenshots"

# ==============================================================================
# ANALYSIS COMPLETE
# ==============================================================================
echo ""
echo "======================================"
echo "ANALYSIS PIPELINE COMPLETE!"
echo "======================================"
echo ""
echo "Generated Files:"
echo "  Assembly:"
echo "    - flye_output/assembly.fasta (3 contigs, 5.10 Mb total)"
echo "    - flye_output/assembly_info.txt"
echo ""
echo "  Alignment:"
echo "    - reads_sorted.bam (indexed alignment file)"
echo "    - detailed_coverage.txt (coverage statistics)"
echo ""
echo "  Variants:"
echo "    - reads_variants.vcf (11,465 variants)"
echo "    - variant_data.txt (processed variant data)"
echo ""
echo "  Figures:"
echo "    - figures/figure1_assembly_final.png"
echo "    - figures/figure2_alignment_final.png"
echo "    - figures/figure3_variants_final.png"
echo ""
echo "Key Results:"
echo "  - Assembly: 3 contigs, N50 = 3.32 Mb, mean coverage = 160×"
echo "  - Alignment: 97.8% chromosome coverage, 151× depth"
echo "  - Variants: $TOTAL_VARIANTS total ($SNP_COUNT SNPs, $INS_COUNT insertions, $DEL_COUNT deletions)"
echo "  - Plasmid variant density: 87× higher than chromosome"
echo ""
echo "For full analysis and interpretation, see:"
echo "  BINF6110_Assignment1_Part2_FINAL.md"
echo ""
echo "GitHub Repository:"
echo "  https://github.com/emantahir1/binf6110-assignment1-salmonella-assembly"
echo ""
echo "======================================"
