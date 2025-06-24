#!/bin/bash

# Somatic Variant Calling Pipeline for SRRYYYYYYY using GATK Mutect2

# Step 0: Setup directories
mkdir -p fastqc_output trimmed_reads

# Step 1: Quality Control
fastqc SRRYYYYYYY_1.fastq.gz SRRYYYYYYY_2.fastq.gz -o fastqc_output

# Step 2: Trimming
trimmomatic PE -threads 4 \
  SRRYYYYYYY_1.fastq.gz SRRYYYYYYY_2.fastq.gz \
  trimmed_reads/forward_paired.fq.gz trimmed_reads/forward_unpaired.fq.gz \
  trimmed_reads/reverse_paired.fq.gz trimmed_reads/reverse_unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36

# Step 3: Reference Genome Setup (skip if already done)
wget -c https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
gatk CreateSequenceDictionary \
  -R Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -O Homo_sapiens.GRCh38.dna.primary_assembly.dict

# Step 4: Alignment
bwa mem -t 8 -R "@RG\tID:SRRYYYYYYY\tSM:SRRYYYYYYY\tPL:ILLUMINA" \
  Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  SRRYYYYYYY_1.fastq.gz SRRYYYYYYY_2.fastq.gz > SRRYYYYYYY.sam

samtools view -Sb SRRYYYYYYY.sam > SRRYYYYYYY.bam
samtools sort SRRYYYYYYY.bam -o SRRYYYYYYY_sorted.bam
samtools index SRRYYYYYYY_sorted.bam

# Step 5: Mark Duplicates
gatk MarkDuplicatesSpark \
  -I SRRYYYYYYY_sorted.bam \
  -O SRRYYYYYYY_sorted_dedup.bam \
  -M SRRYYYYYYY_dedup_metrics.txt
samtools index SRRYYYYYYY_sorted_dedup.bam

# Step 6: Somatic Variant Calling with Mutect2
gatk Mutect2 \
  -R Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -I SRRYYYYYYY_sorted_dedup.bam \
  -I SRRXXXXXXX_sorted_dedup.bam \
  -tumor SRRYYYYYYY \
  -normal SRRXXXXXXX \
  -O SRRYYYYYYY_somatic_raw.vcf

# Step 7: Filter Somatic Variants
gatk FilterMutectCalls \
  -R Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -V SRRYYYYYYY_somatic_raw.vcf \
  -O SRRYYYYYYY_somatic_filtered.vcf

# Step 8: Annotation with VEP
vep_install -a c -s homo_sapiens -y GRCh38
vep -i SRRYYYYYYY_somatic_filtered.vcf \
  -o SRRYYYYYYY_somatic_annotated.vcf \
  --vcf --offline --cache \
  --species homo_sapiens \
  --assembly GRCh38 \
  --force_overwrite

# Step 9: Convert to Table
gatk VariantsToTable \
  -V SRRYYYYYYY_somatic_annotated.vcf \
  -F CHROM -F POS -F REF -F ALT -F AC -F AN -F DP -F AF -F ANN \
  -O SRRYYYYYYY_somatic_variants.table
