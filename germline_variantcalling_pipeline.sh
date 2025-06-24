#!/bin/bash

# Germline Variant Calling Pipeline for SRRYYYYYYY using GATK Haplotypecaller

# Step 0: Setup directories
mkdir -p fastqc_output trimmed_reads

# Step 1: Quality Control
fastqc SRRXXXXXXX_1.fastq.gz SRRXXXXXXX_2.fastq.gz -o fastqc_output

# Step 2: Trimming with Trimmomatic
trimmomatic PE -threads 4 \
  SRRXXXXXXX_1.fastq.gz SRRXXXXXXX_2.fastq.gz \
  trimmed_reads/forward_paired.fq.gz trimmed_reads/forward_unpaired.fq.gz \
  trimmed_reads/reverse_paired.fq.gz trimmed_reads/reverse_unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36

# Step 3: Reference Genome Setup (run only once)
wget -c https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
gatk CreateSequenceDictionary \
  -R Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -O Homo_sapiens.GRCh38.dna.primary_assembly.dict

# Step 4: Alignment
bwa mem -t 8 -R "@RG\tID:SRRXXXXXXX\tSM:SRRXXXXXXX\tPL:ILLUMINA" \
  Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  SRRXXXXXXX_1.fastq.gz SRRXXXXXXX_2.fastq.gz > SRRXXXXXXX.sam

samtools view -Sb SRRXXXXXXX.sam > SRRXXXXXXX.bam
samtools sort SRRXXXXXXX.bam -o SRRXXXXXXX_sorted.bam
samtools index SRRXXXXXXX_sorted.bam

# Step 5: Mark Duplicates
gatk MarkDuplicatesSpark \
  -I SRRXXXXXXX_sorted.bam \
  -O SRRXXXXXXX_sorted_dedup.bam \
  -M SRRXXXXXXX_dedup_metrics.txt
samtools index SRRXXXXXXX_sorted_dedup.bam

# Step 6: Germline Variant Calling
gatk HaplotypeCaller \
  -R Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -I SRRXXXXXXX_sorted_dedup.bam \
  -O SRRXXXXXXX_raw_germline.vcf

# Step 7: Filtering Variants
REF=../Homo_sapiens.GRCh38.dna.primary_assembly.fa
INPUT_VCF=SRRXXXXXXX_raw_germline.vcf

# SNPs
gatk SelectVariants -R $REF -V $INPUT_VCF --select-type-to-include SNP -O SRRXXXXXXX_raw_snps.vcf
gatk VariantFiltration -R $REF -V SRRXXXXXXX_raw_snps.vcf -O SRRXXXXXXX_filtered_snps.vcf \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filter-name "snp_filter"

# INDELs
gatk SelectVariants -R $REF -V $INPUT_VCF --select-type-to-include INDEL -O SRRXXXXXXX_raw_indels.vcf
gatk VariantFiltration -R $REF -V SRRXXXXXXX_raw_indels.vcf -O SRRXXXXXXX_filtered_indels.vcf \
  --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" \
  --filter-name "indel_filter"

# Merge filtered SNPs and INDELs
gatk MergeVcfs \
  -I SRRXXXXXXX_filtered_snps.vcf \
  -I SRRXXXXXXX_filtered_indels.vcf \
  -O SRRXXXXXXX_filtered_variants.vcf

# Step 8: Annotation with VEP
vep_install -a c -s homo_sapiens -y GRCh38
vep -i SRRXXXXXXX_filtered_variants.vcf \
  -o SRRXXXXXXX_annotated.vcf \
  --vcf --offline --cache \
  --species homo_sapiens \
  --assembly GRCh38 \
  --force_overwrite

# Step 9: Convert Annotated VCF to Table
gatk VariantsToTable \
  -V SRRXXXXXXX_annotated.vcf \
  -F CHROM -F POS -F REF -F ALT -F AC -F AN -F DP -F AF -F ANN \
  -O SRRXXXXXXX_variants.table
