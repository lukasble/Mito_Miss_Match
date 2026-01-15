#!/bin/bash
#SBATCH -A uppmax2025-2-119      # <-- use this account
#SBATCH -p core                  # core partition on Rackham
#SBATCH -n 2
#SBATCH -t 05:00:00
#SBATCH -J dog_mito_vcf
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se

module load bioinfo-tools
module load bcftools

VCF=../../../data/ng_VCFs/AutoAndXPAR.SNPs.vcf.vqsr99.gz
BED=../filter_dog_mito_genes/mito_nuclear_genes_final_tab.bed

bcftools view -S ^coyote.txt -R "$BED" -Oz "$VCF" \
  -o dog_mitonuclear.noCoyotes.vcf.gz

bcftools index dog_mitonuclear.noCoyotes.vcf.gz
