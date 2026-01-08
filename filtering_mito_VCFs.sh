#!/bin/bash
#SBATCH -A uppmax2025-2-452     
#SBATCH -p pelle                  
#SBATCH -M pelle                  
#SBATCH -n 1                   
#SBATCH -t 12:00:00              
#SBATCH -J mito_filtered          
#SBATCH -o mito_filter_%j.out     
#SBATCH -e mito_filter_%j.err     
#SBATCH --array=1-1972            

# ------------------------
# Downloading needed modules
# ------------------------
module load BCFtools/1.22-GCC-13.3.0


# ------------------------
# Path and directories
# ------------------------
VCF_DIR="/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/mg_VCFs"
OUT_DIR="/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/filtered_mg_VCFs"
mkdir -p $OUT_DIR


# ---------------------------------------------------------
# Filtering out all variants that have not passed filters
# ---------------------------------------------------------
FILES=($VCF_DIR/*.vcf.gz)
VCF=${FILES[$SLURM_ARRAY_TASK_ID-1]}

BASE=$(basename $VCF .mitoMerged.vcf.gz)

bcftools view \
    -i 'ALT!="" && FMT/AD!="."' \
    -f PASS \
    -O z \
    -o $OUT_DIR/${BASE}_filtered.vcf.gz \
    $VCF

bcftools index $OUT_DIR/${BASE}_filtered.vcf.gz
