#!/bin/bash
#SBATCH -A uppmax2025-2-452      # active project
#SBATCH -p pelle                  # valid partition
#SBATCH -M pelle                  # cluster
#SBATCH -n 8                      # number of cores
#SBATCH -t 12:00:00               # walltime
#SBATCH -J mito_filtered          # job name
#SBATCH -o mito_filter_%j.out     # stdout
#SBATCH -e mito_filter_%j.err     # stderr
#SBATCH --array=1-1972            # array jobs for all files


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
