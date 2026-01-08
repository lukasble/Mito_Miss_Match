#!/bin/bash
#SBATCH -A uppmax2025-2-119
#SBATCH -p pelle
#SBATCH -M pelle
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J filter_AF_vcf
#SBATCH -o /proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/S_Code/out_err/filter_AF_%j.out
#SBATCH -e /proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/S_Code/out_err/filter_AF_%j.err
#SBATCH --array=1-3944   # adjust based on number of files

# ------------------------
# Downloading needed modules
# ------------------------
module load BCFtools/1.22-GCC-13.3.0

# ------------------------
# Paths and directories
# ------------------------
VCF_DIR="/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/minimal_split_VCFs"
OUT_DIR="/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/4_MAF_filtered_VCFs"

mkdir -p "$OUT_DIR"
mkdir -p "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/S_Code/out_err"

# ------------------------
# Loading files
# ------------------------
FILES=($VCF_DIR/*.vcf.gz)
VCF=${FILES[$SLURM_ARRAY_TASK_ID-1]}
BASE=$(basename "$VCF" .vcf.gz)

echo "Filtering: $BASE"

# ------------------------
# Initial MAF threshold 4% 
# ------------------------
bcftools view -i 'INFO/AF>=0.04' -Oz -o "$OUT_DIR/${BASE}_AF_filtered.vcf.gz" "$VCF"

bcftools index "$OUT_DIR/${BASE}_AF_filtered.vcf.gz"

echo "Finished: $BASE"
