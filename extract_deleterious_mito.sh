#!/bin/bash
#SBATCH -A uppmax2025-2-119
#SBATCH -p pelle
#SBATCH -M pelle
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --job-name=deleterious_mtDNA
#SBATCH --output=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/out_err/deleterious_%j.out
#SBATCH --error=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/out_err/deleterious_%j.err

# ------------------------
# Downloading needed modules
# ------------------------
module load BCFtools/1.22-GCC-13.3.0

# ------------------------
# Paths an directories
# ------------------------
PROJ="/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina"
VCF_DIR="$PROJ/Results/8_2_snpeff_annotated"
OUT_TSV="$PROJ/Results/9_deleterious_variants/mito_deleterious_counts.tsv"

# ------------------------
# Creating header
# ------------------------
echo -e "dog_id\tdeleterious_count" > "$OUT_TSV"

# ------------------------------------------------------
# Looping over all VCFs and counting deleterious variants
# ------------------------------------------------------
for VCF in "$VCF_DIR"/*.snpeff.vcf.gz; do
    BASENAME=$(basename "$VCF")

    # Extract dog ID (everything before first dot)
    DOG_ID=${BASENAME%%.*}

    COUNT=$(bcftools view \
        -i 'INFO/ANN~"HIGH" || INFO/ANN~"MODERATE"' \
        "$VCF" \
        | bcftools view -H \
        | wc -l)

    echo -e "${DOG_ID}\t${COUNT}" >> "$OUT_TSV"
done

echo "Finished deleterious counts: $(date)"
echo "Output written to: $OUT_TSV"
