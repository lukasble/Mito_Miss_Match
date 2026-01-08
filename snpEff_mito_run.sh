#!/bin/bash
#SBATCH -A uppmax2025-2-119
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH -t 04:00:00
#SBATCH --job-name=snpeff_mito_all
#SBATCH --output=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/out_err/snpeff_mito_all_%j.out
#SBATCH --error=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/out_err/snpeff_mito_all_%j.err

set -euo pipefail

# ------------------------
# Downloading needed modules 
# ------------------------
module purge
module load bioinfo-tools
module load snpEff/5.2
module load bcftools
module load tabix

# ------------------------
# Path and directories
# ------------------------
PROJ="/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina"

VCF_DIR="$PROJ/Results/1_results_raw_files/filtered_PASS"
RENAMED_DIR="$PROJ/Results/7_renamed_vcfs"
OUTDIR="$PROJ/Results/8_2_snpeff_annotated"
LOGDIR="$OUTDIR/logs"

CONFIG="$PROJ/snpEff/snpEff.config"
DATADIR="$PROJ/snpEff/data"
GENOME="dogMito"

mkdir -p "$RENAMED_DIR" "$OUTDIR" "$LOGDIR"

echo "Started: $(date)"
echo "VCF directory: $VCF_DIR"
echo

# ------------------------
# Looping over all VCFs
# ------------------------
for VCF_IN in "$VCF_DIR"/*.vcf.gz; do
    BASENAME=$(basename "$VCF_IN" .vcf.gz)
    echo "Processing: $BASENAME"

    RENAMED_VCF="$RENAMED_DIR/${BASENAME}.renamed.vcf.gz"
    OUT_VCF="$OUTDIR/${BASENAME}.snpeff.vcf.gz"
    ERRLOG="$LOGDIR/${BASENAME}.snpeff.stderr.log"

    # ------------------------
    # Renaming chromosomes
    # ------------------------
    bcftools annotate \
      --rename-chrs <(echo -e "NC_002008.4\tCM022001.1") \
      -Oz -o "$RENAMED_VCF" \
      "$VCF_IN"

    tabix -p vcf "$RENAMED_VCF"

    # ------------------------
    # Annotating with snpEff
    # ------------------------
    snpEff \
      -c "$CONFIG" \
      -dataDir "$DATADIR" \
      -noStats \
      "$GENOME" \
      "$RENAMED_VCF" 2> "$ERRLOG" \
    | bcftools view -Oz -o "$OUT_VCF"

    bcftools index -f "$OUT_VCF"

    # ------------------------
    # Sanity check
    # ------------------------
    bcftools view -h "$OUT_VCF" | grep -q "##INFO=<ID=ANN" \
      || echo "WARNING: ANN header missing in $BASENAME"

    echo "Finished $BASENAME"
    echo
done

echo "All files completed: $(date)"
