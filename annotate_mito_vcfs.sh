[lubl5753@rackham3 SNPeff]$ cat annotate_mito_vcfs.sh
#!/bin/bash
#SBATCH -A uppmax2025-2-119
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -J snpeff_mito
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/lukas/SNPeff/snpeff_mito_%j.out
#SBATCH -e /proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/lukas/SNPeff/snpeff_mito_%j.err

module load bioinfo-tools
module load snpEff
module load bcftools

# ---- paths ----
PROJ=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project

# Custom SnpEff DB directory + config
SNPEFF_DB_DIR=$PROJ/SNPEff/canFam4_Y.NCBIAnnotation_SNPeffV3
SNPEFF_CONFIG=$SNPEFF_DB_DIR/snpEff.config
GENOME_ID=canFam4_y    # <- from 'canFam4_y.genome: canFam4_y' in snpEff.config

# Input VCFs (already filtered mito nuclear + no coyotes)
VCF_MAIN=$PROJ/lukas/vcf_filtering/final_dog_mito_vcf/dog_mitonuclear.noCoyotes.vcf.gz
VCF_CHRX=$PROJ/lukas/vcf_filtering_chrx/final_dog_mito_vcf_chrx/chrx_dog_mitonuclear.noCoyotes.vcf.gz

# Output directory for SnpEff results
OUTDIR=$PROJ/lukas/SNPeff
mkdir -p "$OUTDIR"

OUT_MAIN=$OUTDIR/dog_mitonuclear.noCoyotes.snpeff.vcf.gz
OUT_CHRX=$OUTDIR/chrx_dog_mitonuclear.noCoyotes.snpeff.vcf.gz

echo "Using SnpEff config: $SNPEFF_CONFIG"
echo "Genome ID:          $GENOME_ID"
echo "Data dir:           $SNPEFF_DB_DIR"
echo "Output directory:   $OUTDIR"
echo "Started at:         $(date)"

cd "$SNPEFF_DB_DIR"

# ---- annotate autosomes + XPAR ----
echo "Annotating main mito VCF..."
snpEff -c "$SNPEFF_CONFIG" -v -dataDir "$SNPEFF_DB_DIR" "$GENOME_ID" "$VCF_MAIN" \
  | bcftools view -Oz -o "$OUT_MAIN"

bcftools index -f "$OUT_MAIN"

# ---- annotate chrX NONPAR ----
echo "Annotating chrX NONPAR mito VCF..."
snpEff -c "$SNPEFF_CONFIG" -v -dataDir "$SNPEFF_DB_DIR" "$GENOME_ID" "$VCF_CHRX" \
  | bcftools view -Oz -o "$OUT_CHRX"

bcftools index -f "$OUT_CHRX"

echo "Finished SnpEff at: $(date)"

# ---- quick sanity checks ----
echo "Main VCF sample count:"
bcftools query -l "$OUT_MAIN" | wc -l

echo "chrX VCF sample count:"
bcftools query -l "$OUT_CHRX" | wc -l

echo "Check ANN header (main):"
bcftools view "$OUT_MAIN" | head -40 | grep -i ANN || echo "No ANN header found?"

echo "First annotated variant (main):"
bcftools view -H "$OUT_MAIN" | head -1
