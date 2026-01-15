#!/bin/bash
#SBATCH -A uppmax2025-2-119
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH -J genetic_load_bcf
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/lukas/genetic_load/genetic_load_bcf_%j.out
#SBATCH -e /proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/lukas/genetic_load/genetic_load_bcf_%j.err

module load bioinfo-tools
module load bcftools

PROJ=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project

MAIN_SNPEFF=$PROJ/lukas/SNPeff/dog_mitonuclear.noCoyotes.snpeff.vcf.gz
CHRX_SNPEFF=$PROJ/lukas/SNPeff/chrx_dog_mitonuclear.noCoyotes.snpeff.vcf.gz

OUTDIR=$PROJ/lukas/genetic_load
mkdir -p "$OUTDIR"

AUTOSOME_DMG=$OUTDIR/autosomes_damaging.vcf.gz
CHRX_DMG=$OUTDIR/chrx_damaging.vcf.gz
MERGED_DMG=$OUTDIR/damaging_merged.vcf.gz
SAMPLES=$OUTDIR/samples.txt
GT_TABLE=$OUTDIR/genotypes_for_load.tsv

echo "Started bcftools-only job at: $(date)"
echo "Main SNPeff VCF:  $MAIN_SNPEFF"
echo "chrX SNPeff VCF:  $CHRX_SNPEFF"
echo "Output dir:       $OUTDIR"
echo

echo "Filtering autosomes + XPAR for HIGH/MODERATE impact..."
bcftools view \
    -i 'INFO/ANN ~ "HIGH" || INFO/ANN ~ "MODERATE"' \
    "$MAIN_SNPEFF" \
    -Oz -o "$AUTOSOME_DMG"

bcftools index -f "$AUTOSOME_DMG"
echo "Number of autosome+XPAR damaging variants:"
bcftools view -H "$AUTOSOME_DMG" | wc -l
echo

echo "Filtering chrX NONPAR for HIGH/MODERATE impact..."
bcftools view \
    -i 'INFO/ANN ~ "HIGH" || INFO/ANN ~ "MODERATE"' \
    "$CHRX_SNPEFF" \
    -Oz -o "$CHRX_DMG"

bcftools index -f "$CHRX_DMG"
echo "Number of chrX NONPAR damaging variants:"
bcftools view -H "$CHRX_DMG" | wc -l
echo

echo "Concatenating damaging variants (autosomes + chrX NONPAR)..."
bcftools concat -a \
    "$AUTOSOME_DMG" \
    "$CHRX_DMG" \
    -Oz -o "$MERGED_DMG"

bcftools index -f "$MERGED_DMG"

echo "Merged damaging VCF: $MERGED_DMG"
echo "Number of merged damaging variants:"
bcftools view -H "$MERGED_DMG" | wc -l
echo

echo "Extracting sample IDs..."
bcftools query -l "$MERGED_DMG" > "$SAMPLES"
NSAMP=$(wc -l < "$SAMPLES")
echo "Number of samples: $NSAMP"
echo

echo "Creating genotype table (this may take a bit)..."
bcftools query -f '%CHROM\t%POS[\t%GT]\n' "$MERGED_DMG" > "$GT_TABLE"

echo "Genotype table written to: $GT_TABLE"
echo "Number of lines (variants) in genotype table:"
wc -l "$GT_TABLE"

echo "bcftools-only job finished at: $(date)"
