#Code for filtering out variants that have not passed the quality control. 

#!/bin/bash
#SBATCH -A uppmax2025-2-119      
#SBATCH -p pelle                  
#SBATCH -M pelle                 
#SBATCH -n 8                      
#SBATCH -t 12:00:00               
#SBATCH -J mito_filtered          
#SBATCH -o mito_filter_%j.out     
#SBATCH -e mito_filter_%j.err     
#SBATCH --array=1-1972            


module load BCFtools/1.22-GCC-13.3.0

VCF_DIR="/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/mg_VCFs"
OUT_DIR="/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/filtered_mg_VCFs"
mkdir -p $OUT_DIR

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

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\n' /proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/filtered_mg_VCFs/ACKR000003_filtered.vcf.gz | head

