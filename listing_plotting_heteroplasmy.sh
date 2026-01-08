#!/bin/bash
#SBATCH -A uppmax2025-2-119
#SBATCH -p pelle
#SBATCH -M pelle
#SBATCH -n 12
#SBATCH -t 12:00:00
#SBATCH --output=mito_pipeline.out
#SBATCH --error=mito_pipeline.err

# ------------------------
# Downloading needed modules
# ------------------------
module load GCC/13.3.0
module load HTSlib/1.22-GCC-13.3.0
module load BCFtools/1.22-GCC-13.3.0
module load Python/3.11.5-GCCcore-13.3.0

# ------------------------
# Paths and directories
# ------------------------
VCF_DIR="/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/1_results_raw_files/filtered_PASS"
OUTDIR="/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/5_results_raw_files"
mkdir -p $OUTDIR

LIST="$OUTDIR/vcf_list.txt"
MERGED="$OUTDIR/merged_mito.vcf.gz"
EXTRACT="$OUTDIR/mito_AD_AF.tsv"
PYTHON_SCRIPT="$OUTDIR/analyze_mito_AF.py"

# ------------------------
# Setting thresholds
# ------------------------
MIN_DP=50
MIN_ALT=10
AF_MIN=0.02
AF_MAX=0.95

# ------------------------
# Listing the VCF files
# ------------------------
ls ${VCF_DIR}/*.vcf.gz > $LIST
echo "Found $(wc -l < $LIST) VCFs"

# ------------------------
# Fixing malinformed INFO fiels
# ------------------------
mkdir -p $OUTDIR/cleaned_vcfs
for v in ${VCF_DIR}/*.vcf.gz; do
    base=$(basename $v)
    bcftools annotate -x INFO/AS_FilterStatus \
        -O z -o $OUTDIR/cleaned_vcfs/$base $v
    tabix -p vcf $OUTDIR/cleaned_vcfs/$base
done
ls $OUTDIR/cleaned_vcfs/*.vcf.gz > $LIST

# ------------------------
# Merging
# ------------------------
bcftools merge --file-list $LIST \
    --missing-to-ref \
    --force-samples \
    -Oz \
    -o $MERGED
tabix -p vcf $MERGED


# ------------------------
# Extracting AD/DP/AF
# ------------------------
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE\t%AD\t%DP\t%AF]\n' \
    $MERGED > $EXTRACT

# ------------------------
# Starting the python analysis script
# ------------------------
cat <<EOF > $PYTHON_SCRIPT
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict

# ------------------------
# Setting thresholds
# ------------------------
MIN_ALT = 10
MIN_DP  = 50
AF_MIN  = 0.02
AF_MAX  = 0.95

infile = "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/5_results_raw_files/mito_AD_AF.tsv"
afs = []
burden = defaultdict(int)
het_list = []


# ------------------------
# Getting all sample IDs
# ------------------------
all_samples = set()

with open(infile) as f:
    for line in f:
        p = line.strip().split("\t")
        if len(p) < 8:
            continue
        for i in range(4, len(p), 4):
            all_samples.add(p[i])


# ------------------------
# Processing heteroplasmies
# ------------------------
with open(infile) as f:
    for line in f:
        p = line.strip().split("\t")
        if len(p) < 8:
            continue

        chrom, pos, ref, alt = p[0], p[1], p[2], p[3]

        for i in range(4, len(p), 4):
            sample, ad, dp, af = p[i], p[i+1], p[i+2], p[i+3]

            if ad == ".":
                continue
            try:
                ads = list(map(int, ad.split(",")))
            except:
                continue

            ref_ct = ads[0]
            alt_cts = ads[1:]
            alt_sum = sum(alt_cts)
            dp_val = sum(ads)

            if dp_val < MIN_DP:
                continue

            af_val = alt_sum / dp_val if dp_val > 0 else 0

            if alt_sum >= MIN_ALT:
                afs.append(af_val)

            if (alt_sum >= MIN_ALT) and (AF_MIN <= af_val <= AF_MAX):
                burden[sample] += 1
                het_list.append({
                    "Sample": sample,
                    "Chrom": chrom,
                    "Pos": pos,
                    "Ref": ref,
                    "Alt": alt,
                    "DP": dp_val,
                    "ALT_Count": alt_sum,
                    "AF": round(af_val, 5)
                })

# ------------------------------
# Saving burden table
# ------------------------------
full_burden = {s: burden.get(s, 0) for s in sorted(all_samples)}

burden_df = pd.DataFrame({
    "Sample": list(full_burden.keys()),
    "Heteroplasmy_Count": list(full_burden.values())
})

burden_df.to_csv("/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/5_results_raw_files/heteroplasmy_burden.tsv",
                 sep="\t", index=False)

# ------------------------------
# Saving full list of heteroplasmies
# ------------------------------
het_df = pd.DataFrame(het_list)
het_df = het_df.sort_values(by=["Sample", "Chrom", "Pos"])
het_df.to_csv("/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/5_results_raw_files/heteroplasmy_list.tsv",
              sep="\t", index=False)

# ------------------------------
# Plotting histogram
# ------------------------------
plt.figure(figsize=(10,4))
plt.hist(afs, bins=2000, edgecolor='none')
plt.xlabel("Allele fraction (AF)")
plt.ylabel("Count")
plt.title("Mitochondrial allele fraction distribution")
plt.tight_layout()
plt.savefig("/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/5_results_raw_files/AF_histogram.png", dpi=150)

# ------------------------------
# Printing bin summary
# ------------------------------
bins = [0,0.01,0.05,0.2,0.8,0.95,1.0]
labels = ["<1%","1–5%","5–20%","20–80%","80–95%","95–100%"]
counts = [0]*6

for af in afs:
    for i in range(len(bins)-1):
        if bins[i] <= af < bins[i+1]:
            counts[i] += 1
            break

print("AF bin counts:")
for lab,c in zip(labels,counts):
    print(f"{lab}: {c}")

EOF


# ------------------------
# Running Python
# ------------------------

pip install --user matplotlib pandas
python3 $PYTHON_SCRIPT


echo "DONE"
echo "All output files saved in: $OUTDIR"


