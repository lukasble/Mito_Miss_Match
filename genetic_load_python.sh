#!/bin/bash
#SBATCH -A uppmax2025-2-119
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -J genetic_load_py
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/lukas/genetic_load/python_run/genetic_load_py_%j.out
#SBATCH -e /proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/lukas/genetic_load/python_run/genetic_load_py_%j.err

module load bioinfo-tools
module load python/3.9.5

# --- Directories ---
BASEDIR=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/lukas/genetic_load
OUTDIR=$BASEDIR/python_run
mkdir -p "$OUTDIR"

# --- Input files copied from the bcftools run ---
cp "$BASEDIR/samples.txt" "$OUTDIR/samples.txt"
cp "$BASEDIR/genotypes_for_load.tsv" "$OUTDIR/genotypes_for_load.tsv"

SAMPLES=$OUTDIR/samples.txt
GENO_TSV=$OUTDIR/genotypes_for_load.tsv
LOAD_TABLE=$OUTDIR/genetic_load_per_sample.tsv

echo "Started Python-only genetic load job at: $(date)"
echo "OUTDIR:       $OUTDIR"
echo "Samples:      $SAMPLES"
echo "Genotypes:    $GENO_TSV"
echo "Output table: $LOAD_TABLE"
echo

echo "File checks:"
wc -l "$SAMPLES"
wc -l "$GENO_TSV"
ls -lh "$GENO_TSV"
echo

export OUTDIR

python3 - << 'EOF'
import os
import sys

outdir = os.environ["OUTDIR"]
samples_file = os.path.join(outdir, "samples.txt")
geno_file = os.path.join(outdir, "genotypes_for_load.tsv")
load_out = os.path.join(outdir, "genetic_load_per_sample.tsv")

# Load samples
with open(samples_file) as f:
    samples = [line.strip() for line in f if line.strip()]

n = len(samples)
loads = [0] * n

sys.stderr.write(f"Loaded {n} samples.\n")
sys.stderr.write(f"Reading genotypes from: {geno_file}\n")

line_count = 0

with open(geno_file) as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        gts = parts[2:]

        if len(gts) != n:
            sys.stderr.write(f"WARNING: mismatch #GT={len(gts)} vs samples={n} at variant {line_count+1}\n")

        for i, gt in enumerate(gts):
            if i >= n:
                break
            if gt in ("./.", ".|.", ".", ""):
                continue
            alleles = gt.replace("|", "/").split("/")
            alt_count = sum(1 for a in alleles if a not in ("0", "."))
            loads[i] += alt_count

        line_count += 1
        if line_count % 10000 == 0:
            sys.stderr.write(f"Processed {line_count} variants...\n")

sys.stderr.write(f"Finished processing {line_count} variants.\n")

with open(load_out, "w") as out:
    out.write("Sample\tDamagingAltAlleleCount\n")
    for s, l in zip(samples, loads):
        out.write(f"{s}\t{l}\n")

sys.stderr.write(f"Output written to: {load_out}\n")
EOF

echo
echo "Python-only genetic load finished at: $(date)"
