#!/bin/bash
#SBATCH -A uppmax2025-2-452
#SBATCH -p pelle
#SBATCH -M pelle
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --job-name=plot_mito_weight
#SBATCH --output=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/out_err/plot_mito_weight_%j.out
#SBATCH --error=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/out_err/plot_mito_weight_%j.err

module purge
module load BCFtools/1.22-GCC-13.3.0
module load GCC/13.3.0
module load HTSlib/1.22-GCC-13.3.0
module load BCFtools/1.22-GCC-13.3.0
module load Python/3.11.5-GCCcore-13.3.0



python plot_deleterious_vs_weight.py
