#!/bin/bash
#SBATCH -A uppmax2025-2-452 
#SBATCH -p pelle
#SBATCH -M pelle
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --job-name=plot_mito_weight
#SBATCH --output=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/out_err/plot_mito_box_%j.out
#SBATCH --error=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/out_err/plot_mito_box_%j.err

module purge
module load Python/3.11.5-GCCcore-13.3.0
module load BCFtools/1.22-GCC-13.3.0
module load BCFtools

python plot_categories_mito_deleterious_boxplots.py
