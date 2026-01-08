# Dog Aging and Mito-Nuclear Mismatch

This repository contains all code used for the study. Below, the files belonging to each analysis are listed.  

---

## Heteroplasmy Analysis

*(Add the scripts and files relevant to heteroplasmy analysis here. Example placeholders:)*

- [`heteroplasmy_pipeline.sh`](./heteroplasmy_pipeline.sh) – Main pipeline to call mitochondrial variants.  
- [`filter_heteroplasmy.py`](./filter_heteroplasmy.py) – Script to filter heteroplasmic sites.  
- [`plot_heteroplasmy_distribution.R`](./plot_heteroplasmy_distribution.R) – Generates plots of heteroplasmy levels across samples.  

---

## Mitochondrial Genetic Load

- **Creating a custom snpEff database:** [`snpEff_mito_config`](./snpEff_mito_config)  
- **Running snpEff on all VCFs:** [`snpEff_mito_run.sh`](./snpEff_mito_run.sh)  
- **Listing all samples with the count of deleterious variants:** [`extract_deleterious_mito.sh`](./extract_deleterious_mito.sh)  
- **Plotting deleterious variant count vs body weight:**  
  - [`plot_geneticload_vs_weights.py`](./plot_geneticload_vs_weights.py)  
  - [`plot_geneticload_vs_weights.sh`](./plot_geneticload_vs_weights.sh)  
- **Plotting deleterious variants in three categories:**  
  - [`plot_groups_mito_geneticload.py`](./plot_groups_mito_geneticload.py)  
  - [`plot_groups_mito_geneticload.sh`](./plot_groups_mito_geneticload.sh)  
