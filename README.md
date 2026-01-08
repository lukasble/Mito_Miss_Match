# Dog Aging and Mito-Nuclear Mismatch

This repository contains all code used for the study. Below, the files belonging to each analysis are listed.  

---

## Heteroplasmy Analysis

- **Filtering to only keep variants that have passed QC (works for both raw and filtered files):** [`filtering_mito_VCFs.sh`](./filtering_mito_VCFs.sh)  
- **Filtering with initial threshold 4%:** [`4_threshold_heteroplasmy.sh`](./4_threshold_heteroplasmy.sh)  
- **Listing and plotting histograms (works for the different threshold):**  
  - [`listing_plotting_heteroplasmy.py`](./listing_plotting_heteroplasmy.py)  
  - [`listing_plotting_heteroplasmy.sh`](./listing_plotting_heteroplasmy.sh)  

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
