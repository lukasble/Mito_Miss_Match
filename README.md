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

---
    
## Mito-Nuclear Genetic Load

### Core Genetic Load Pipeline

- **Build target region + subset VCF:** [`vcf_filtering_translation.txt`](./vcf_filtering_translation.txt)  
- **Subsetting genome-wide VCF to mito-nuclear gene regions (excluding coyotes):** [`bcftools.sh`](./bcftools.sh)
- **Subsetting ChrX (non-PAR) VCF to mito-nuclear genes (excluding coyotes):** [`bcftools_chrx.sh`](./bcftools_chrx.sh)
- **Annotating mito-nuclear subset VCFs with snpEff (autosomes + XPAR and ChrX NONPAR):** [`annotate_mito_vcfs.sh`](./annotate_mito_vcfs.sh)
- **Extracting damaging mito-nuclear variants and generating genotype tables:** [`genetic_load.sh`](./genetic_load.sh) 
- **Computing per-sample genetic load from damaging variant genotypes:** [`genetic_load_python.sh`](./genetic_load_python.sh)

### Core Pathway-level Analysis and Data Examination:

- **Annotation of VCF file for Nuclear Autosome + Chrx-NonPar SNPs:** annotation_SNPeff.sh
- **Pathway-level genetic load analyses:** [`vcf_filtering_translation.txt`](./vcf_filtering_translation.txt)  
- **Filtering of annotated VCF files, data extraction and conversion (General Impact, Synonymous and Nonsynonymous Variant Analysis):**
  - **Data Extraction**:
    - [`syn_extract.sh`](./syn_extract.sh)
    - [`nonsyn_extract.sh`](./nonsyn_extract.sh)
  - **Data Conversion**:
    - [`syn_extract.sh`](./syn_extract.sh)
    - [`nonsyn_extract.sh`](./nonsyn_extract.sh)
- **Data Analysis and Visualization study of genetic variation and genetic load metrics:** [`genetic_load.sh`](./genetic_load.sh)
