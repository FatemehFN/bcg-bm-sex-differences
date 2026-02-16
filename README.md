# Sex-Biased Gene Expression Analysis in Human Bone Marrow

## Overview

This repository contains a comprehensive analysis of sex differences in CD34+ bone marrow scRNA-seq data from Sun et al. (Immunity 2024; GEO accessions GSE248728/30). The analysis identifies sex-biased genes, modules, and sample-level trajectories in the context of BCG vaccination response.

## Key Findings

- **Global sex effects**: Identified significant sex-associated genes (FDR < 0.05) using both OLS and mixed-effects models
- **Sex-stratified vaccination response**: Males and females show differential transcriptional responses to BCG vaccination
- **Gene co-expression modules**: Built KNN-based gene networks to identify modules with sex-biased activity
- **Trajectory analysis**: Females occupy a more "advanced" state along the PC1-defined transcriptional trajectory

## Repository Structure

```
.
├── sex_diff_analysis.ipynb           # Main analysis notebook
├── donor_sex.csv                     # Inferred donor sex labels (M/F)
├── donor_sex_inference_detailed.csv  # Detailed sex inference metrics
├── infer_donor_sex_from_scRNA.py     # Script to infer donor sex from gene expression
├── simple_sex_inference.py           # Simplified sex inference approach
├── make_scRNA_metadata_summary.py    # Cell metadata aggregation script
├── download_captures.py              # GEO data download utility
├── scRNA_donor_timepoint_summary.csv # Per-donor/timepoint summary statistics
├── analysis_outputs/                 # Generated results
│   ├── donor_summary_with_sex.csv
│   ├── capture_summary_with_sex.csv
│   ├── global_sex_effects_log2cpm.csv
│   ├── global_sex_effects_log2cpm_mixedLM.csv
│   ├── pseudobulk_pb_filt_counts.csv
│   ├── pseudobulk_pb_log2cpm.csv
│   └── sex_biased_module_*_genes.txt
├── BCG-humanBM-code/                 # Original R analysis pipeline
└── GSE248728_downloads/              # Raw 10x matrices (14 captures)
```

## Analysis Workflow

### 1. Data Processing

- **Input**: Raw 10x matrices from GSE248728 (14 scRNA-seq captures)
- **QC filters**: 
  - `percent.mt < 15%`
  - `nCount_RNA > 500`
- **Pseudobulk aggregation**: Counts summed per sample (UNIQUE_ID)
- **Gene filtering**: Mean CPM > 1 across samples
- **Normalization**: Log2(CPM + 1) transformation

### 2. Sex Inference

Donor sex was inferred using expression of sex-linked marker genes:
- **Female markers**: XIST, TSIX
- **Male markers**: RPS4Y1, RPS4Y2, DDX3Y, KDM5D, UTY, USP9Y, EIF1AY, ZFY, TXLNGY, TMSB4Y

See [`infer_donor_sex_from_scRNA.py`](infer_donor_sex_from_scRNA.py) for implementation details.

### 3. Statistical Models

#### Simple linear models (OLS)
```
log2CPM ~ group + timepoint + sex
```

#### Mixed-effects models
```
log2CPM ~ group + timepoint + sex + (1 | donor)
```
- Random intercept per donor to account for repeated measures
- Fitted with statsmodels MixedLM (method='lbfgs', REML=False)

#### Sex-stratified vaccination models
```
log2CPM ~ group + timepoint + group:timepoint + (1 | donor)
```
Stratified by sex to identify sex-specific BCG responses.

### 4. Network Analysis

- **Gene-gene KNN graph**: Built using top 2000 most variable genes
- **Clustering**: Leiden community detection (resolution=1.0)
- **Module eigengenes**: Mean z-scored expression per module
- **Testing**: Welch's t-test for M vs F module activity

### 5. Trajectory Analysis

- **PC1 trajectory**: Normalized PC1 from pseudobulk log2CPM PCA
- **Differential abundance**: Logistic regression of sex on PC1
- **Result**: Females enriched at later positions (β < 0, p ≈ 0.02)

## Key Results Files

| File | Description |
|------|-------------|
| `global_sex_effects_log2cpm.csv` | Per-gene sex coefficients, p-values, FDR (OLS model) |
| `global_sex_effects_log2cpm_mixedLM.csv` | Per-gene sex coefficients with donor random effects |
| `pseudobulk_pb_log2cpm.csv` | Filtered log2CPM matrix (genes × samples) |
| `pseudobulk_pb_filt_counts.csv` | Filtered count matrix (genes × samples) |
| `donor_summary_with_sex.csv` | Donor-level statistics with sex labels |
| `sex_biased_module_*_genes.txt` | Gene lists for sex-associated modules |

## Requirements

### Python packages
```bash
pip install numpy pandas matplotlib seaborn scipy scikit-learn statsmodels python-igraph leidenalg
```

### R packages (for BCG-humanBM-code comparison)
```r
install.packages(c("Seurat", "limma", "emmreml", "lme4"))
```

## Usage

### Run the main analysis
```bash
jupyter notebook sex_diff_analysis.ipynb
```

### Infer donor sex from scratch
```bash
python infer_donor_sex_from_scRNA.py
```

### Download GEO data
```bash
python download_captures.py
```

## Data Sources

- **Primary data**: Sun et al., Immunity 2024
- **GEO accessions**: GSE248728 (scRNA-seq), GSE248730 (scATAC-seq)
- **Original code**: [BCG-humanBM-code](https://github.com/sarahjiesun/BCG-humanBM-code)

## Methods Summary

1. **Sample alignment**: Paired donors (present at both Td0 and Tm3) used for modeling
2. **CPM normalization**: Library size normalization with log2 transformation
3. **Multiple testing correction**: Benjamini-Hochberg FDR control
4. **Visualization**: Volcano plots, t-SNE embeddings, trajectory plots
5. **Module detection**: KNN + Leiden clustering on gene correlation network

## Interpretation

### Sex-biased genes
- Y-chromosome genes (e.g., KDM5D, RPS4Y1) show expected male bias
- X-chromosome escapees (e.g., XIST) show expected female bias
- Autosomal genes with sex bias may reflect hormonal or immune differences

### Vaccination response
- Sex-stratified models reveal differential BCG responses
- Interaction effects quantify how vaccination timing affects males vs females differently

### Trajectory positioning
- PC1 captures major axis of transcriptional variation
- Female samples skew toward "activated" end of trajectory
- Suggests sex-dependent bone marrow state differences

## Limitations

- Sex inferred from expression (not verified against metadata)
- Single study/cohort (Sun et al. 2024)
- CD34+ cells only (not whole bone marrow)
- Cross-sectional timepoints (limited longitudinal resolution)

## Author

Analysis performed by Fatemeh Fatemi for computational biology interview.

## References

Sun et al. (2024). "BCG vaccination induces long-term epigenetic remodeling in human bone marrow hematopoietic stem cells." *Immunity*.

## License

Analysis code is provided for review purposes. Original data subject to GEO terms of use.
