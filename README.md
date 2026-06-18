# GEO Gene Expression Analysis — GSE10072

> **Lung Adenocarcinoma vs. Normal Tissue**  
> A reproducible end-to-end RNA expression analysis pipeline built on public NCBI GEO data.

-----

## Table of Contents

- [Background](#background)
- [Dataset](#dataset)
- [Repository Structure](#repository-structure)
- [Pipeline Overview](#pipeline-overview)
- [Key Results](#key-results)
- [Quick Start](#quick-start)
- [Dependencies](#dependencies)
- [Methods](#methods)
- [Figures](#figures)
- [Citation](#citation)

-----

## Background

Lung adenocarcinoma (LUAD) is the most common subtype of non-small-cell lung cancer and a leading cause of cancer mortality worldwide. Identifying gene expression signatures that distinguish tumour from normal tissue can reveal biomarkers for diagnosis and potential therapeutic targets.

This project downloads, processes, and analyses the publicly available **GSE10072** dataset using a fully reproducible Python pipeline.

-----

## Dataset

|Field            |Value                                                                  |
|-----------------|-----------------------------------------------------------------------|
|**GEO Accession**|[GSE10072](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10072)|
|**Platform**     |GPL96 — Affymetrix Human Genome U133A Array                            |
|**Samples**      |107 (58 lung adenocarcinoma tumours, 49 adjacent normal tissue)        |
|**Organism**     |*Homo sapiens*                                                         |
|**Publication**  |Landi et al., *Cancer Research*, 2008                                  |

-----

## Repository Structure

```
geo-expression-analysis/
│
├── analysis.py              # Main pipeline script
├── requirements.txt         # Python dependencies
├── .gitignore
├── README.md
│
├── data/                    # Downloaded GEO files (auto-created; gitignored)
│
└── results/
    ├── figures/
    │   ├── 01_sample_boxplot.png     # Per-sample QC distribution
    │   ├── 02_pca.png                # PCA — tumour vs. normal separation
    │   ├── 03_volcano.png            # Volcano plot of DEGs
    │   └── 04_heatmap_top_DEGs.png  # Heatmap of top 50 DEGs
    └── tables/
        ├── DE_results.csv            # Full differential expression table
        └── summary.csv               # Run summary statistics
```

-----

## Pipeline Overview

```
Download GEO (GEOparse)
        │
        ▼
Parse expression matrix & metadata
        │
        ▼
Log₂ transform + low-expression filter
        │
        ▼
Quantile normalisation (between-sample)
        │
        ▼
Welch's t-test per probe (Tumour vs. Normal)
        │
        ▼
FDR correction (Benjamini-Hochberg)
        │
        ▼
Visualisations: QC boxplot, PCA, Volcano, Heatmap
        │
        ▼
Export tables (DE_results.csv, summary.csv)
```

-----

## Key Results

|Metric                                     |Value  |
|-------------------------------------------|-------|
|Probes after filtering                     |~10,000|
|Significant DEGs (FDR < 0.05, |log2FC| ≥ 1)|~1,200 |
|Up-regulated in Tumour                     |~700   |
|Down-regulated in Tumour                   |~500   |

*Exact numbers vary slightly with GEOparse version and platform availability.*

-----

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/<your-username>/geo-expression-analysis.git
cd geo-expression-analysis
```

### 2. Create a virtual environment (recommended)

```bash
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

### 4. Run the analysis

```bash
python analysis.py
```

GEO data is downloaded automatically on first run (~150 MB). Subsequent runs use the cached files in `data/`.

-----

## Dependencies

|Package       |Purpose                      |
|--------------|-----------------------------|
|`GEOparse`    |Download & parse GEO datasets|
|`pandas`      |Data manipulation            |
|`numpy`       |Numerical computation        |
|`scipy`       |Welch’s t-test               |
|`statsmodels` |FDR correction (BH method)   |
|`scikit-learn`|PCA                          |
|`matplotlib`  |Base plotting                |
|`seaborn`     |Statistical visualisations   |

See [`requirements.txt`](requirements.txt) for pinned versions.

-----

## Methods

### Normalisation

Raw probe-level signal intensities were log₂-transformed after clipping values below 1 to avoid undefined logarithms. Probes with a median log₂ expression below log₂(50) across all samples were removed as unexpressed. Remaining probes were quantile-normalised to make the distribution of expression values comparable across samples.

### Differential Expression

Differentially expressed probes were identified using **Welch’s two-sample t-test** (assuming unequal variances) comparing 58 tumour samples against 49 normal samples. Raw p-values were corrected for multiple testing using the **Benjamini-Hochberg false discovery rate (FDR)** procedure. Probes with FDR-adjusted p-value < 0.05 **and** |log₂ fold-change| ≥ 1 (i.e., ≥ 2-fold change) were considered statistically and biologically significant.

### Visualisations

- **Boxplot**: per-sample expression distributions post-normalisation — confirms successful normalisation.
- **PCA**: principal component analysis of all samples coloured by tissue type — assesses global separation.
- **Volcano plot**: log₂FC vs. −log₁₀(adjusted p-value) — overview of DEG magnitude and significance.
- **Heatmap**: row-clustered z-scored expression of the top 50 DEGs by absolute fold-change — reveals expression patterns across sample groups.

-----

## Figures

### PCA — Tumour vs. Normal Separation

![PCA](results/figures/02_pca.png)

### Volcano Plot

![Volcano](results/figures/03_volcano.png)

### Heatmap — Top 50 DEGs

![Heatmap](results/figures/04_heatmap_top_DEGs.png)

-----

## Citation

**Dataset:**  
Landi MT et al. Gene expression signature of cigarette smoking and its role in lung adenocarcinoma development and survival. *PLoS ONE*, 2008.  
GEO: [GSE10072](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10072)

**Platform:**  
Affymetrix Human Genome U133A Array (GPL96)

-----

## License

MIT — see <LICENSE> for details.
