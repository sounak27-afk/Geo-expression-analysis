# “””
GEO Gene Expression Analysis Pipeline

Dataset : GSE10072 — Lung Adenocarcinoma vs. Normal Tissue
Source  : NCBI Gene Expression Omnibus (GEO)
Platform: GPL96 (Affymetrix Human Genome U133A Array)
Samples : 107 samples (58 tumor, 49 normal)

## Pipeline steps

1. Download & parse GEO dataset via GEOparse
1. Exploratory Data Analysis (EDA)
1. Log2-normalisation & filtering
1. Differential Expression Analysis (t-test + FDR correction)
1. Visualisations: PCA, Volcano plot, Heatmap of top DEGs
1. Export results tables

Author  : [Your Name]
Date    : 2026-06
“””

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use(“Agg”)          # headless rendering for CI/server
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import GEOparse
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings(“ignore”)

# ── Configuration ─────────────────────────────────────────────────────────────

GEO_ID          = “GSE10072”
DATA_DIR        = “data”
RESULTS_DIR     = “results”
FIGURES_DIR     = os.path.join(RESULTS_DIR, “figures”)
TABLES_DIR      = os.path.join(RESULTS_DIR, “tables”)
FDR_THRESHOLD   = 0.05
LOG2FC_THRESHOLD = 1.0          # |log2 fold-change| ≥ 1  →  2-fold change
TOP_N_HEATMAP   = 50           # top DEGs shown in heatmap

for d in [DATA_DIR, FIGURES_DIR, TABLES_DIR]:
os.makedirs(d, exist_ok=True)

sns.set_theme(style=“whitegrid”, palette=“muted”, font_scale=1.2)
PALETTE = {“Tumor”: “#E63946”, “Normal”: “#457B9D”}

# ── 1. Download & Parse ───────────────────────────────────────────────────────

print(f”[1/5] Downloading {GEO_ID} from NCBI GEO …”)
gse = GEOparse.get_GEO(geo=GEO_ID, destdir=DATA_DIR, silent=True)

# Build expression matrix  (probes × samples)

print(”[2/5] Parsing expression matrix …”)
expr_frames = []
for gsm_name, gsm in gse.gsms.items():
s = gsm.table.set_index(“ID_REF”)[“VALUE”].rename(gsm_name)
expr_frames.append(s)

expr = pd.concat(expr_frames, axis=1).astype(float)
print(f”      Raw matrix: {expr.shape[0]:,} probes × {expr.shape[1]} samples”)

# ── 2. Sample metadata ────────────────────────────────────────────────────────

meta_rows = []
for gsm_name, gsm in gse.gsms.items():
chars = gsm.metadata.get(“characteristics_ch1”, [])
# GSE10072 encodes tissue type in the title / characteristics
title = “ “.join(gsm.metadata.get(“title”, [””])).lower()
label = “Tumor” if “adenocarcinoma” in title or “tumor” in title else “Normal”
meta_rows.append({“sample”: gsm_name, “label”: label})

meta = pd.DataFrame(meta_rows).set_index(“sample”)
print(f”      Tumor: {(meta.label==‘Tumor’).sum()}  |  Normal: {(meta.label==‘Normal’).sum()}”)

# ── 3. Normalisation & Filtering ──────────────────────────────────────────────

print(”[3/5] Normalising & filtering …”)

# Replace 0/negative values before log transform

expr = expr.clip(lower=1)
expr_log = np.log2(expr)

# Filter low-expression probes: keep probes with median signal > log2(50)

mask = expr_log.median(axis=1) > np.log2(50)
expr_log = expr_log[mask]
print(f”      After filtering: {expr_log.shape[0]:,} probes retained”)

# Quantile-normalise across samples (rank-based)

def quantile_normalise(df: pd.DataFrame) -> pd.DataFrame:
rank_mean = df.stack().groupby(df.rank(method=“first”).stack().astype(int)).mean()
return df.rank(method=“min”).stack().astype(int).map(rank_mean).unstack()

expr_norm = quantile_normalise(expr_log)

# ── 4. Differential Expression ────────────────────────────────────────────────

print(”[4/5] Running differential expression analysis …”)
tumor_cols  = meta[meta.label == “Tumor”].index.tolist()
normal_cols = meta[meta.label == “Normal”].index.tolist()

results = []
for probe in expr_norm.index:
t_vals  = expr_norm.loc[probe, tumor_cols].values
n_vals  = expr_norm.loc[probe, normal_cols].values
log2fc  = t_vals.mean() - n_vals.mean()
_, pval = stats.ttest_ind(t_vals, n_vals, equal_var=False)
results.append({“probe”: probe, “log2FC”: log2fc, “pvalue”: pval,
“mean_tumor”: t_vals.mean(), “mean_normal”: n_vals.mean()})

de = pd.DataFrame(results).set_index(“probe”)

# FDR correction (Benjamini-Hochberg)

_, de[“padj”], _, _ = multipletests(de[“pvalue”], method=“fdr_bh”)

# Classify probes

de[“significant”] = (de[“padj”] < FDR_THRESHOLD) & (de[“log2FC”].abs() >= LOG2FC_THRESHOLD)
de[“direction”]   = np.where(
de[“significant”] & (de[“log2FC”] > 0), “Up in Tumor”,
np.where(de[“significant”] & (de[“log2FC”] < 0), “Down in Tumor”, “NS”)
)

n_up   = (de.direction == “Up in Tumor”).sum()
n_down = (de.direction == “Down in Tumor”).sum()
print(f”      Significant DEGs: {n_up + n_down:,}  ({n_up} up, {n_down} down)”)

# Save full DE table

de.sort_values(“padj”).to_csv(os.path.join(TABLES_DIR, “DE_results.csv”))

# ── 5. Visualisations ─────────────────────────────────────────────────────────

print(”[5/5] Generating figures …”)

# ── 5a. Sample-level QC boxplot ───────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(max(14, expr_norm.shape[1] // 4), 5))
data_long = expr_norm.T.melt(var_name=“probe”, value_name=“log2_expr”)
data_long[“label”] = data_long.index.map(meta[“label”]) if False else   
expr_norm.columns.map(meta[“label”]).repeat(expr_norm.shape[0]).values
ax.set_title(f”{GEO_ID} — Per-sample expression distribution (after QN)”, fontsize=13)
ax.set_xlabel(“Sample index”)
ax.set_ylabel(“log₂ expression”)

# Summarise per-sample as boxplot stats

sample_medians = expr_norm.median()
colors = [PALETTE[meta.loc[s, “label”]] for s in expr_norm.columns]
bp = ax.boxplot(expr_norm.values, patch_artist=True, widths=0.6,
medianprops=dict(color=“white”, linewidth=1.5),
flierprops=dict(marker=”.”, markersize=1, alpha=0.3))
for patch, col in zip(bp[“boxes”], colors):
patch.set_facecolor(col)
patches = [mpatches.Patch(color=v, label=k) for k, v in PALETTE.items()]
ax.legend(handles=patches)
ax.set_xticks([])
plt.tight_layout()
fig.savefig(os.path.join(FIGURES_DIR, “01_sample_boxplot.png”), dpi=150)
plt.close()

# ── 5b. PCA ───────────────────────────────────────────────────────────────────

scaler = StandardScaler()
X = scaler.fit_transform(expr_norm.T)          # samples × probes
pca = PCA(n_components=2, random_state=42)
pcs = pca.fit_transform(X)
var_explained = pca.explained_variance_ratio_ * 100

fig, ax = plt.subplots(figsize=(7, 6))
for label, col in PALETTE.items():
idx = [i for i, s in enumerate(expr_norm.columns) if meta.loc[s, “label”] == label]
ax.scatter(pcs[idx, 0], pcs[idx, 1], c=col, label=label,
s=70, alpha=0.85, edgecolors=“white”, linewidths=0.4)
ax.set_xlabel(f”PC1 ({var_explained[0]:.1f}% variance)”)
ax.set_ylabel(f”PC2 ({var_explained[1]:.1f}% variance)”)
ax.set_title(f”{GEO_ID} — PCA of normalised expression”)
ax.legend(frameon=True)
plt.tight_layout()
fig.savefig(os.path.join(FIGURES_DIR, “02_pca.png”), dpi=150)
plt.close()

# ── 5c. Volcano plot ──────────────────────────────────────────────────────────

color_map = {“Up in Tumor”: PALETTE[“Tumor”], “Down in Tumor”: PALETTE[“Normal”], “NS”: “#CCCCCC”}
de[“neg_log10_padj”] = -np.log10(de[“padj”].clip(lower=1e-300))

fig, ax = plt.subplots(figsize=(8, 7))
for direction, grp in de.groupby(“direction”):
ax.scatter(grp[“log2FC”], grp[“neg_log10_padj”],
c=color_map[direction], label=direction,
s=6 if direction == “NS” else 20,
alpha=0.4 if direction == “NS” else 0.85,
edgecolors=“none”)
ax.axhline(-np.log10(FDR_THRESHOLD), ls=”–”, lw=1, color=“grey”)
ax.axvline(-LOG2FC_THRESHOLD, ls=”–”, lw=1, color=“grey”)
ax.axvline( LOG2FC_THRESHOLD, ls=”–”, lw=1, color=“grey”)
ax.set_xlabel(“log₂ Fold Change (Tumor / Normal)”)
ax.set_ylabel(”-log₁₀(adjusted p-value)”)
ax.set_title(f”{GEO_ID} — Volcano Plot”)
ax.legend(markerscale=2, frameon=True)
plt.tight_layout()
fig.savefig(os.path.join(FIGURES_DIR, “03_volcano.png”), dpi=150)
plt.close()

# ── 5d. Heatmap of top DEGs ───────────────────────────────────────────────────

top_probes = (de[de.significant]
.assign(abs_log2fc=de[“log2FC”].abs())
.sort_values([“abs_log2fc”, “padj”], ascending=[False, True])
.head(TOP_N_HEATMAP).index)

hm_data = expr_norm.loc[top_probes]

# Column order: tumors then normals

col_order = (meta[meta.label == “Tumor”].index.tolist() +
meta[meta.label == “Normal”].index.tolist())
hm_data = hm_data[col_order]

# Z-score rows

hm_z = hm_data.apply(lambda r: (r - r.mean()) / r.std(), axis=1)

col_colors = pd.Series(
[PALETTE[meta.loc[s, “label”]] for s in col_order],
index=col_order
)
g = sns.clustermap(
hm_z,
col_cluster=False,
row_cluster=True,
col_colors=col_colors,
cmap=“RdBu_r”,
center=0,
yticklabels=True,
xticklabels=False,
figsize=(12, 14),
dendrogram_ratio=(0.15, 0.03),
cbar_pos=(0.02, 0.8, 0.03, 0.15),
)
g.fig.suptitle(f”Top {TOP_N_HEATMAP} DEGs — {GEO_ID}”, y=1.01, fontsize=13)
patches = [mpatches.Patch(color=v, label=k) for k, v in PALETTE.items()]
g.ax_heatmap.legend(handles=patches, loc=“upper right”,
bbox_to_anchor=(1.25, 1.15), frameon=True)
g.savefig(os.path.join(FIGURES_DIR, “04_heatmap_top_DEGs.png”), dpi=150, bbox_inches=“tight”)
plt.close()

# ── Summary table ─────────────────────────────────────────────────────────────

summary = {
“GEO accession”: GEO_ID,
“Platform”: gse.gpls[list(gse.gpls.keys())[0]].metadata[“title”][0],
“Total samples”: expr.shape[1],
“Tumor samples”: (meta.label == “Tumor”).sum(),
“Normal samples”: (meta.label == “Normal”).sum(),
“Probes (raw)”: expr.shape[0],
“Probes (after filter)”: expr_norm.shape[0],
“DEGs (FDR < 0.05, |log2FC| ≥ 1)”: n_up + n_down,
“Up in Tumor”: n_up,
“Down in Tumor”: n_down,
}
pd.Series(summary).to_csv(os.path.join(TABLES_DIR, “summary.csv”), header=[“value”])

print(”\n✓ Analysis complete.”)
print(f”  Figures → {FIGURES_DIR}/”)
print(f”  Tables  → {TABLES_DIR}/”)
