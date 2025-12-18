# Physics-Based Therapeutic Strategy for Triple-Negative Breast Cancer (TNBC) Resistance

## Abstract

Acquired drug resistance in Triple-Negative Breast Cancer (TNBC) remains a fatal clinical hurdle, often attributed to stochastic genetic mutations. However, recent evidence suggests resistance may emerge as a non-genetic phase transition into a high-entropy attractor state. In this study, we utilized single-cell RNA sequencing (scRNA-seq) data from the GSE113197 cohort to map the “Energy Landscape” of TNBC progression.

Using Topological Data Analysis (TDA) and Diffusion Pseudotime (DPT) inference, we identified a critical bifurcation point at Pseudotime 0.7, marking a sudden, first-order phase transition from a proliferative epithelial state to a dormant, drug-tolerant mesenchymal state (Leiden Cluster 9). This transition was not driven by FOXP1 alone, but was preceded by the upregulation of **COL17A1**, a hemidesmosome anchor protein, and **MT2A**, a stress-response metallothionein. Survival analysis of the TCGA-BRCA cohort (n = 1,082) confirmed that high **COL17A1** expression is significantly associated with poor overall survival (*p* = 0.040), validating it as a clinically relevant target.

Finally, we performed an *in silico* gene knockout simulation targeting the COL17A1/MT2A regulatory module. The simulation resulted in a 100% collapse of the resistant phenotype score (Shift: 1.5 → 0.0), effectively destabilizing the resistant attractor. These findings suggest that resistance in TNBC is a mechanically anchored thermodynamic state that can be reversed by targeting the structural proteins maintaining cellular dormancy, offering a novel “physics-based” therapeutic strategy.

---

## Table of contents

1. [Introduction / Background](#introduction)
2. [Data](#data)
3. [Requirements](#requirements)
4. [Quick start (copy & paste)](#quick-start)
5. [Pipeline (detailed)](#pipeline)
6. [Core functions & notebooks](#core)
7. [Reproducing the *in silico* knockout (virtual cure)](#virtual-cure)
8. [Figures and outputs](#figures)
9. [Limitations and interpretation](#limitations)
10. [Reproducibility notes](#reproducibility)
11. [How to cite / contact / contribute](#cite)
12. [License](#license)

---

## 1. Introduction / Background <a name="introduction"></a>

Triple-Negative Breast Cancer (TNBC) lacks ER/PR/HER2 expression and is clinically aggressive. Conventional explanations for therapy failure emphasize genetic evolution; this work explores a complementary hypothesis: that drug resistance can arise as a *non-genetic thermodynamic phase transition* in cellular transcriptional state space. By computing entropy and energy-landscape proxies on single-cell transcriptomes and inferring trajectories (DPT), we identify a discontinuous bifurcation into a resistant attractor anchored by structural adhesion machinery (notably COL17A1).

---

## 2. Data <a name="data"></a>

Two primary datasets were used:

* **scRNA-seq (GSE113197)** — single-cell transcriptomes used to build the cellular energy landscape, compute entropy, cluster (Leiden), and infer trajectory (DPT).
* **TCGA-BRCA (bulk RNA + clinical)** — used for independent clinical validation (survival analysis) of candidate driver genes (COL17A1).

**Data bundle (optional):** The analysis referenced a Google Drive bundle with preprocessed files:
`https://drive.google.com/drive/folders/1l0CpWbRYPziaagnd19w1uboj-lXtnKZR`

> If you do not have the Drive bundle, the notebook contains helper cells to download raw GEO files and TCGA expression/clinical data as appropriate.

---

## 3. Requirements <a name="requirements"></a>

This analysis is implemented in Python. The minimal set of libraries required:

* Python 3.8+
* `scanpy`
* `anndata`
* `numpy`
* `pandas`
* `scipy`
* `matplotlib`
* `seaborn` (optional, for nicer plots)
* `umap-learn` (if not provided through `scanpy`)
* `python-igraph` & `leidenalg` (for clustering)
* `mygene` (optional, for ID mapping)
* `lifelines` (for Kaplan–Meier / log-rank tests)
* `tqdm` (progress bars)

A `requirements.txt` is included in the repository. For exact reproducibility, pin versions (see Reproducibility).

---

## 4. Quick start (copy & paste) <a name="quick-start"></a>

### Clone repository

```bash
git clone https://github.com/sayandeep520/Thermodynamic-Resistance-TNBC.git
cd Thermodynamic-Resistance-TNBC
```

### Create virtual environment and install dependencies

```bash
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

> If you want pinned versions, first create `requirements-pinned.txt` by running `pip freeze` in the environment used to generate the repository (recommendation in Reproducibility).

### Launch Jupyter and open the master notebook

```bash
jupyter lab
# then open Triple_Negative_Breast_Cancer_Master_File.ipynb
```

### Run in Google Colab

1. Upload the notebook to your Google Drive or open the notebook from the GitHub link in Colab.
2. Mount Drive (notebook provides a cell) and set `DATA_DIR`.
3. Install dependencies in a Colab cell:

```python
!pip install -r /content/drive/MyDrive/path/to/requirements.txt
```

4. Run the notebook cells sequentially (or “Runtime → Run all”).

---

## 5. Pipeline (detailed) <a name="pipeline"></a>

The master notebook (`Triple_Negative_Breast_Cancer_Master_File.ipynb`) implements the following canonical pipeline:

1. **Data loading / AnnData assembly**

   * Import raw counts or load a preprocessed `breast_cancer_atlas.h5ad`.

2. **Quality control (QC)**

   * Filter cells with low counts, high mitochondrial fraction, and genes expressed in very few cells.

3. **Normalization & log transform**

   * `scanpy.pp.normalize_total(...); scanpy.pp.log1p(...)`.

4. **Highly variable genes (HVG) & scaling**

   * Select HVGs and scale the matrix (`scanpy.pp.highly_variable_genes`, `scanpy.pp.scale`).

5. **Dimensionality reduction**

   * PCA → neighbor graph (`scanpy.pp.neighbors`) → UMAP (`scanpy.tl.umap`).

6. **Clustering**

   * Leiden clustering (`scanpy.tl.leiden`), annotate clusters; resistance identified as Leiden cluster 9.

7. **Entropy / physics metrics**

   * Compute Shannon entropy per cell and an energy proxy (functions in `physics_metrics.py`).

8. **Trajectory inference**

   * Diffusion pseudotime (DPT), identify a bifurcation at pseudotime ≈ 0.7.

9. **Marker detection & resistance network**

   * `scanpy.tl.rank_genes_groups` for markers; compute pairwise correlations to build resistance network.

10. **TCGA-BRCA survival validation**

    * Map COL17A1 expression to patient overall survival and compute Kaplan–Meier curves & log-rank p-value.

11. **In silico knockout (virtual cure)**

    * Zero or scale down expression of target genes (COL17A1, MT2A, partners) in resistant cells, recompute resistance score and project results.

---

## 6. Core functions & notebooks <a name="core"></a>

* `Triple_Negative_Breast_Cancer_Master_File.ipynb` — canonical pipeline and figure generation (run top-to-bottom).
* `physics_metrics.py` — compute Shannon entropy, resistance identity score, and related metrics. Inspect docstrings for function signatures.
* `visualization.py` — plotting utilities (UMAP overlays, pseudotime profiles, KM-plot wrappers).
* `Cluster9_Resistance_Signature.csv` — signature gene list used to compute resistance scores.
* `Resistance_Network_Correlations.csv` — gene–gene correlations for the resistance module.

**Usage examples (copy & paste):**

Compute entropy and resistance score:

```python
from physics_metrics import compute_shannon_entropy, resistance_identity_score
import pandas as pd

# compute entropy
adata.obs['entropy'] = compute_shannon_entropy(adata, layer='X')

# load signature and compute score
sig_genes = pd.read_csv('Cluster9_Resistance_Signature.csv')['gene'].tolist()
adata.obs['resistance_score'] = resistance_identity_score(adata, sig_genes)
```

Plot COL17A1 pseudotime profile:

```python
from visualization import plot_pseudotime_profile
plot_pseudotime_profile(adata, gene='COL17A1', pseudotime_key='dpt_pseudotime', save='phase_transition_COL17A1.png')
```

---

## 7. Reproducing the *in silico* knockout (virtual cure) <a name="virtual-cure"></a>

This procedure is implemented in the notebook and illustrated here for reproducibility.

**Steps (copy & paste):**

```python
# assume adata is your AnnData and sig_genes is signature list
import numpy as np
adata_treated = adata.copy()
target_genes = ['COL17A1', 'MT2A']  # extend this list with correlated partners from Resistance_Network_Correlations.csv

# Hard knockout (set expression to zero for target genes)
for g in target_genes:
    if g in adata_treated.var_names:
        adata_treated[:, g].X = np.zeros(adata_treated[:, g].X.shape)

# recompute resistance score on treated copy
adata_treated.obs['resistance_score_treated'] = resistance_identity_score(adata_treated, sig_genes)

# compare distributions (example: Kolmogorov–Smirnov test)
from scipy.stats import ks_2samp
ks_stat, ks_p = ks_2samp(adata.obs['resistance_score'], adata_treated.obs['resistance_score_treated'])
print(f'KS p-value: {ks_p:.4g}')
```

**Expected outcome in the repository analysis:** Near-complete collapse of the resistance identity score (example reported shift: 1.5 → 0.0). The magnitude depends on score definition and normalization; verify exact scoring implementation in `physics_metrics.py`.

---

## 8. Figures and outputs <a name="figures"></a>

The notebook generates and saves the primary figures included in the repo:

* `trajectory_pseudotime.png` — UMAP colored by DPT pseudotime.
* `energy_landscape_entropy.png` — entropy mapped on UMAP (energy landscape view).
* `phase_transition_COL17A1.png` — expression profile of COL17A1 along pseudotime (phase transition).
* `phase_transition_FOXP1.png` — FOXP1 pseudotime profile.
* `survival_km_plot_COL17A1.png` — Kaplan–Meier curve for COL17A1 (TCGA-BRCA).
* `virtual_cure_histogram.png` — resistance score distribution before/after virtual cure.
* `virtual_cure_umap_collapse.png` — UMAP showing collapse of resistant cluster after in silico perturbation.

To regenerate figures, run the corresponding notebook cells or call the plotting helpers in `visualization.py`.

---

## 9. Limitations and interpretation <a name="limitations"></a>

* **Model dependence:** Trajectory and pseudotime inference depend strongly on preprocessing choices (HVG selection, neighbors, root selection). Results must be interpreted in that context.
* **In silico perturbations are idealized:** Setting transcript counts to zero is not equivalent to biological inhibition and ignores post-transcriptional regulation and protein stability.
* **Correlation vs causation:** The resistance network is correlation-based; experimental perturbation (CRISPR, RNAi, small molecules) is required to establish causal roles.
* **Bulk validation caveat:** TCGA analyses are bulk; they validate association with clinical outcome but cannot resolve single-cell heterogeneity directly.

---

## 10. Reproducibility notes <a name="reproducibility"></a>

* Fix random seeds at the start of the notebook:

```python
import numpy as np, random, scanpy as sc
np.random.seed(0)
random.seed(0)
sc.settings.seed = 0
```

* Provide `random_state` to UMAP / Leiden where supported for deterministic runs.
* Save intermediate processed `AnnData` (e.g., `breast_cancer_atlas.h5ad`) after QC & normalization to avoid rerunning stochastic steps.
* Pin package versions: run `pip freeze > requirements-pinned.txt` in the environment you use to generate the repo results. Commit `requirements-pinned.txt` to the repo for exact reproducibility.
* Document critical parameter values (HVG count, PCA components, Leiden resolution, DPT root) in a `params.yml` or in the notebook header.

---

## 11. How to cite / contact / contribute <a name="cite"></a>

**Suggested citation:**

> Bera SD. *Thermodynamic-Resistance-TNBC*. GitHub repository, 2025. [https://github.com/sayandeep520/Thermodynamic-Resistance-TNBC](https://github.com/sayandeep520/Thermodynamic-Resistance-TNBC)

**Contact & issues:** Use the GitHub Issues page to report problems, request data, or propose changes.

**Contributing:** Fork → branch (`feature/...`) → add tests/notebook → open a Pull Request with a clear description and reproducible example.

---

## 12. License <a name="license"></a>

This repository is distributed under the **MIT License**. See `LICENSE` for full text.

---

## Final notes

* The master notebook is the canonical source of analysis logic and parameters. When in doubt, follow the notebook cells in order.
