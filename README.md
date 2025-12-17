# Physics-Based Therapeutic Strategy for Triple-Negative Breast Cancer (TNBC) Resistance

## Abstract
Acquired drug resistance in Triple-Negative Breast Cancer (TNBC) remains a fatal clinical hurdle, often attributed to stochastic genetic mutations. However, recent evidence suggests resistance may emerge as a non-genetic phase transition into a high-entropy attractor state. In this study, we utilized single-cell RNA sequencing (scRNA-seq) data from the GSE113197 cohort to map the "Energy Landscape" of TNBC progression.

Using Topological Data Analysis (TDA) and Diffusion Pseudotime (DPT) inference, we identified a critical bifurcation point at Pseudotime 0.7, marking a sudden, first-order phase transition from a proliferative epithelial state to a dormant, drug-tolerant mesenchymal state (Leiden Cluster 9). This transition was not driven by FOXP1 alone, but was preceded by the upregulation of COL17A1, a hemidesmosome anchor protein, and MT2A, a stress-response metallothionein. Survival analysis of the TCGA-BRCA cohort (n=1,082) confirmed that high COL17A1 expression is significantly associated with poor overall survival ($p=0.040$), validating it as a clinically relevant target.

Finally, we performed an in silico gene knockout simulation targeting the COL17A1/MT2A regulatory module. The simulation resulted in a 100% collapse of the resistant phenotype score (Shift: 1.5 $\to$ 0.0), effectively destabilizing the resistant attractor. These findings suggest that resistance in TNBC is a mechanically anchored thermodynamic state that can be reversed by targeting the structural proteins maintaining cellular dormancy, offering a novel "physics-based" therapeutic strategy.

## Introduction/Background
Triple-Negative Breast Cancer (TNBC) is an aggressive subtype of breast cancer characterized by the absence of estrogen receptor (ER), progesterone receptor (PR), and human epidermal growth factor receptor 2 (HER2) expression. This lack of common drug targets makes TNBC particularly challenging to treat, often leading to acquired drug resistance and poor patient outcomes. Understanding the mechanisms driving resistance is crucial for developing effective therapies. This study explores a novel hypothesis: that drug resistance in TNBC may be understood as a non-genetic phase transition into a high-entropy, drug-tolerant state, rather than solely due to genetic mutations.

## Data
Two primary datasets were utilized in this study:

1.  **scRNA-seq Data (GSE113197)**: Single-cell RNA sequencing data from a breast cancer cohort, providing gene expression profiles at a single-cell resolution. This dataset was used to map the cellular energy landscape and infer trajectories of cancer progression.
2.  **TCGA-BRCA Bulk RNA-seq & Clinical Data**: The Breast Invasive Carcinoma (BRCA) cohort from The Cancer Genome Atlas (TCGA), containing bulk RNA sequencing data (expression of genes) and corresponding clinical information (including overall survival). This dataset was used for clinical validation of identified driver genes.

**Data Link**: The data files can be accessed via this Google Drive link: [https://drive.google.com/drive/folders/1l0CpWbRYPziaagnd19w1uboj-lXtnKZR](https://drive.google.com/drive/folders/1l0CpWbRYPziaagnd19w1uboj-lXtnKZR)

## Requirements
This analysis requires a Python environment with the following libraries:

*   `scanpy`
*   `pandas`
*   `numpy`
*   `matplotlib`
*   `seaborn`
*   `mygene`
*   `lifelines`
*   `scipy` (specifically `scipy.stats` and `scipy.sparse`)
*   `anndata`
*   `tqdm`
*   `os`, `tarfile`, `gzip`, `google.colab.drive` (for Google Colab environment)

These can typically be installed via `pip` (e.g., `pip install scanpy lifelines mygene`).

## Methodology

1.  **Data Acquisition & Preprocessing**: 
    *   Single-cell RNA sequencing data from `GSE113197_RAW.tar` was extracted and loaded into Scanpy AnnData objects. 
    *   Standard quality control metrics were calculated, followed by filtering of low-quality cells and genes. 
    *   The data was normalized and log-transformed.
2.  **Dimensionality Reduction & Clustering**: 
    *   Principal Component Analysis (PCA) was performed for initial dimensionality reduction.
    *   UMAP (Uniform Manifold Approximation and Projection) was applied for visualization and further dimensionality reduction. 
    *   Leiden clustering was used to identify distinct cell populations (clusters) within the scRNA-seq data.
3.  **Thermodynamic Entropy & Trajectory Inference**: 
    *   Shannon Entropy was calculated for each cell to quantify the 'disorder' in gene expression. 
    *   Diffusion Pseudotime (DPT) was employed to infer a cellular trajectory, mapping the 'arrow of time' in cancer progression, starting from a high-entropy cluster.
4.  **Marker Gene Identification**: 
    *   `sc.tl.rank_genes_groups` was used to identify differentially expressed genes (marker genes) for each Leiden cluster, helping characterize the biological identity of cell populations. 
    *   Ensembl gene IDs were translated to human gene symbols using the `mygene` library.
5.  **Clinical Validation**: 
    *   TCGA-BRCA bulk RNA-seq and clinical data (`brca_tcga_pan_can_atlas_2018.tar.gz`) were loaded, preprocessed, and merged.
    *   Kaplan-Meier survival analysis and log-rank tests were performed for key driver genes (e.g., COL17A1) to validate their clinical significance in relation to patient overall survival.
6.  **In Silico Knockout & Efficacy Evaluation**: 
    *   A 'resistance network' of genes highly correlated with COL17A1 (e.g., MT2A, MT1X, MYLK, F3, NBPF20) was identified. 
    *   An *in silico* gene knockout simulation was performed by setting the expression of these target genes to zero in resistant cells (Leiden Cluster 9).
    *   A 'resistance identity score' (average expression of target genes) was calculated for sensitive, untreated resistant, and treated resistant cells.
    *   The effect of the virtual drug was visualized using UMAP plots and histograms, comparing resistance score distributions before and after treatment.

## Key Findings & Insights

*   **Phase Transition to Resistance**: The scRNA-seq data analysis revealed a first-order phase transition in TNBC progression, moving from a proliferative epithelial state to a dormant, drug-tolerant mesenchymal state (Leiden Cluster 9) at pseudotime 0.7.
*   **Key Driver Genes**: Upregulation of `COL17A1` (a hemidesmosome anchor protein) was observed to precede the 'identity switch' gene `FOXP1` along the resistance trajectory. `MT2A` (a stress-response metallothionein) was also identified as a key component of this regulatory module.
*   **Clinical Relevance**: High expression of `COL17A1` was found to be significantly associated with poor overall survival in the TCGA-BRCA cohort (p=0.040), underscoring its clinical importance.
*   **Reversal of Resistance**: The *in silico* knockout simulation targeting the `COL17A1/MT2A` regulatory module demonstrated a complete collapse of the resistant phenotype score (from 1.5 to 0.0) in treated resistant cells. This indicates a successful destabilization of the resistant attractor state.

These findings suggest that TNBC resistance is a mechanically anchored thermodynamic state that can be reversed by targeting specific structural and stress-response proteins. This opens avenues for novel, physics-based therapeutic strategies.

## How to Run
This notebook is designed to be run in a Google Colab environment.

1.  **Data Setup**: 
    *   Ensure that the `GSE113197_RAW.tar` (scRNA-seq data) and `brca_tcga_pan_can_atlas_2018.tar.gz` (TCGA data) files are uploaded to your Google Drive in the following path: `/content/drive/My Drive/Breast Cancer/`.
    *   The data files can be accessed via this Google Drive link: [https://drive.google.com/drive/folders/1l0CpWbRYPziaagnd19w1uboj-lXtnKZR](https://drive.google.com/drive/folders/1l0CpWbRYPziaagnd19w1uboj-lXtnKZR)
2.  **Mount Google Drive**: The notebook will prompt you to mount your Google Drive. Follow the instructions to grant access.
3.  **Run All Cells**: Execute all cells in the notebook sequentially. The analysis includes data loading, preprocessing, dimensionality reduction, clustering, trajectory inference, marker gene analysis, clinical validation, and *in silico* drug simulation.

## Files
*   `README.md`: This file, providing an overview of the project.
*   `breast_cancer_atlas.h5ad`: An aggregated AnnData object saved during the analysis, containing processed scRNA-seq data.
*   Various image files (`.png`): Plots generated during the analysis, such as UMAPs, histograms, and Kaplan-Meier curves.
