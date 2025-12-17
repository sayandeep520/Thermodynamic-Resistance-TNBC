
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from lifelines import KaplanMeierFitter

def plot_umap_clusters(adata, color_key='leiden', title='UMAP: Clusters'):
    """Plots UMAP colored by clusters."""
    plt.figure(figsize=(8, 8))
    sc.pl.umap(adata, color=color_key, title=title, show=False)
    plt.tight_layout()
    plt.show()

def plot_gene_expression_on_umap(adata, gene_symbol, gene_symbols_key='gene_symbol', title_prefix='UMAP: '):
    """Plots single gene expression on UMAP, translating symbol to ensembl ID if needed."""
    if gene_symbol in adata.var[gene_symbols_key].values:
        gene_id_for_plot = adata.var[adata.var[gene_symbols_key] == gene_symbol].index[0]
        plt.figure(figsize=(8, 8))
        sc.pl.umap(adata, color=gene_id_for_plot, title=f'{title_prefix}{gene_symbol} Expression', show=False)
        plt.tight_layout()
        plt.show()
    else:
        print(f"Gene {gene_symbol} not found in adata.var['{gene_symbols_key}'].")

def plot_resistance_scores(adata_untreated, adata_treated, cluster_id='9'):
    """Plots side-by-side UMAPs of resistance scores and a histogram comparison."""
    sns.set_style('whitegrid')

    plt.figure(figsize=(14, 7))

    plt.subplot(1, 2, 1)
    sc.pl.umap(adata_untreated, color='resistance_score_untreated', title='Untreated: Resistance Score', ax=plt.gca(), show=False, cmap='viridis')
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')

    plt.subplot(1, 2, 2)
    sc.pl.umap(adata_treated, color='resistance_score_treated', title='Treated: Resistance Score', ax=plt.gca(), show=False, cmap='viridis')
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')

    plt.tight_layout()
    plt.show()

    print("
--- Comparing Resistance Score Distributions in Cluster 9 ---")
    untreated_resistant_scores = adata_untreated.obs[adata_untreated.obs['leiden'] == cluster_id]['resistance_score_untreated']
    treated_resistant_scores = adata_treated.obs[adata_treated.obs['leiden'] == cluster_id]['resistance_score_treated']

    plt.figure(figsize=(10, 6))
    sns.histplot(untreated_resistant_scores, bins=30, color='blue', label='Untreated Resistant Cells', kde=True, alpha=0.6)
    sns.histplot(treated_resistant_scores, bins=30, color='red', label='Treated Resistant Cells', kde=True, alpha=0.6)
    plt.title(f'Distribution of Resistance Scores in Leiden Cluster {cluster_id} (Before vs. After Treatment)')
    plt.xlabel('Resistance Score (Average Target Gene Expression)')
    plt.ylabel('Number of Cells')
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_kaplan_meier(merged_df, gene_symbol, T_col='OS_MONTHS', E_col='OS_STATUS'):
    """Plots Kaplan-Meier survival curves for a given gene."""
    kmf = KaplanMeierFitter()
    median_expression = merged_df[gene_symbol].median()
    merged_df[f'{gene_symbol}_high_expression'] = merged_df[gene_symbol] > median_expression

    high_expression_group = merged_df[merged_df[f'{gene_symbol}_high_expression'] == True]
    low_expression_group = merged_df[merged_df[f'{gene_symbol}_high_expression'] == False]

    if not high_expression_group.empty and not low_expression_group.empty:
        kmf.fit(merged_df[T_col][low_expression_group.index], event_observed=merged_df[E_col][low_expression_group.index], label=f'Low {gene_symbol} Expression')
        ax = kmf.plot(figsize=(8, 6), show_censors=True)
        kmf.fit(merged_df[T_col][high_expression_group.index], event_observed=merged_df[E_col][high_expression_group.index], label=f'High {gene_symbol} Expression')
        kmf.plot(ax=ax, show_censors=True)
        plt.title(f'Kaplan-Meier Survival Curves for {gene_symbol} Expression')
        plt.xlabel('Time (Months)')
        plt.ylabel('Survival Probability')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
    else:
        print(f"Skipping Kaplan-Meier for {gene_symbol}: groups are empty.")

