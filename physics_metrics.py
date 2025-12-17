
import numpy as np
import scipy.sparse
import scipy.stats

def calculate_shannon_entropy(adata):
    # Get the raw expression matrix (cells x genes)
    expression_matrix = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X

    # Normalize so probabilities sum to 1 for each cell (p_i)
    row_sums = expression_matrix.sum(axis=1)
    row_sums[row_sums == 0] = 1e-9 # Avoid division by zero
    probs = expression_matrix / row_sums[:, np.newaxis]

    # Calculate S = -sum(p_i * log(p_i))
    entropy = scipy.stats.entropy(probs.T)

    return entropy
