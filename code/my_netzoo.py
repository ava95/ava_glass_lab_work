import numpy as np
import pandas as pd
from scipy.stats import zscore
import math
from scipy import stats

def t_function(x, y=None):
    if y is None:
        a_matrix = np.dot(x, x.T)
        s = np.square(x).sum(axis=1)
        denom = s + s.reshape(-1, 1) - np.abs(a_matrix)
        denom = np.clip(denom, 1e-12, None) 
        a_matrix /= np.sqrt(denom)
    else:
        a_matrix = np.dot(x, y)
        denom = np.square(y).sum(axis=0) + np.square(x).sum(axis=1).reshape(-1, 1) - np.abs(a_matrix)
        denom = np.clip(denom, 1e-12, None) 
        a_matrix /= np.sqrt(denom)
    return a_matrix


def update_diagonal(diagonal_matrix, num, alpha, step):
    m = diagonal_matrix.copy()
    mask = ~np.eye(m.shape[0], dtype=bool)
    diagonal_std = np.std(m[mask].reshape(m.shape[0], -1), axis=1, ddof=0)
    diagonal_fill = diagonal_std * num * np.exp(2 * alpha * step)
    np.fill_diagonal(m, diagonal_fill)
    return m


def normalize_network(x):
    norm_col = zscore(x, ddof=0, axis=0)
    if x.shape[0] == x.shape[1] and np.allclose(x, x.T):
        norm_row = norm_col.T
    else:
        norm_row = zscore(x, ddof=0, axis=1)
    normalized_matrix = (norm_col + norm_row) / np.sqrt(2)
    
    normalized_matrix = (norm_col + norm_row) / math.sqrt(2)
    norm_total = (x - np.mean(x)) / np.std(x, ddof=1) 
    nan_col = np.isnan(norm_col)
    nan_row = np.isnan(norm_row)
    normalized_matrix[nan_col] = (norm_row[nan_col] + norm_total[nan_col]) / math.sqrt(
        2
    )
    normalized_matrix[nan_row] = (norm_col[nan_row] + norm_total[nan_row]) / math.sqrt(
        2
    )
    normalized_matrix[nan_col & nan_row] = (
        2 * norm_total[nan_col & nan_row] / math.sqrt(2)
    )
    
    return normalized_matrix

def compute_panda(correlation_matrix, ppi_matrix, motif_matrix, alpha=0.1, threshold=0.001):
    num_tfs, num_genes = motif_matrix.shape
    step = 0
    hamming = 1.0
    motif = motif_matrix.copy()
    ppi = ppi_matrix.copy()
    corr = correlation_matrix.copy()
    while hamming > threshold:
        W = 0.5 * (t_function(ppi, motif) + t_function(motif, corr))
        hamming = np.abs(motif - W).mean()
        print(f"step: {step}, hamming: {hamming:.6f}")
        motif = (1-alpha)*motif + alpha*W
        ppi_next = t_function(motif)
        ppi_next = update_diagonal(ppi_next, num_tfs, alpha, step)
        ppi = (1-alpha)*ppi + alpha*ppi_next
        corr_next = t_function(motif.T)
        corr_next = update_diagonal(corr_next, num_genes, alpha, step)
        corr = (1-alpha)*corr + alpha*corr_next
        step += 1
        if not np.isfinite(hamming): 
            print("Warning: hamming is not finite. Exiting.")
            break
    return motif

def my_panda(expression_df, motif_df, ppi_df=None, correlation_matrix=None, alpha=0.1, threshold=0.001):
    # Prepare correlation matrix (PANDA: Pearson by default, unless user provides one)
    if correlation_matrix is None:
        expr = expression_df.values if isinstance(expression_df, pd.DataFrame) else np.array(expression_df)
        correlation_matrix = np.corrcoef(expr)
    elif isinstance(correlation_matrix, pd.DataFrame):
        correlation_matrix = correlation_matrix.values
    else:
        correlation_matrix = np.array(correlation_matrix)
    correlation_matrix = normalize_network(correlation_matrix)
    
    if motif_df is None:
        expr = expression_df.values if isinstance(expression_df, pd.DataFrame) else np.array(expression_df)
        return np.corrcoef(expr)
    
    # Motif matrix: original PANDA is "TF, Gene, Weight" as prior 
    motifs = pd.pivot_table(motif_df, values="Weight", index="TF", columns="Gene", fill_value=0)
    motif_matrix = normalize_network(motifs.values)


    # PPI matrix
    if ppi_df is not None:
        ppi_table = pd.pivot_table(ppi_df, values="Weight", index="TF1", columns="TF2", fill_value=0)
        tf_names = motifs.index
        ppi_table = ppi_table.reindex(index=tf_names, columns=tf_names, fill_value=0)
        ppi_matrix = normalize_network(ppi_table.values)
    else:
        ppi_matrix = np.identity(motif_matrix.shape[0])

    # Call main loop
    panda_network = compute_panda(correlation_matrix, ppi_matrix, motif_matrix, alpha=alpha, threshold=threshold)
    return panda_network  # [n_TF x n_Gene]






def my_lioness(expression_df, method, motif_df=None, ppi_df=None, alpha=0.1):
    """
    Returns:
      - In PANDA mode: dict mapping sample names to DataFrames (TF x Gene)
      - In Pearson mode: dict mapping sample names to DataFrames (Gene x Gene)
    """
    expr = expression_df.values if isinstance(expression_df, pd.DataFrame) else np.array(expression_df)
    n_genes, n_samples = expr.shape
    N = n_samples
    
    # Get gene and sample names for later
    gene_names = list(expression_df.index)
    sample_names = list(expression_df.columns)
    
    if method == "panda" and motif_df is not None:
        # Motif matrix for TF/gene ordering
        motifs = pd.pivot_table(motif_df, values="Weight", index="TF", columns="Gene", fill_value=0)
        tf_names = list(motifs.index)
        gene_order = list(motifs.columns)
        fullnet = my_panda(expression_df, motif_df, ppi_df, alpha=alpha)
        lioness_networks = {}
        for q, sample in enumerate(sample_names):
            expr_minus_q = np.delete(expr, q, axis=1)
            expr_minus_q_df = pd.DataFrame(expr_minus_q, index=gene_names, columns=[c for i, c in enumerate(sample_names) if i != q])
            panda_minus_q = my_panda(expr_minus_q_df, motif_df, ppi_df, alpha=alpha)
            lioness_net = N * fullnet - (N - 1) * panda_minus_q
            lioness_networks[sample] = pd.DataFrame(lioness_net, index=tf_names, columns=gene_order)
        return lioness_networks

    else: 
        print("Running in Pearson Mode")
        fullnet = np.corrcoef(expr)
        lioness_networks = {}
        for q, sample in enumerate(sample_names):
            expr_minus_q = np.delete(expr, q, axis=1)
            net_minus_q = np.corrcoef(expr_minus_q)
            lioness_net = N * fullnet - (N - 1) * net_minus_q
            lioness_networks[sample] = pd.DataFrame(lioness_net, index=gene_names, columns=gene_names)
        return lioness_networks
    

def bare_ols_with_pvalue(X, y):
    """
    X: n x p matrix (with intercept column included if desired)
    y: n-vector
    Returns: beta, SE(beta), t-stat, p-value
    """
    n, p = X.shape
    beta = np.linalg.inv(X.T @ X) @ X.T @ y
    y_hat = X @ beta
    resid = y - y_hat
    df = n - p
    s2 = (resid @ resid) / df
    var_beta = s2 * np.linalg.inv(X.T @ X)
    se_beta = np.sqrt(np.diag(var_beta))
    t_stats = beta / se_beta
    p_values = 2 * stats.t.sf(np.abs(t_stats), df)
    return beta, se_beta, t_stats, p_values