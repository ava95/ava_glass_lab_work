import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import my_netzoo
from scipy.stats import pearsonr
import pickle
from datetime import datetime



# LIONESS -> PANDA
def l2p(expression_df, motif, method='pearson'):
    """_summary_

    Args:
        expression_df (DataFrame): DataFrame of miRNA x sample expression data
        motif (DataFrame): DataFrame of known edge x TF/gene/weight
        tf_names (_type_): _description_
        gene_names (_type_): _description_
        method (str, optional): _description_. Defaults to 'pearson'.
    """
    
    # Establish tf/gene names
    motifs = pd.pivot_table(motif, values="Weight", index="TF", columns="Gene", fill_value=0)
    tf_names = list(motifs.index)
    gene_names = list(motifs.columns)
    

    # Generate correlation matricies
    l2p_corr = my_netzoo.my_lioness(expression_df, method)

    # Run PANDA on each correlation matrix
    l2p_networks = {}
    for key in l2p_corr:
        l2p_network = my_netzoo.my_panda(expression_df=None, motif_df=motif, correlation_matrix=l2p_corr[key])
        l2p_networks[key] = l2p_network
        
    # Assign names to genes/TFs
    l2p_networks_df = {
        key: pd.DataFrame(mat.T, index=tf_names, columns=gene_names) 
        for key, mat in l2p_networks.items()
    }
    
    return(l2p_networks_df)

# PANDA -> LIONESS
def p2l(expression_df, motif, method="panda"):
    return(my_netzoo.my_lioness(expression_df, method, motif))


# Turn raw single sample network dict into long dataframe (row = single edge in one single sample)
def make_long(network_dict, tf_names, gene_names):
    long_format = []
    for sample, arr in network_dict.items():
        df = pd.DataFrame(arr, index=tf_names, columns=gene_names)
        melted = df.reset_index().melt(id_vars='index', var_name='Gene', value_name='Score')
        melted = melted.rename(columns={'index': 'TF'})
        melted['sample'] = sample
        long_format.append(melted)
        
    all_long = pd.concat(long_format, ignore_index=True)
    
    return(all_long)


# Turn raw single sample network dict into long dataframe (row = sample, col = edge, value = weight)
def make_wide(network_dict):
    wide_rows = []
    i = 1
    for sample, df in network_dict.items():
        print(sample, i)
        # df: index = TF, columns = Gene
        row = {'sample': sample}
        for tf in df.index:
            for gene in df.columns:
                row[f"{tf}_{gene}"] = df.loc[tf, gene]
        wide_rows.append(row)
        i+=1
    network_df = pd.DataFrame(wide_rows)
    
    new_cols = []
    # Rename columns to TF_Gene
    for c in network_df.columns:
        if "_" in c:
            a, b = c.split("_")
            new_cols.append(f"{b}_{a}")
        else:
            new_cols.append(c)
    network_df.columns = new_cols
    
    return(network_df)
    
   
def main(): 
    date = datetime.today().strftime('%Y-%m-%d-%H:%M:%S')
    # Load in miRNA expression data, sample phenotype data, and motif file
    final_expr_full = pd.read_csv("../data/final_expr.txt", sep ="\t")
    pheno_df = pd.read_csv("../data/pheno_df.csv", index_col=0)
    motif = pd.read_csv("../data/motif.txt", sep = "\t")

    # Generate single-sample L2P and P2L networks
    l2p_networks = l2p(final_expr_full, motif)
    p2l_networks = p2l(final_expr_full, motif, "panda")

    with open("../data/l2p_raw.pkl", "wb") as f:
        pickle.dump(l2p_networks, f)
        
    with open("../data/p2l_raw.pkl", "wb") as f:
        pickle.dump(p2l_networks, f)

    # Generate wide df of networks (sample x edge)
    l2p_wide = make_wide(l2p_networks)
    l2p_wide.to_feather("../data/" + date + "_l2p_networks_wide.feather")

    # The below script had to be run separately in code/generate_p2l_merged.py due to memory issues.
    p2l_wide = make_wide(p2l_networks)
    p2l_wide.to_feather("../data/" + date + "_p2l_networks_wide.feather")
    
    # The below script had to be run separately in code/generate_merged_networks.py due to memory issues.
    l2p_merged = l2p_wide.merge(pheno_df, left_on='sample', right_on='title', how='left')
    l2p_merged.to_feather("../data/" + date + "_l2p_merged.feather")

    p2l_merged = p2l_wide.merge(pheno_df, left_on='sample', right_on='title', how='left')
    p2l_merged.to_feather("../data/" + date + "_p2l_merged.feather")

main()







# # Save Pearson correlation data
# l2p_corr = my_netzoo.my_lioness(final_expr_full, "pearson")
# long_format = []

# for sample, df in l2p_corr.items():
#     melted = df.reset_index().melt(id_vars=df.index.name or 'index', var_name='gene2', value_name='correlation')
#     melted = melted.rename(columns={df.index.name or 'index': 'gene1'})
#     melted['sample'] = sample
#     long_format.append(melted)

# all_long = pd.concat(long_format, ignore_index=True)
# all_long.to_feather("../data/l2p_corr.feather")


# # Save L2P networks
# long_format = []
# for sample, arr in l2p_networks.items():
#     df = pd.DataFrame(arr, index=tf_names, columns=gene_names)
#     melted = df.reset_index().melt(id_vars='index', var_name='Gene', value_name='Score')
#     melted = melted.rename(columns={'index': 'tf'})
#     melted['sample'] = sample
#     long_format.append(melted)

# all_long = pd.concat(long_format, ignore_index=True)
# all_long.to_feather('l2p_networks_long.feather')


# # Save P2L networks
# long_format = []
# for sample, arr in p2l_networks.items():
#     df = pd.DataFrame(arr, index=tf_names, columns=gene_names)
#     melted = df.reset_index().melt(id_vars='index', var_name='Gene', value_name='Score')
#     melted = melted.rename(columns={'index': 'TF'})
#     melted['sample'] = sample
#     long_format.append(melted)

# all_long = pd.concat(long_format, ignore_index=True)
# all_long.to_feather('p2l_networks_long.feather')