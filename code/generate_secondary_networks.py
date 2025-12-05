import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import my_netzoo
from scipy.stats import pearsonr
import os


# Load in miRNA expression data, sample phenotype data, and motif file
final_expr_full = pd.read_csv("../data/final_expr.txt", sep ="\t")
pheno_df = pd.read_csv("../data/pheno_df.csv", index_col=0)
motif = pd.read_csv("../data/motif.txt", sep = "\t")

# Establish tf/gene names
motifs = pd.pivot_table(motif, values="Weight", index="TF", columns="Gene", fill_value=0)
tf_names = list(motifs.index)
gene_names = list(motifs.columns)

# Create COPD and CTRL phenotype dataframe and miRNA expression dataframe
pheno_df_copd = pheno_df[pheno_df['disease'] == "2-COPD/Emphysema"]
pheno_df_ctrl = pheno_df[pheno_df['disease'] == "3-Control"]
pheno_df_copd.index = pheno_df_copd['title']
pheno_df_ctrl.index = pheno_df_ctrl['title']

input_df_copd = final_expr_full.T.merge(pheno_df_copd, left_index=True, right_index=True)
input_df_copd = input_df_copd[final_expr_full.T.columns]
input_df_ctrl = final_expr_full.T.merge(pheno_df_ctrl, left_index=True, right_index=True)
input_df_ctrl = input_df_ctrl[final_expr_full.T.columns]

# Run PANDA on COPD and CTRL input dataframes and save as feather file
copd_panda_network = my_netzoo.my_panda(input_df_copd.T, motif)
ctrl_panda_network = my_netzoo.my_panda(input_df_ctrl.T, motif)

# Save control and case full panda dataframes
ctrl_panda_network_df = pd.DataFrame(ctrl_panda_network, index=motif["TF"].unique(), columns = input_df_ctrl.columns)
ctrl_panda_network_df.to_feather("../data/ctrl_panda_network.feather")
copd_panda_network_df = pd.DataFrame(copd_panda_network, index=motif["TF"].unique(), columns = input_df_ctrl.columns)
copd_panda_network_df.to_feather("../data/copd_panda_network.feather")

# Just get one main PANDA network
p_network = my_netzoo.my_panda(final_expr_full, motif_df=motif)
p_network_df = pd.DataFrame(p_network, index=tf_names, columns=gene_names)
p_network_df.to_feather("../data/panda_network.feather")

print(final_expr_full.shape)
print(final_expr_full.head())

