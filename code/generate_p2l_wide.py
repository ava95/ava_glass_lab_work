
# This code recapitulates a part of code/pipeline01_generate_ss_networks.py 
# which failed locally after hitting a memory limit.

import pickle
import pandas as pd
from datetime import datetime


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
    
    # Rename columns to TF_Gene
    new_cols = []
    for c in network_df.columns:
        if "_" in c:
            a, b = c.split("_")
            new_cols.append(f"{b}_{a}")
        else:
            new_cols.append(c)
    network_df.columns = new_cols
    
    return(network_df)
    
    
date = datetime.today().strftime('%Y-%m-%d-%H:%M:%S')

p2l_networks = pd.read_pickle("../data/p2l_raw.pkl")

p2l_wide = make_wide(p2l_networks)
p2l_wide.to_feather("../data/" + date + "_p2l_networks_wide.feather")

