import pandas as pd
from datetime import datetime

# date = datetime.today().strftime('%Y-%m-%d-%H:%M:%S')

# p2l_wide = pd.read_feather("../data/2025-12-04-15:23:29_p2l_networks_wide.feather")
# pheno_df = pd.read_csv("../data/pheno_df.csv", index_col=0)


# l2p_merged = l2p_wide.merge(pheno_df, left_on='sample', right_on='title', how='left')
# l2p_merged.to_feather("../data/" + date + "_l2p_merged.feather")

# p2l_merged = p2l_wide.merge(pheno_df, left_on='sample', right_on='title', how='left')
# p2l_merged.to_feather("../data/" + date + "_p2l_merged.feather")



# new = pd.read_feather("../data/2025-12-04-18:54:47_p2l_merged.feather")
# old = pd.read_feather("../data/p2l_merged.feather")

# print(new["percent.emphysema"])
# print(old.head())
