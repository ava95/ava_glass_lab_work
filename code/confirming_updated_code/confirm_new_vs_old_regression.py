import pandas as pd

new_regression_p2l = pd.read_feather("../data/2025-12-04-19:27:56_p2l_regression_gender_and_age.feather")
new_regression_l2p = pd.read_feather("../data/2025-12-04-19:27:56_l2p_regression_gender_and_age.feather")

old_regression_p2l = pd.read_feather("../data/p2l_gender_and_age_regression.feather")
old_regression_l2p = pd.read_feather("../data/l2p_gender_and_age_regression.feather")

print(new_regression_l2p.head())
print(old_regression_l2p.head())
