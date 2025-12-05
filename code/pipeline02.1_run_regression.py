import pandas as pd
from scipy.stats import pearsonr
from scipy import stats
import my_netzoo
from datetime import datetime
import numpy as np

# Analyze correlation between edges and fev1/fvc
def fev_correlation_test(merged_network):
    """
    For each edge, calculate the correlation between the edge score and the fev score

    Args:
        merged_network (dataframe): samples x (edges + pheno) dataframe with column fev1.fvc.ratio

    Returns:
        results (dataframe): dataframe with row = edge, correlation and pvalue
    """

    # Drop any samples with an NA value for fev1/fvc
    merged_network = merged_network.dropna(subset=['fev1.fvc.ratio'])
    
    # Columns to exclude from correlation (these may or may not be present)
    exclude_cols = ['title', 'fev1.fvc.ratio', 'gender', 'age', 'disease', 'smoking.status', 
                    'race', 'percent.emphysema', 'pack.years', 'sample', 'disease_group'] 
    edge_cols = [col for col in merged_network.columns if col not in exclude_cols]
    
    x = merged_network["fev1.fvc.ratio"]
    
    # Correlate per edge
    results = []
    for col in edge_cols:
        print(col)
        y = merged_network[col].values
        c, p = pearsonr(x, y)
        results.append({'edge': col, 'correlation': c, 'p_value': p})
    return pd.DataFrame(results)


# Conduct ttest COPD L to P
def copd_ttest(merged_network):
    """_summary_

    Args:
        merged_network (DataFrame): samples x (edges + pheno) dataframe with column disease_group

    Returns:
        results (DataFrame): rows = edge, values for mean difference and pvalue
    """
    
    # Drop sample rows with no disease value
    merged_network = merged_network.dropna(subset=['disease_group'])
    
    # Exclude columns
    exclude_cols = ['title', 'fev1.fvc.ratio', 'gender', 'age', 'disease', 'smoking.status', 
                    'race', 'percent.emphysema', 'pack.years', 'sample', 'disease_group'] 
    edge_cols = [col for col in merged_network.columns if col not in exclude_cols]
    
    results = []
    for edge in edge_cols:
        print(edge)
        group0 = merged_network.loc[merged_network['disease_group'] == 0, edge].dropna()
        group1 = merged_network.loc[merged_network['disease_group'] == 1, edge].dropna()
        t_stat,p_value = stats.ttest_ind(group0, group1)
        mean_diff = group0.mean() - group1.mean()
        results.append({'edge': edge, 'mean_diff': mean_diff, 'p_value': p_value})
    return pd.DataFrame(results)


# Conduct linear regression of each singular edge
def edgewise_multiple_regression(
        merged_network,
        covariates=['gender', 'age'],
        outcome_col='fev1.fvc.ratio'):
    """
    Run a linear regression on each edge on specified covariates and otucome

    Args:
        merged_network (DataFrame): DataFrame with one row per sample, columns = [sample, covariates, edges...]
        covariates (list, optional): Covariates to incorporate in linear regression. Defaults to ['gender', 'age'].
        outcome_col (str, optional): Outcome value of regression model. Defaults to 'fev1.fvc.ratio'.

    Returns:
        results (DataFrame): Row = edge, beta score and p-value of edge in relation to outcome
    """

    merged_network = merged_network.copy()
    
    # Encode gender if needed
    if 'gender' in covariates and merged_network['gender'].dtype == object:
        merged_network['gender'] = merged_network['gender'].map({'Male': 0, 'Female': 1})
        
    if 'race' in covariates:
        # Binary encoding: 0 = White, 1 = Non-White
        def encode_race(x):
            if pd.isna(x):
                return np.nan
            elif '1-White (Caucasian)' in str(x):
                return 0
            else:
                return 1  # All other races
    
        merged_network['race'] = merged_network['race'].apply(encode_race)
        merged_network['race'] = pd.to_numeric(merged_network['race'], errors='coerce')
        
    if 'pack.years' in covariates:
        merged_network['pack.years'] = pd.to_numeric(merged_network['pack.years'], errors='coerce')
        print(merged_network['pack.years'])
        
    # Only include from pheno_df edge columns that aren't in the exclude_cols list
    exclude_cols = ['geo_accession', 'source_name_ch1', 'organism_ch1', 
                    'molecule_ch1', 'extract_protocol_ch1', 'extract_protocol_ch1.1', 
                    'taxid_ch1', 'data_processing', 'data_processing.1', 
                    'data_processing.2', 'data_processing.3', 'data_processing.4', 
                    'platform_id', 'contact_zip.postal_code', 'instrument_model', 
                    'library_selection', 'library_source', 'library_strategy', 
                    'supplementary_file_1']
    
    edge_cols = [
    col for col in merged_network.columns
    if "_" in col and col not in exclude_cols]

    results = []
    count = 1
    
    
    for edge in edge_cols:
        if count % 10000 == 0:
            print(count) 
        
        subset_cols = [edge] + covariates + [outcome_col]
        df = merged_network[subset_cols].copy()
        
        for col in subset_cols:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        df = df.dropna()
        
        if len(df) < len(covariates) + 2:  # Need enough samples
            continue
        
        y = df[outcome_col].values
        
        if df[edge].std() < 1e-8:
            continue
        
        # Build design matrix
        X = df[[edge] + covariates].values
        
        # **ENSURE X IS FLOAT**
        X = X.astype(np.float64)
        
        # Add intercept column
        X = np.column_stack((np.ones(X.shape[0]), X))

        # Run bare OLS regression
        beta, se_beta, t_stats, p_values = my_netzoo.bare_ols_with_pvalue(X, y)
        results.append({
            'edge': edge,
            'beta_score': beta[1],
            'p_value': p_values[1]
        })
        
        count+=1
        
    return pd.DataFrame(results)


def main():
    
    date = datetime.today().strftime('%Y-%m-%d-%H:%M:%S')

    # Read in data
    l2p_merged = pd.read_feather("../data/2025-12-04-16:39:25_l2p_merged.feather")
    p2l_merged = pd.read_feather("../data/2025-12-04-18:54:47_p2l_merged.feather")
    new_pheno_df = pd.read_csv("../data/pheno_df.csv")
    
    print("Data successfully loaded...")

    disease_map = {'3-Control': 0, '2-COPD/Emphysema': 1}

    l2p_merged['diseasegroup'] = l2p_merged['disease'].map(disease_map)
    p2l_merged['diseasegroup'] = l2p_merged['disease'].map(disease_map)
    
    
    # Create a mapping dictionary from pheno_df
    pack_years_map = new_pheno_df.set_index('title')['pack.years'].to_dict()

    # Update l2p_merged using the mapping
    l2p_merged['pack.years'] = l2p_merged['sample'].map(pack_years_map)
    p2l_merged['pack.years'] = p2l_merged['sample'].map(pack_years_map)
    
    
    print("Calculating P2L regression")
    p2l_regression = edgewise_multiple_regression(p2l_merged, covariates=['gender', 'age', 'pack.years', 'race'])
    p2l_regression.to_feather("../data/" + date + "_p2l_regression_gender_age_pckyr_race.feather")
    
    print("Calculating L2P regresssion")
    l2p_regression = edgewise_multiple_regression(l2p_merged, covariates=['gender', 'age', 'pack.years', 'race'])
    l2p_regression.to_feather("../data/" + date + "_l2p_regression_gender_age_pckyr_race.feather")

    # plt.hist(l2p_regression.p_value, bins=20)
    # plt.hist(p2l_regression.p_value, bins=20)

    # l2p_fev_test = fev_test(l2p_merged)
    # p2l_fev_test = fev_test(p2l_merged)
    # plt.hist(p2l_fev_test.p_value, bins=20)
    # l2p_copd_ttest = copd_ttest(l2p_merged)
    # p2l_copd_ttest = copd_ttest(p2l_networks)
    # plt.hist(l2p_copd_ttest.p_value, bins=20)

    # regression.to_feather('../data/p2l_regression.feather')
    
main()