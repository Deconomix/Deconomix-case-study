# Libraries
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import deconomix as dcx
from scipy.stats import spearmanr, pearsonr

# Functions
def create_composition_overview(data_train, data_test):
    """
    Creates an overview table of cell type counts and percentages for train and test datasets.
    
    Parameters
    ----------
    data_train : pd.DataFrame
        Dataframe where columns are cell types for the 'train' set.
    data_test : pd.DataFrame
        Dataframe where columns are cell types for the 'test' set.
        
    Returns
    -------
    pd.DataFrame
        Summary table with columns: 'cell type', 'train number', 'train percentage', 
        'test number', 'test percentage', and additional summary rows.
    """

    def summarize_cell_types(data, count_col, pct_col):
        cell_type_counts = pd.Series(data.columns).value_counts().sort_index()
        relative_abundance = cell_type_counts / cell_type_counts.sum()
        summary_table = pd.DataFrame({
            'cell type': cell_type_counts.index,
            count_col: cell_type_counts.values,
            pct_col: (relative_abundance.values * 100).round(3)
        })
        return summary_table

    table_train = summarize_cell_types(data_train, 'train number', 'train percentage')
    table_test = summarize_cell_types(data_test, 'test number', 'test percentage')

    merged_table = pd.merge(
        table_train,
        table_test,
        on='cell type',
        how='outer'
    )

    # Add a row "sum immune cell types" which sums all columns except for 'hidden'
    immune_rows = merged_table[merged_table['cell type'] != 'hidden']
    immune_sum = pd.Series({
        'cell type': 'sum immune cell types',
        'train number': immune_rows['train number'].dropna().sum(),
        'train percentage': immune_rows['train percentage'].dropna().sum(),
        'test number': immune_rows['test number'].dropna().sum(),
        'test percentage': immune_rows['test percentage'].dropna().sum()
    })
    merged_table = pd.concat([merged_table, pd.DataFrame([immune_sum])], ignore_index=True)

    # Add a row "total" which sums up the rows 'hidden' and 'sum immune cell types'
    total_rows = merged_table[merged_table['cell type'].isin(['hidden', 'sum immune cell types'])]

    total_sum = pd.Series({
        'cell type': 'total',
        'train number': total_rows['train number'].dropna().sum(),
        'train percentage': total_rows['train percentage'].dropna().sum(),
        'test number': total_rows['test number'].dropna().sum(),
        'test percentage': total_rows['test percentage'].dropna().sum()
    })

    merged_table = pd.concat([merged_table, pd.DataFrame([total_sum])], ignore_index=True)

    # Set numeric columns type
    merged_table['train number'] = merged_table['train number'].astype('Int64')
    merged_table['test number'] = merged_table['test number'].astype('Int64')

    # Rename the 'hidden' row to 'hidden contribution'
    merged_table['cell type'] = merged_table['cell type'].replace({'hidden': 'hidden contribution'})

    # Get the indices for relevant rows
    idx_sum = merged_table.index[merged_table['cell type'] == 'sum immune cell types'].tolist()
    idx_hidden = merged_table.index[merged_table['cell type'] == 'hidden contribution'].tolist()
    idx_total = merged_table.index[merged_table['cell type'] == 'total'].tolist()

    # Drop these rows temporarily
    main_table = merged_table.drop(idx_sum + idx_hidden + idx_total)

    # Concatenate in the new order: (main rows) + sum + hidden + total
    ordered_table = pd.concat([
        main_table,
        merged_table.loc[idx_sum],
        merged_table.loc[idx_hidden],
        merged_table.loc[idx_total]
    ], ignore_index=True)

    return ordered_table


def create_performance_overview(data_train, data_test):
    # Simulate train and test bulks
    X_ref, Y_train, C_train = dcx.utils.simulate_data(data_train, n_mixtures=20000, n_cells_in_mix=100)
    X_test, Y_test, C_test = dcx.utils.simulate_data(data_test, n_mixtures=10000, n_cells_in_mix=100)

    # Split Y_test and C_test in 10 parts of 1000 samples each
    Y_test_splits = [Y_test.iloc[:,i*1000:(i+1)*1000] for i in range(10)]
    C_test_splits = [C_test.iloc[:,i*1000:(i+1)*1000] for i in range(10)]

    # DTD
    model_DTD = dcx.methods.DTD(X_ref, Y_train, C_train)
    model_DTD.run(iterations=1000)
    gamma_naive = pd.Series(np.ones(5000))
    gamma_naive.index = model_DTD.gamma.index

    # Training performance
    C_est_train = dcx.utils.calculate_estimated_composition(X_ref, Y_train, model_DTD.gamma)
    C_est_train_naive = dcx.utils.calculate_estimated_composition(X_ref, Y_train, gamma_naive)
    corr_train_naive = dcx.utils.calculate_corr(C_train, C_est_train_naive)
    corr_train = dcx.utils.calculate_corr(C_train, C_est_train)

    # Test performance
    corr_test = []
    corr_test_naive = []
    corr_test_ADTD = []
    corr_x_ADTD = []
    corr_x_ADTD_pearson = []
    x_true = X_test.loc[:, 'hidden'].values.flatten()

    for Y_test_split, C_test_split in zip(Y_test_splits, C_test_splits):
        # DTD inference for test split
        C_est_test = dcx.utils.calculate_estimated_composition(X_ref, Y_test_split, model_DTD.gamma)
        C_est_test_naive = dcx.utils.calculate_estimated_composition(X_ref, Y_test_split, gamma_naive)
        corr_test.append(dcx.utils.calculate_corr(C_test_split, C_est_test))
        corr_test_naive.append(dcx.utils.calculate_corr(C_test_split, C_est_test_naive))
        # ADTD model
        model_adtd = dcx.methods.ADTD(X_ref, Y_test_split, model_DTD.gamma, C_static=False, Delta_static=True)
        model_adtd.run()
        corr_adtd = dcx.utils.calculate_corr(C_test_split, model_adtd.C_est, hidden_ct='hidden', c_est=model_adtd.c_est)
        corr_test_ADTD.append(corr_adtd)
        x_corr, _ = spearmanr(x_true, model_adtd.x_est.values.flatten())
        corr_x_ADTD.append(x_corr)
        x_corr_pearson, _ = pearsonr(x_true, model_adtd.x_est.values.flatten())
        corr_x_ADTD_pearson.append(x_corr_pearson)

    # Post-processing for Performance Summary Table
    # Combine the series in corr_test into a DataFrame, with split number as column names
    corr_test_df = pd.concat(corr_test, axis=1)
    corr_test_df.columns = [f"split {i+1}" for i in range(len(corr_test))]

    # Do the same for the naive results
    corr_test_naive_df = pd.concat(corr_test_naive, axis=1)
    corr_test_naive_df.columns = [f"split {i+1}" for i in range(len(corr_test_naive))]

    # Only use the splits, not the 'Average' column:
    split_cols = [col for col in corr_test_df.columns if col.startswith("split ")]
    corr_test_df["AVG"] = corr_test_df[split_cols].mean(axis=1)
    corr_test_df["STD"] = corr_test_df[split_cols].std(axis=1)
    corr_test_df = corr_test_df.round(3)

    split_cols_naive = [col for col in corr_test_naive_df.columns if col.startswith("split ")]
    corr_test_naive_df["AVG"] = corr_test_naive_df[split_cols_naive].mean(axis=1)
    corr_test_naive_df["STD"] = corr_test_naive_df[split_cols_naive].std(axis=1)
    corr_test_naive_df = corr_test_naive_df.round(3)

    # ADTD performance consensus profile x
    corr_x_ADTD_pearson_AVG = np.array(corr_x_ADTD_pearson).mean()
    corr_x_ADTD_pearson_STD = np.array(corr_x_ADTD_pearson).std()
    corr_x_ADTD_spearman_AVG = np.array(corr_x_ADTD).mean()
    corr_x_ADTD_spearman_STD = np.array(corr_x_ADTD).std()

    corr_ADTD_df = pd.concat(corr_test_ADTD, axis=1)
    corr_ADTD_df.columns = [f"split {i+1}" for i in range(len(corr_test_ADTD))]
    split_cols = [col for col in corr_ADTD_df.columns if col.startswith("split ")]
    corr_ADTD_df["AVG"] = corr_ADTD_df[split_cols].mean(axis=1)
    corr_ADTD_df["STD"] = corr_ADTD_df[split_cols].std(axis=1)
    corr_ADTD_df = corr_ADTD_df.round(3)

    summary_df = pd.DataFrame({
        'Train Naive': corr_train_naive.round(3),
        'Train DTD': corr_train.round(3),
        'Test Naive': corr_test_naive_df['AVG'].astype(str) + " ± " + corr_test_naive_df['STD'].astype(str),
        'Test DTD': corr_test_df['AVG'].astype(str) + " ± " + corr_test_df['STD'].astype(str),
        'Test ADTD': corr_ADTD_df['AVG'].astype(str) + " ± " + corr_ADTD_df['STD'].astype(str)
    })

    train_naive_avg = corr_train_naive.mean().round(3)
    train_DTD_avg = corr_train.mean().round(3)
    test_naive_avg = corr_test_naive_df['AVG'].mean().round(3)
    test_naive_std = corr_test_naive_df['STD'].mean().round(3)
    test_DTD_avg = corr_test_df['AVG'].mean().round(3)
    test_DTD_std = corr_test_df['STD'].mean().round(3)
    test_ADTD_avg = corr_ADTD_df.loc[corr_ADTD_df.index != 'hidden', 'AVG'].mean().round(3)
    test_ADTD_std = corr_ADTD_df.loc[corr_ADTD_df.index != 'hidden', 'STD'].mean().round(3)

    # Append new row
    summary_df.loc['mean (excl. hidden)'] = [
        train_naive_avg,
        train_DTD_avg,
        str(test_naive_avg) + " ± " + str(test_naive_std),
        str(test_DTD_avg) + " ± " + str(test_DTD_std),
        str(test_ADTD_avg) + " ± " + str(test_ADTD_std)
    ]

    # Swap the last two rows of the dataframe if needed
    if list(summary_df.index[-2:]) == ['hidden', 'mean (excl. hidden)']:
        summary_df = pd.concat([summary_df.iloc[:-2], summary_df.iloc[[-1, -2]]])
        summary_df.index = list(summary_df.index[:-2]) + [summary_df.index[-2], summary_df.index[-1]]

    summary_df.loc['hidden reference', 'Test ADTD'] = f"{corr_x_ADTD_spearman_AVG.round(3).item()} ± {corr_x_ADTD_spearman_STD.round(3).item()}"

    summary_df = summary_df.reset_index().rename(columns={'index': 'cell type'})

    return summary_df