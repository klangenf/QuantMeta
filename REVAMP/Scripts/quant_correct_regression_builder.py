#!/usr/bin/env python
"""
Quantification Correction Regression Builder

Purpose: Detecting and correcting quantification errors caused by non-specific 
mapping or assembly errors requires limitations on the acceptable read depth 
variability across target sequences. The acceptable read depth variability may 
differ depending on the library preparation, if PCR amplification was performed, 
and the sequencing technology used as each may introduce different bias and 
increase or decrease how much read depth may be an intrinsic result of sequencing. 
Langenfeld et al. (2025) used Swift 1S Plus library prep for simultaneous sequencing 
of dsDNA and ssDNA with Illumina NovaSeq on SP flowcells to produce 251-bp paired-end 
reads. It is recommended that specific read depth variability thresholds are developed 
for each sequencing protocol.

There are several variables that may need to be manually adjusted in the read depth
variability regressions. Ideal read depth variability regressions will have (1) 
minimize RMSE outliers in each regression and (2) alter regressions so the Q-Q plot 
follows the 1:1 line reasonably closely (will prevent over or under fitting the data).

In Langenfeld et al. (2025), the standards were binned based on the average read depth 
across a target sequence in 0-10, 10-100, 100-1000, and >1000 reads/bp. These bins may 
need to be adjusted depending on the sequencing depth and sample diversity.

The terms in each equation may vary slightly. Based on Browne et al. (2020), Illumina 
technologies introduce as quadratic GC bias in sequencing results. Additionally, the 
total average read depth and E_rel may also significantly impact the observed read 
depth variability. In Langenfeld et al. (2025), E_rel was only a significant factor 
for the least abundant standards with incomplete coverage (i.e., 0-10 reads/bp bin).
"""

import argparse
import pandas as pd
import numpy as np
import json
import pickle
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from scipy.stats import linregress
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

def fit_read_var_regression(X, y, bin_name):
    """
    Fit polynomial regression with appropriate terms for each depth bin.
    
    Parameters:
    -----------
    X : pd.DataFrame
        Feature matrix with columns: avg_GC, total_avg_depth, E_rel (if applicable)
    y : np.array
        Response variable (avg_depth or log-transformed)
    bin_name : str
        Bin identifier ('bin1', 'bin2', 'bin3', or 'bin4')
    
    Returns:
    --------
    dict : Model with coefficients, type, and bin information
    """
    
    if bin_name == 'bin4':  # >= 1000 reads/bp
        # avg_depth ~ poly(avg_GC, 2) + log(gene_copies)
        X_fit = np.column_stack([
            X['avg_GC']**2, X['avg_GC'], np.ones(len(X)),  # poly(avg_GC, 2)
            np.log(X['gene_copies'])
        ])
        response_type = 'linear'
        
    elif bin_name == 'bin3':  # 100-1000 reads/bp
        # log(avg_depth) ~ poly(avg_GC, 2) + log(gene_copies)
        y = np.log(y)
        X_fit = np.column_stack([
            X['avg_GC']**2, X['avg_GC'], np.ones(len(X)),  # poly(avg_GC, 2)
            np.log(X['gene_copies'])
        ])
        response_type = 'log'
        
    elif bin_name == 'bin2':  # 10-100 reads/bp
        # log(avg_depth+1) ~ poly(avg_GC, 2) + poly(log(gene_copies), 2)
        y = np.log(y + 1)
        log_depth = np.log(X['gene_copies'])
        X_fit = np.column_stack([
            X['avg_GC']**2, X['avg_GC'], np.ones(len(X)),  # poly(avg_GC, 2)
            log_depth**2, log_depth  # poly(log(gene_copies), 2)
        ])
        response_type = 'log'
        
    elif bin_name == 'bin1':  # 0-10 reads/bp
        # log(avg_depth+1) ~ poly(avg_GC, 2) + poly(log(gene_copies), 2) + poly(E_rel, 2)
        y = np.log(y + 1)
        log_depth = np.log(X['gene_copies'])
        X_fit = np.column_stack([
            X['avg_GC']**2, X['avg_GC'], np.ones(len(X)),  # poly(avg_GC, 2)
            log_depth**2, log_depth,  # poly(log(gene_copies), 2)
            X['E_rel']**2, X['E_rel']  # poly(E_rel, 2)
        ])
        response_type = 'log'
    
    # Fit OLS regression
    coeffs = np.linalg.lstsq(X_fit, y, rcond=None)[0]
    residuals = y - X_fit @ coeffs
    sse = np.sum(residuals**2)
    mse = sse / (len(y) - len(coeffs))
    rmse = np.sqrt(mse)
    
    return {
        'coefficients': coeffs,
        'response_type': response_type,
        'bin': bin_name,
        'X_test': X,
        'y_original': y,
        'residuals': residuals,
        'rmse': rmse
    }


def predict_read_var(model, X_new):
    """
    Predict using fitted read depth variability model.
    
    Parameters:
    -----------
    model : dict
        Fitted model from fit_read_var_regression
    X_new : pd.DataFrame
        New features with same structure
    
    Returns:
    --------
    np.array : Predictions
    """
    bin_name = model['bin']
    coeffs = model['coefficients']
    
    if bin_name == 'bin4':
        X_fit = np.column_stack([
            X_new['avg_GC']**2, X_new['avg_GC'], np.ones(len(X_new)),
            np.log(X_new['total_avg_depth'])
        ])
        preds = X_fit @ coeffs
        
    elif bin_name == 'bin3':
        X_fit = np.column_stack([
            X_new['avg_GC']**2, X_new['avg_GC'], np.ones(len(X_new)),
            np.log(X_new['total_avg_depth'])
        ])
        preds = X_fit @ coeffs
        preds = np.exp(preds)
        
    elif bin_name == 'bin2':
        log_depth = np.log(X_new['total_avg_depth'])
        X_fit = np.column_stack([
            X_new['avg_GC']**2, X_new['avg_GC'], np.ones(len(X_new)),
            log_depth**2, log_depth
        ])
        preds = X_fit @ coeffs
        preds = np.exp(preds) - 1
        
    elif bin_name == 'bin1':
        log_depth = np.log(X_new['total_avg_depth'])
        X_fit = np.column_stack([
            X_new['avg_GC']**2, X_new['avg_GC'], np.ones(len(X_new)),
            log_depth**2, log_depth,
            X_new['E_rel']**2, X_new['E_rel']
        ])
        preds = X_fit @ coeffs
        preds = np.exp(preds) - 1
    
    return preds

def calculate_rmse(sliding_window, reg1, reg2, reg3, reg4):
    """
    Vectorized RMSE calculation using groupby operations.
    """
    # Add bin assignment
    def assign_bin(depth):
        if depth >= 1000:
            return "≥1,000 reads/bp"
        elif 100 <= depth < 1000:
            return "100-1,000 reads/bp"
        elif 10 <= depth < 100:
            return "10-100 reads/bp"
        elif depth < 10:
            return "0-10 reads/bp"
        else:
            return None
    
    sliding_window = sliding_window.copy()
    sliding_window['bin'] = sliding_window['gene_copies'].apply(assign_bin)
    
    # Group by unique_ID and calculate predictions
    def predict_for_group(group):
        depth = group['gene_copies'].iloc[0]
        bin_label = group['bin'].iloc[0]
        
        # Select appropriate regression
        if depth >= 1000 and reg4 is not None:
            pred = predict_read_var(reg4, group[['avg_GC', 'gene_copies', 'E_rel']])
        elif 100 <= depth < 1000 and reg3 is not None:
            pred = predict_read_var(reg3, group[['avg_GC', 'gene_copies', 'E_rel']])
        elif 10 <= depth < 100 and reg2 is not None:
            pred = predict_read_var(reg2, group[['avg_GC', 'gene_copies', 'E_rel']])
        elif depth < 10 and reg1 is not None:
            pred = predict_read_var(reg1, group[['avg_GC', 'gene_copies', 'E_rel']])
        else:
            return pd.Series(dtype=float)
        
        # Calculate RMSE for this group
        rmse = np.sqrt(np.sum((pred - group['avg_depth'].values)**2) / len(group))
        norm_rmse = rmse / depth
        
        return pd.Series({
            'RMSE': rmse,
            'norm_RMSE': norm_rmse,
            'gene_copies': depth,
            'bin': bin_label
        })
    
    # Apply predictions group-wise
    rmse_results = sliding_window.groupby('unique_ID').apply(predict_for_group).reset_index()
    
    # Add sample info
    sample_info = sliding_window[['unique_ID', 'sample', 'ID']].drop_duplicates('unique_ID')
    rmse_results = rmse_results.merge(sample_info, on='unique_ID', how='left')
    
    # Set categorical bin order
    rmse_results['bin'] = pd.Categorical(
        rmse_results['bin'],
        categories=["0-10 reads/bp", "10-100 reads/bp", "100-1,000 reads/bp", "≥1,000 reads/bp"],
        ordered=True
    )
    
    return rmse_results

def calculate_standardized_residuals(model):
    """
    Calculate standardized residuals from fitted model.
    Approximation of R's rstandard() function.
    
    Parameters:
    -----------
    model : dict
        Fitted model
    
    Returns:
    --------
    np.array : Standardized residuals
    """
    residuals = model['residuals']
    rmse = model['rmse']
    return residuals / rmse


def main():
    parser = argparse.ArgumentParser(description='Quantification Correction Regression Builder.')
    parser.add_argument('--sample-names', required=True, default='Config/sample_list.txt', help='List of samples')
    
    args = parser.parse_args()

    # Create output directories
    output_rdv = 'Regressions/read_depth_variability'
    Path(output_rdv).mkdir(parents=True, exist_ok=True)
    output_trdv = 'Regressions/threshold_read_depth_variability'
    Path(output_trdv).mkdir(parents=True, exist_ok=True)

    # Load standards data
    sample_names = pd.read_csv(args.sample_names, header=None)[0].tolist()

    for _, sample in enumerate(sample_names):
        mapping_path = f'Mapping/{sample}/standards_mapping_analysis.txt'
        window_path = f'Mapping/{sample}/standards_sliding_window.txt'
        if not Path(mapping_path).exists():
            print(f"Warning: Mapping analysis file not found for sample {sample} at {mapping_path}")
            continue
        
        mapping_df = pd.read_csv(mapping_path, sep='\t', header=0)
        mapping_df = pd.DataFrame(mapping_df[['ID', 'E_rel', 'E_detect', 'gene_copies']])
        mapping_df = mapping_df[mapping_df['E_rel'] > mapping_df['E_detect']].copy()  # Filter out standards below detection

        window_df = pd.read_csv(window_path, sep='\t', header=0)
        window_df = pd.DataFrame(window_df[['ID', 'avg_GC', 'avg_depth']])
        window_df = pd.merge(window_df, mapping_df[['ID', 'E_rel', 'gene_copies']], on='ID', how='inner')
        window_df['sample'] = sample

        if _ == 0:
            sliding_window = window_df
        else:
            sliding_window = pd.concat([sliding_window, window_df], ignore_index=True)
    
    
    # Fit regressions for each depth bin
    # Bin 4: >= 1000 reads/bp
    bin4_data = sliding_window[sliding_window['gene_copies'] >= 1000].copy()
    if len(bin4_data) > 5:
        reg4 = fit_read_var_regression(
            bin4_data[['avg_GC', 'gene_copies', 'E_rel']],
            bin4_data['avg_depth'].values,
            'bin4'
        )
    else:
        print("Warning: Insufficient data for bin4 (>= 1000 reads/bp)")
        reg4 = None
    
    # Bin 3: 100-1000 reads/bp
    bin3_data = sliding_window[
        (sliding_window['gene_copies'] >= 100) & 
        (sliding_window['gene_copies'] < 1000)
    ].copy()
    if len(bin3_data) > 5:
        reg3 = fit_read_var_regression(
            bin3_data[['avg_GC', 'gene_copies', 'E_rel']],
            bin3_data['avg_depth'].values,
            'bin3'
        )
    else:
        print("Warning: Insufficient data for bin3 (100-1000 reads/bp)")
        reg3 = None
    
    # Bin 2: 10-100 reads/bp
    bin2_data = sliding_window[
        (sliding_window['gene_copies'] >= 10) & 
        (sliding_window['gene_copies'] < 100)
    ].copy()
    if len(bin2_data) > 5:
        reg2 = fit_read_var_regression(
            bin2_data[['avg_GC', 'gene_copies', 'E_rel']],
            bin2_data['avg_depth'].values,
            'bin2'
        )
    else:
        print("Warning: Insufficient data for bin2 (10-100 reads/bp)")
        reg2 = None
    
    # Bin 1: 0-10 reads/bp
    bin1_data = sliding_window[sliding_window['gene_copies'] < 10].copy()
    if len(bin1_data) > 5:
        reg1 = fit_read_var_regression(
            bin1_data[['avg_GC', 'gene_copies', 'E_rel']],
            bin1_data['avg_depth'].values,
            'bin1'
        )
    else:
        print("Warning: Insufficient data for bin1 (0-10 reads/bp)")
        reg1 = None
    
    # RMSE analysis by unique standard ID
    sliding_window['unique_ID'] = sliding_window['sample'] + "_" + sliding_window['ID'].astype(str)
    
    RMSE_results = calculate_rmse(sliding_window, reg1, reg2, reg3, reg4)
    
    # Boxplot of normalized RMSE
    fig, ax = plt.subplots(figsize=(10, 6))
    RMSE_results.boxplot(column='norm_RMSE', by='bin', ax=ax)
    ax.set_xlabel('Read Depth Bin')
    ax.set_ylabel('Normalized RMSE')
    ax.set_title('Normalized RMSE by Read Depth Bin')
    plt.suptitle('')  # Remove automatic title
    plt.tight_layout()
    plt.savefig(f"{output_rdv}/boxplot_regressions.png", dpi=400, bbox_inches='tight')
    plt.close()
    
    # RMSE scatter plot with thresholds
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.scatter(RMSE_results['gene_copies'], RMSE_results['RMSE'], alpha=0.6")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Average Read Depth (reads/bp)')
    ax.set_ylabel('RMSE')
    ax.legend()
    ax.axvline(x=10, linestyle='--', color='gray', alpha=0.5)
    ax.axvline(x=100, linestyle='--', color='gray', alpha=0.5)
    ax.axvline(x=1000, linestyle='--', color='gray', alpha=0.5)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_trdv}/RMSE_scatter.png", dpi=400, bbox_inches='tight')
    plt.close()
    
    # Q-Q plots for standardized residuals
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    
    if reg1 is not None:
        stdres1 = calculate_standardized_residuals(reg1)
        stats.probplot(stdres1, dist="norm", plot=axes[0])
        axes[0].set_title('0-10 reads/bp')
    
    if reg2 is not None:
        stdres2 = calculate_standardized_residuals(reg2)
        stats.probplot(stdres2, dist="norm", plot=axes[1])
        axes[1].set_title('10-100 reads/bp')
    
    if reg3 is not None:
        stdres3 = calculate_standardized_residuals(reg3)
        stats.probplot(stdres3, dist="norm", plot=axes[2])
        axes[2].set_title('100-1,000 reads/bp')
    
    if reg4 is not None:
        stdres4 = calculate_standardized_residuals(reg4)
        stats.probplot(stdres4, dist="norm", plot=axes[3])
        axes[3].set_title('≥1,000 reads/bp')
    
    for ax in axes:
        ax.set_xlabel('Normal Scores')
        ax.set_ylabel('Standardized Residuals')
    
    plt.tight_layout()
    plt.savefig(f"{output_rdv}/qq_plots_regressions.png", dpi=400, bbox_inches='tight')
    plt.close()
    
    # Save models
    models = {
        'reg1': reg1,
        'reg2': reg2,
        'reg3': reg3,
        'reg4': reg4
    }
    
    for name, model in models.items():
        if model is not None:
            with open(f"{output_rdv}/{name}.pkl", 'wb') as f:
                pickle.dump(model, f)
    
    # Build RMSE threshold functions
    # Aggregate RMSE by rounded depth
    RMSE_limits = RMSE_results.copy()
    RMSE_limits['gene_copies_rounded'] = np.round(RMSE_limits['gene_copies'], 1)
    RMSE_limits = RMSE_limits.groupby('gene_copies_rounded').agg({
        'RMSE': 'max'
    }).reset_index()
    RMSE_limits.columns = ['gene_copies', 'RMSE']
    
    # Fit RMSE limit functions for each bin
    rmse_limit_funcs = {}
    
    # Bin 1: 0-10
    bin1_limits = RMSE_limits[RMSE_limits['gene_copies'] < 10]
    if len(bin1_limits) > 2:
        log_rmse = np.log(bin1_limits['RMSE'].values) + 0.25
        log_depth = np.log(bin1_limits['gene_copies'].values)
        slope, intercept, _, _, _ = linregress(log_depth, log_rmse)
        rmse_limit_funcs['func1'] = {'intercept': intercept, 'slope': slope, 'type': 'log_linear'}
    
    # Bin 2: 10-100
    bin2_limits = RMSE_limits[(RMSE_limits['gene_copies'] >= 10) & 
                               (RMSE_limits['gene_copies'] < 100)]
    if len(bin2_limits) > 2:
        log_rmse = np.log(bin2_limits['RMSE'].values) + 0.25
        log_depth = np.log(bin2_limits['gene_copies'].values)
        slope, intercept, _, _, _ = linregress(log_depth, log_rmse)
        rmse_limit_funcs['func2'] = {'intercept': intercept, 'slope': slope, 'type': 'log_linear'}
    
    # Bin 3: 100-1000
    bin3_limits = RMSE_limits[(RMSE_limits['gene_copies'] >= 100) & 
                               (RMSE_limits['gene_copies'] < 1000)]
    if len(bin3_limits) > 2:
        log_rmse = np.log(bin3_limits['RMSE'].values) + 0.25
        log_depth = np.log(bin3_limits['gene_copies'].values)
        slope, intercept, _, _, _ = linregress(log_depth, log_rmse)
        rmse_limit_funcs['func3'] = {'intercept': intercept, 'slope': slope, 'type': 'log_linear'}
    
    # Bin 4: >= 1000 (use constant)
    bin4_limits = RMSE_limits[RMSE_limits['gene_copies'] >= 1000]
    if len(bin4_limits) > 0:
        rmse_limit_funcs['func4'] = {'value': float(bin4_limits['RMSE'].max() + np.exp(0.25)), 'type': 'constant'}
    
    # Save RMSE limit functions
    for name, func in rmse_limit_funcs.items():
        if func is not None:
            with open(f"{output_trdv}/{name}.pkl", 'wb') as f:
                pickle.dump(func, f)
    
    # Visualize RMSE thresholds
    fig, ax = plt.subplots(figsize=(10, 7))

    ax.scatter(RMSE_results['gene_copies'], RMSE_results['RMSE'], alpha=0.6)
    
    # Plot threshold lines
    depths_log = np.logspace(0, 4, 100)
    
    # Bin 1
    if 'func1' in rmse_limit_funcs:
        func = rmse_limit_funcs['func1']
        preds = np.exp(func['intercept'] + func['slope'] * np.log(np.minimum(depths_log, 10)))
        ax.plot(depths_log[depths_log < 10], preds[depths_log < 10], 'r-', linewidth=2)
    
    # Bin 2
    if 'func2' in rmse_limit_funcs:
        func = rmse_limit_funcs['func2']
        mask = (depths_log >= 10) & (depths_log < 100)
        preds = np.exp(func['intercept'] + func['slope'] * np.log(depths_log[mask]))
        ax.plot(depths_log[mask], preds, 'r-', linewidth=2)
    
    # Bin 3
    if 'func3' in rmse_limit_funcs:
        func = rmse_limit_funcs['func3']
        mask = (depths_log >= 100) & (depths_log < 1000)
        preds = np.exp(func['intercept'] + func['slope'] * np.log(depths_log[mask]))
        ax.plot(depths_log[mask], preds, 'r-', linewidth=2)
    
    # Bin 4
    if 'func4' in rmse_limit_funcs:
        func = rmse_limit_funcs['func4']
        ax.axhline(y=func['value'], xmin=0.75, color='r', linestyle='-', linewidth=2)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Average Read Depth (reads/bp)')
    ax.set_ylabel('RMSE')
    ax.axvline(x=10, linestyle='--', color='gray', alpha=0.5)
    ax.axvline(x=100, linestyle='--', color='gray', alpha=0.5)
    ax.axvline(x=1000, linestyle='--', color='gray', alpha=0.5)
    ax.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_trdv}/RMSE_thresholds.png", dpi=400, bbox_inches='tight')
    plt.close()
    
    return {
        'regressions': models,
        'RMSE_results': RMSE_results,
        'RMSE_limits': RMSE_limits,
        'RMSE_limit_functions': rmse_limit_funcs
    }


if __name__ == '__main__':
    main()
